"""Conditional coordinates and projection onto the frozen reference tree."""
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import logsumexp, softmax
from scipy.stats import norm


def reference_transform(values, detected, positive, nzero, nref):
    """Apply the existing detection-aware Ctr midrank transform to each group."""
    z = np.empty_like(values, dtype=float)
    low = np.zeros_like(detected, dtype=bool)
    high = np.zeros_like(detected, dtype=bool)
    for j, ref in enumerate(positive):
        seen = detected[:, j]
        u = np.full(len(values), max(nzero[j] / (2 * nref), .5 / nref))
        x = values[seen, j]
        left = np.searchsorted(ref, x, side='left')
        right = np.searchsorted(ref, x, side='right')
        ranks = np.clip((left + right) / 2, .5, len(ref) - .5)
        u[seen] = (nzero[j] + ranks) / nref
        low[seen, j] = x < ref[0]
        high[seen, j] = x > ref[-1]
        z[:, j] = norm.ppf(np.clip(u, .5 / nref, 1 - .5 / nref))
    return z, low, high


def reduced_objective(z, v, gram, nodes, gamma, sigma):
    """Profile r out; omit the constant x'x when comparing candidate positions."""
    dist2 = np.sum((z[..., None, :] - nodes) ** 2, axis=-1)
    return (np.einsum('...i,ij,...j->...', z, gram, z)
            - 2 * np.sum(z * v, axis=-1)
            - gamma * sigma * logsumexp(-dist2 / sigma, axis=-1))


def value_grad(z, v, gram, nodes, gamma, sigma):
    dist2 = np.sum((nodes - z) ** 2, axis=1)
    r = softmax(-dist2 / sigma)
    f = z @ gram @ z - 2 * z @ v - gamma * sigma * logsumexp(-dist2 / sigma)
    grad = 2 * (gram @ z - v + gamma * (z - r @ nodes))
    return f, grad


def solve_coordinates(x, w, nodes, edge_indices, gamma, sigma):
    """Multi-start alternating minimization, followed by gradient checks."""
    gram = w.T @ w
    inv_a = np.linalg.inv(gram + gamma * np.eye(2))
    v = x @ w
    mapped_nodes = (v[:, None, :] + gamma * nodes[None, :, :]) @ inv_a
    edge_mid = (nodes[edge_indices[:, 0]] + nodes[edge_indices[:, 1]]) / 2
    linear = v @ np.linalg.inv(gram)
    seeds = np.concatenate([
        mapped_nodes, np.broadcast_to(edge_mid, (len(x), len(edge_mid), 2)),
        linear[:, None, :]], axis=1)
    initial_values = reduced_objective(seeds, v[:, None, :], gram, nodes, gamma, sigma)
    z = seeds.copy()
    # Coordinate descent is monotone in the full conditional objective.
    for iteration in range(1000):
        r = softmax(-np.sum((z[:, :, None, :] - nodes) ** 2, axis=-1) / sigma, axis=-1)
        new_z = (v[:, None, :] + gamma * (r @ nodes)) @ inv_a
        step = np.linalg.norm(new_z - z, axis=-1)
        z = new_z
        if step.max() < 1e-11:
            break
    # Only repair starts that have not settled, using the same exact objective.
    repairs = np.argwhere(step >= 1e-8)
    for i, j in repairs:
        opt = minimize(value_grad, z[i, j], args=(v[i], gram, nodes, gamma, sigma), jac=True,
                       method='BFGS', options={'gtol': 1e-9, 'maxiter': 1000})
        z[i, j] = opt.x
    values = reduced_objective(z, v[:, None, :], gram, nodes, gamma, sigma)
    assert np.max(values - initial_values) < 1e-7
    best = values.argmin(axis=1)
    row = np.arange(len(x))
    winner = z[row, best]
    f = values[row, best]
    far = np.linalg.norm(z - winner[:, None, :], axis=-1) > 1e-5
    distinct_values = np.where(far, values, np.inf)
    second = distinct_values.argmin(axis=1)
    second_value = distinct_values[row, second]
    alternative = z[row, second].copy()
    alternative[~np.isfinite(second_value)] = np.nan
    r = softmax(-np.sum((winner[:, None, :] - nodes) ** 2, axis=-1) / sigma, axis=-1)
    mean = r @ nodes
    grad = 2 * (winner @ gram - v + gamma * (winner - mean))
    centered = nodes[None, :, :] - mean[:, None, :]
    covariance = np.einsum('nk,nki,nkj->nij', r, centered, centered)
    hessian = 2 * (gram + gamma * np.eye(2))[None] - 4 * gamma / sigma * covariance
    hessmin = np.linalg.eigvalsh(hessian)[:, 0]
    assert np.max(np.linalg.norm(grad, axis=1)) < 1e-5
    assert np.min(hessmin) > -1e-5
    # Exact bound from log-sum-exp, reported without claiming global optimality.
    node_hard_min = np.min(
        np.einsum('nki,ij,nkj->nk', mapped_nodes, gram, mapped_nodes)
        - 2 * np.sum(mapped_nodes * v[:, None, :], axis=-1)
        + gamma * np.sum((mapped_nodes - nodes[None]) ** 2, axis=-1), axis=1)
    lower_bound = node_hard_min - gamma * sigma * np.log(len(nodes))
    assert np.all(f >= lower_bound - 1e-7)
    return winner, pd.DataFrame({
        'objective_without_x2': f,
        'gradient_norm': np.linalg.norm(grad, axis=1),
        'hessian_min_eigenvalue': hessmin,
        'alternate_objective_gap': np.where(np.isfinite(second_value), second_value - f, np.nan),
        'alternative_Z1': alternative[:, 0], 'alternative_Z2': alternative[:, 1],
        'objective_lower_bound_gap': f - lower_bound,
        'max_node_weight': r.max(axis=1),
        'multistart_count': z.shape[1],
        'max_alternating_iterations': iteration + 1,
        'starts_repaired_with_bfgs': len(repairs),
    })


def freeze_tree(nodes_frame, edges, root_xy, native_projection):
    coordinates={str(k):v for k,v in zip(nodes_frame.index,nodes_frame.to_numpy(float))}
    adjacency={k:[] for k in coordinates}
    for a,b in edges[['from','to']].itertuples(index=False,name=None):
        adjacency[a].append(b);adjacency[b].append(a)
    assert len(edges)==len(coordinates)-1
    leaves=sorted(k for k,v in adjacency.items() if len(v)==1)
    original_branches=sorted(k for k,v in adjacency.items() if len(v)>2)
    root_leaf=min(leaves,key=lambda k:np.linalg.norm(coordinates[k]-root_xy))
    extensions=[]
    root_node=root_leaf
    for leaf in leaves:
        neighbor=adjacency[leaf][0]
        outward=coordinates[leaf]-coordinates[neighbor]
        outward=outward/np.linalg.norm(outward)
        native=native_projection.loc[native_projection.nearest_center.eq(leaf),['P1','P2']].to_numpy()
        extension=max(0.,float(np.max((native-coordinates[leaf])@outward))) if len(native) else 0.
        bound=''
        if extension>1e-12:
            bound=f'Ctr_bound_{leaf}'
            coordinates[bound]=coordinates[leaf]+extension*outward
            adjacency[bound]=[leaf];adjacency[leaf].append(bound)
            if leaf==root_leaf:root_node=bound
        extensions.append(dict(learned_leaf=leaf,bound_node=bound,extension_length=extension,native_cells=len(native),root_side=leaf==root_leaf))
    # Each segment ends at a fork or a leaf. Numbering is deterministic from the root.
    seen=set();distance={root_node:0.};parent={root_node:None};rows=[];segment_count=0
    def walk(node,segment):
        nonlocal segment_count
        assert node not in seen
        seen.add(node)
        children=sorted(k for k in adjacency[node] if k!=parent[node])
        for child in children:
            next_segment=segment
            if segment is None or len(children)>1:
                segment_count+=1;next_segment=f'S{segment_count}'
            parent[child]=node
            length=float(np.linalg.norm(coordinates[child]-coordinates[node]))
            assert length>0
            distance[child]=distance[node]+length
            rows.append(dict(graph_edge=len(rows),from_node=node,to_node=child,
                x=coordinates[node][0],y=coordinates[node][1],xend=coordinates[child][0],yend=coordinates[child][1],
                length=length,root_distance_from=distance[node],root_distance_to=distance[child],tree_segment=next_segment,
                terminal_edge=len(adjacency[child])==1))
            walk(child,next_segment)
    walk(root_node,None)
    assert seen==set(coordinates)
    edge_frame=pd.DataFrame(rows)
    max_depth=max(distance.values())
    node_frame=pd.DataFrame([dict(node_id=k,x=v[0],y=v[1],root_distance=distance[k],t_graph=distance[k]/max_depth,
        degree=len(adjacency[k]),root_node=k==root_node,learned_node=k in nodes_frame.index) for k,v in coordinates.items()])
    return dict(nodes=node_frame,edges=edge_frame,extensions=pd.DataFrame(extensions),
        root_node=root_node,root_leaf=root_leaf,learned_branch_nodes=original_branches,
        learned_leaf_nodes=leaves,max_root_distance=max_depth,total_tree_length=edge_frame.length.sum())


def project_graph(z,tree):
    edges=tree['edges'];a=edges[['x','y']].to_numpy();b=edges[['xend','yend']].to_numpy()
    delta=b-a;lengths=edges.length.to_numpy()
    fraction=np.clip(np.einsum('nej,ej->ne',z[:,None,:]-a,delta)/lengths[None,:]**2,0,1)
    q=a[None]+fraction[...,None]*delta[None]
    dist2=np.sum((z[:,None]-q)**2,axis=2)
    edge=dist2.argmin(axis=1);row=np.arange(len(z));selected=q[row,edge];f=fraction[row,edge]
    distance=edges.root_distance_from.to_numpy()[edge]+f*lengths[edge]
    terminal=(f>=1-1e-10)&edges.terminal_edge.to_numpy()[edge]
    return pd.DataFrame(dict(graph_x=selected[:,0],graph_y=selected[:,1],graph_edge=edge,
        graph_from=edges.from_node.to_numpy()[edge],graph_to=edges.to_node.to_numpy()[edge],
        tree_segment=edges.tree_segment.to_numpy()[edge],graph_edge_fraction=f,
        root_distance=distance,t_graph=np.clip(distance/tree['max_root_distance'],0,1),
        distance_to_tree=np.sqrt(dist2[row,edge]),at_early_endpoint=distance<=1e-10,
        at_late_endpoint=terminal,terminal_endpoint=np.where(terminal,edges.to_node.to_numpy()[edge],'')))
