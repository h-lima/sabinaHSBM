import_modules <- function(){
return('
import numpy as np
import pandas as pd
import pickle
import itertools
import re
import os
import sys
import gc
from scipy.special import logsumexp
from graph_tool.all import *')

}

add_names_vertex_prop <- function(){

return('
def add_names_vertex_prop(g):
    r_names = []
    for r in g.ep.v1_names:
        if r not in r_names:
            r_names.append(r)
    n_rows = len(r_names)

    vnames =  g.new_vp("string")
    vbipartite =  g.new_vp("int")
    for v in g.vertices():
        for e in v.all_edges():
            edge = g.edge(e.source(), e.target())
            if v <= (n_rows - 1):
                vbipartite[v] = 0
                vnames[v] = g.ep.v1_names[edge]
            else:
                vbipartite[v] = 1
                vnames[v] = g.ep.v2_names[edge]
            break

    g.vertex_properties["names"] = vnames
    g.vertex_properties["bipartite"] = vbipartite')
}

create_graph <- function(){
return('
def create_graph(df,
    eprop_types = ["double", "string", "string", "string","int", "int"],
    add_bipartite = True):
    g = Graph(directed = False)
    eprops = [g.new_ep(t) for t in eprop_types]
    g.add_edge_list(list(df.to_records(index = False)), hashed = False,
               eprops = eprops)
    eprop_names = ["value", "v1_names", "v2_names", "edge_type", "n", "x"]
    for i, p in enumerate(eprops):
        ename = eprop_names[i]
        g.ep[ename] = p

    if add_bipartite:
        add_names_vertex_prop(g)

    return(g)')
}



get_missing_edges <- function(){
return('
def get_missing_edges(df, g):

    edge_types = g.ep.edge_type.get_2d_array(range(1))[0]
    is_doc = [e_doc == "documented" for e_doc in edge_types]
    g_doc = GraphView(g, efilt = is_doc)

    v1s = df.v1.unique()
    v2s = df.v2.unique()

    all_edges = itertools.product(v1s, v2s)

    all_missing = [e for e in all_edges if g_doc.edge(e[0], e[1]) is None]

    return(all_missing)')
}

hsbm_link_prediction <- function(){

return('
def hsbm_predict(g, elist, wait = 1000,
                 niter = 10, force_niter = 10000, rnd_seed = None,
                 n_default = 1, x_default = 0, alpha = 1, beta = 1,
                 all_missing = None,
                 method = "conditional_missing", save_pickle = False):

    if rnd_seed:
        rnd_seed = int(rnd_seed)
        np.random.seed(rnd_seed)
        seed_rng(rnd_seed)

    N = g.num_vertices()
    E = g.num_edges()

    # initial hierarchical partition
    bs = [g.get_vertices()] + [np.zeros(1)] * 5

    state = MeasuredBlockState(g, n = g.ep.n, x = g.ep.x,
                               fn_params = dict(alpha = alpha,
                                                beta = beta),
                               n_default = n_default,
                               x_default = x_default,
                               state_args = dict(bs = bs))

    mcmc_equilibrate(state, wait = wait,
                     mcmc_args = dict(niter = niter))

    # Marginal posterior edge probabilities
    entropy = state.entropy()
    state_min_dl = state.copy()

    if method == "conditional_missing":
        if all_missing == None:
            all_missing = get_missing_edges(df = elist, g = g)

        edges_prob = []
        def collect_marginals(s):
            nonlocal entropy, state_min_dl, edges_prob
            edges_prob.append(s.get_edges_prob(all_missing))
            if s.entropy() < entropy:
                entropy = s.entropy()
                state_min_dl = s.copy()
        print("Using conditional_missing...")
    elif method == "marginal_all":
        print("Computing marginal_all...")
        # Collecting marginal for network reconstruction
        u = None
        def collect_marginals(s):
            nonlocal u, entropy, state_min_dl
            u = s.collect_marginal(u)
            # Save min description length state
            if s.entropy() < entropy:
                entropy = s.entropy()
                state_min_dl = s.copy()

    print("\tCollecting edge probabilities...")
    mcmc_equilibrate(state, force_niter = force_niter,
                     mcmc_args=dict(niter = niter),
                     callback = collect_marginals)

    print("\tGathering data...")
    if method == "conditional_missing":
        ps = []
        for i in range(len(all_missing)):
            ps.append(logsumexp([ep[i] for ep in edges_prob]) - np.log(len(edges_prob)))

        pred_dict = {"v1": [], "v2": [],
                    "p": np.exp(ps),
                    "v1_names": [], "v2_names": [],
                    "edge_type": []}
        col_names = ["v1", "v2", "p", "v1_names", "v2_names",
                    "edge_type"]
        for v1, v2 in all_missing:
            pred_dict["v1"].append(v1)
            pred_dict["v2"].append(v2)
            pred_dict["v1_names"].append(g.vp.names[v1])
            pred_dict["v2_names"].append(g.vp.names[v2])
            if g.edge(v1, v2):
                pred_dict["edge_type"].append(g.ep.edge_type[g.edge(v1, v2)])
            else:
                pred_dict["edge_type"].append("reconstructed")

        df = pd.DataFrame(data = pred_dict, columns = col_names)
        # Add documented edges to df
        doc_df = elist
        doc_df = doc_df.drop(columns = ["n", "x"])
        doc_df = doc_df.rename(columns = {"row_names": "v1_names", "col_names": "v2_names"})
        doc_df = doc_df[doc_df["edge_type"] == "documented"]
        doc_df["p"] = 1
        df = pd.concat([df, doc_df[list(df.columns)]], ignore_index = True)

        res_dict = {"graph": g,
                    "state": state,
                    "state_min_dl": state_min_dl,
                    "pred_probs": df}
    elif method == "marginal_all":
        res_dict = {"graph": g,
                    "state": state,
                    "state_min_dl": state_min_dl,
                    "pred_graph": u}
        res_dict["pred_probs"] = get_predicted_network(res_dict)

    res_dict["min_dl"] = state_min_dl.entropy()

    return(res_dict)')
}


get_groups <- function(){
    return('
def get_groups(state, g, name = None, save_path = None):
    dl = state.entropy()
    nested_state = state.get_block_state()
    bstack = nested_state.get_bstack()

    # Save group membership of nodes
    for i in range(len(bstack)):
        projected_partition = nested_state.project_partition(i, 0).a
        # Get node ids from first projected partition
        if i == 0:
            nodes = np.arange(0, len(projected_partition))
            levels = {"node": nodes}
            col_names = ["node"]

        col_names.append(f"G{i + 1}")
        levels[col_names[i + 1]] = projected_partition

    col_names.append("name")
    levels["name"] = list(g.vp.names)

    df = pd.DataFrame(data = levels, columns = col_names)

    return(df)')
}


get_predicted_network <- function(){
    return('
def not_estimated_edges(g, u):
    missing = []
    for e in g.edges():
        if not u.edge(e.source(), e.target()):
            missing.append((int(e.source()),
                            int(e.target())))
    return(missing)

def get_predicted_network(res_dict):
    u = res_dict["pred_graph"]
    g = res_dict["graph"]
    not_estimated = not_estimated_edges(g, u)
    if len(not_estimated) > 0:
        print(f"\tA total of {len(not_estimated)} defined edges not found on marginal graph.")

    pred_dict = {"v1": [], "v2": [], "p": [],
                 "v1_names": [], "v2_names": [],
                 "edge_type": []}
    col_names = ["v1", "v2", "p", "v1_names", "v2_names",
                 "edge_type"]

    # Data for estimated edges
    eprob = u.ep.eprob
    for e in u.edges():
        s = e.source()
        t = e.target()
        pred_dict["v1"].append(int(s))
        pred_dict["v2"].append(int(t))
        pred_dict["p"].append(eprob[e])
        if not g.edge(s, t):
            # New not defined edge
            pred_dict["v1_names"].append(g.vp.names[s])
            pred_dict["v2_names"].append(g.vp.names[t])
            pred_dict["edge_type"].append("reconstructed")
        else:
            pred_dict["v1_names"].append(g.ep.v1_names[g.edge(s, t)])
            pred_dict["v2_names"].append(g.ep.v2_names[g.edge(s, t)])
            pred_dict["edge_type"].append(g.ep.edge_type[g.edge(s, t)])

    # Not estimated edges
    for e in g.edges():
        s = e.source()
        t = e.target()
        if not u.edge(s, t):
            pred_dict["v1"].append(int(s))
            pred_dict["v2"].append(int(t))
            pred_dict["p"].append(0)
            pred_dict["v1_names"].append(g.ep.v1_names[g.edge(s, t)])
            pred_dict["v2_names"].append(g.ep.v2_names[g.edge(s, t)])
            pred_dict["edge_type"].append(g.ep.edge_type[g.edge(s, t)])

    predicted_df = pd.DataFrame(data = pred_dict, columns = col_names)

    return(predicted_df)')
}

save_pickle <- function(){
return(
'
def save_pickle(res_dict, i):
    pkl_file = f"hsbm_res_fold{i}.pkl"
    with open(pkl_file, "wb") as pkl:
        pickle.dump(res_dict, pkl)
    print(f"\tPython pickle saved as {pkl_file}.")
    return(0)
')

}

save_plot <- function(){
    return(
    '
def save_plot(nested_state, i):
    plt_name = f"hierarchical_plot_fold{i}.pdf"
    nested_state.get_block_state().draw(output = plt_name)
    print(f"\tHierarchical edge bundling plot saved as {plt_name}.")
    return(0)
    ')
}
