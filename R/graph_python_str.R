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
from graph_tool.all import *')

}

add_taxa_vertex_prop <- function(){

return('
def add_taxa_vertex_prop(g):
    r_names = []
    for r in g.ep.v1_names:
        if r not in r_names:
            r_names.append(r)
    n_rows = len(r_names)

    vtaxa =  g.new_vp("string")
    vbipartite =  g.new_vp("int")
    for v in g.vertices():
        for e in v.all_edges():
            edge = g.edge(e.source(), e.target())
            if v <= (n_rows - 1):
                vbipartite[v] = 0
                vtaxa[v] = g.ep.v1_names[edge]
            else:
                vbipartite[v] = 1
                vtaxa[v] = g.ep.v2_names[edge]
            break

    g.vertex_properties["taxa"] = vtaxa
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
        add_taxa_vertex_prop(g)

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
                 niter = 10, force_niter = 10000,
                 n_default = 1, x_default = 0, alpha = 1, beta = 1,
                 all_missing = None,
                 method = "binary_classifier", save_pickle = False):

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

    if method == "binary_classifier":
        if all_missing == None:
            all_missing = get_missing_edges(df = elist, g = g)

        marginal_sums = np.zeros(len(all_missing))
        def collect_marginals(s):
            nonlocal entropy, marginal_sums, state_min_dl
            edges_prob = s.get_edges_prob(all_missing)
            marginal_sums  = np.add(marginal_sums, edges_prob)
            if s.entropy() < entropy:
                entropy = s.entropy()
                state_min_dl = s.copy()
        print("using binary_classifier")
    elif method == "full_reconstruction":
        print("Computing full reconstruction")
        # Collecting marginal for network reconstruction
        u = None
        i = 0
        def collect_marginals(s):
            nonlocal u, i, entropy, state_min_dl
            u = s.collect_marginal(u)
            # Save min description length state
            if s.entropy() < entropy:
                entropy = s.entropy()
                state_min_dl = s.copy()
            i += 1

    print("\tCollect marginals")
    mcmc_equilibrate(state, force_niter = force_niter,
                     mcmc_args=dict(niter = niter),
                     callback = collect_marginals)

    print("\tGather data")
    if method == "binary_classifier":
        v1s = list()
        v2s = list()
        edge_type = list()
        value = list()
        v1_names = list()
        v2_names = list()
        for v1, v2 in all_missing:
            v1s.append(v1)
            v2s.append(v2)
            v1_names.append(g.vp.taxa[v1])
            v2_names.append(g.vp.taxa[v2])
            if g.edge(v1, v2):
                edge_type.append(g.ep.edge_type[g.edge(v1, v2)])
                value.append(g.ep.value[g.edge(v1, v2)])
            else:
                edge_type.append("reconstructed")
                value.append(0)

        df= pd.DataFrame({"v1": v1s,
                          "v2": v2s,
                          "p": np.exp(marginal_sums/force_niter),
                          "value": value,
                          "v1_names": v1_names,
                          "v2_names": v2_names,
                          "edge_type": edge_type})
        # Add documented edges to df
        doc_df = elist
        doc_df = doc_df.drop(columns = ["n", "x"])
        doc_df = doc_df.rename(columns = {"row_names": "v1_names", "col_names": "v2_names"})
        doc_df = doc_df[doc_df["edge_type"] == "documented"]
        doc_df["p"] = 1
        df = pd.concat([df, doc_df[list(df.columns)]])

        res_dict = {"graph": g,
                    "state": state,
                    "state_min_dl": state_min_dl,
                    "pred_probs": df}
    elif method == "full_reconstruction":
        res_dict = {"graph": g,
                    "state": state,
                    "state_min_dl": state_min_dl,
                    "pred_graph": u}
        res_dict["pred_probs"] = get_predicted_network(res_dict)

    res_dict["min_dl"] = state_min_dl.entropy()

    return(res_dict)')
}


get_groups <- function(){
    return('def get_groups(state, g, name = None, save_path = None):
    dl = state.entropy()
    nested_state = state.get_block_state()
    bstack = nested_state.get_bstack()

    # Save group membership of nodes
    for i in range(len(bstack)):
        projected_partition = nested_state.project_partition(i, 0).a
        # Get node ids from first projected partition
        if i == 0:
            nodes = np.arange(0, len(projected_partition))
            levels = {"nodes": nodes}
            col_names = ["nodes"]

        col_names.append(f"G{i + 1}")
        levels[col_names[i + 1]] = projected_partition

    col_names.append("taxa")
    levels["taxa"] = list(g.vp.taxa)

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
        print(f"\tA total of {len(not_estimated)} defined edges not estimated")
        #print("\t\t", not_estimated)

    gmpd_dict = {"v1": [], "v2": [], "p": [],
                 "v1_names": [], "v2_names": [],
                 "edge_type": [], "n": [], "x": []}
    col_names = ["v1", "v2", "p", "v1_names", "v2_names",
                 "edge_type", "n", "x"]

    # Data for estimated edges
    eprob = u.ep.eprob
    for e in u.edges():
        s = e.source()
        t = e.target()
        gmpd_dict["v1"].append(int(s))
        gmpd_dict["v2"].append(int(t))
        gmpd_dict["p"].append(eprob[e])
        if not g.edge(s, t):
            # New not defined edge
            gmpd_dict["v1_names"].append(g.vp.taxa[s])
            gmpd_dict["v2_names"].append(g.vp.taxa[t])
            gmpd_dict["edge_type"].append("reconstructed")
            gmpd_dict["n"].append(0)
            gmpd_dict["x"].append(0)
        else:
            gmpd_dict["v1_names"].append(g.ep.v1_names[g.edge(s, t)])
            gmpd_dict["v2_names"].append(g.ep.v2_names[g.edge(s, t)])
            gmpd_dict["edge_type"].append(g.ep.edge_type[g.edge(s, t)])
            gmpd_dict["n"].append(g.ep.n[g.edge(s, t)])
            gmpd_dict["x"].append(g.ep.x[g.edge(s, t)])

    # Not estimated edges
    for e in g.edges():
        s = e.source()
        t = e.target()
        if not u.edge(s, t):
            gmpd_dict["v1"].append(int(s))
            gmpd_dict["v2"].append(int(t))
            gmpd_dict["p"].append(0)
            gmpd_dict["v1_names"].append(g.ep.v1_names[g.edge(s, t)])
            gmpd_dict["v2_names"].append(g.ep.v2_names[g.edge(s, t)])
            gmpd_dict["edge_type"].append(g.ep.edge_type[g.edge(s, t)])
            gmpd_dict["n"].append(g.ep.n[g.edge(s, t)])
            gmpd_dict["x"].append(g.ep.x[g.edge(s, t)])

    predicted_df = pd.DataFrame(data = gmpd_dict, columns = col_names)

    #prediction_summary(predicted_df, type_model)

    return(predicted_df)')
}

save_pickle <- function(){
return(
'
def save_pickle(res_dict, i):
    pkl_file = f"hsbm_res_fold{i}.pkl"
    with open(pkl_file, "wb") as pkl:
        pickle.dump(res_dict, pkl)
    print(f"\tPython pickle save as {pkl_file}")
    return(0)
')

}
