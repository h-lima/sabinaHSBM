import numpy as np
import pandas as pd
import pickle
import itertools
import re
import os
import sys
from graph_tool.all import *

def get_predicted_network(g, u, reconstructed = True,
                          name = None, save_path = None): 

    if not reconstructed:
        u = g

    gmpd_dict = {"v1": [], "v2": [], "p": [], 
                 "host": [], "parasite": [], 
                 "edge_type": [], "q": []}
    v2_names = ["v1", "v2", "p", "host", "parasite",
                 "edge_type", "q"]

    eprob = u.ep.eprob
    for e in u.edges(): 
        s = e.source()
        t = e.target()
        gmpd_dict["v1"].append(s)
        gmpd_dict["v2"].append(t)
        gmpd_dict["p"].append(eprob[e])
        if not g.edge(s, t): 
            # New not defined edge
            gmpd_dict["host"].append("")
            gmpd_dict["parasite"].append("")
            gmpd_dict["edge_type"].append("reconstructed")
            gmpd_dict["q"].append(0)   
        else:
            gmpd_dict["host"].append(g.ep.host[g.edge(s, t)])
            gmpd_dict["parasite"].append(g.ep.parasite[g.edge(s, t)])   
            gmpd_dict["edge_type"].append(g.ep.edge_type[g.edge(s, t)])    
            gmpd_dict["q"].append(g.ep.q[g.edge(s, t)])

    df = pd.DataFrame(data = gmpd_dict, columns = v2_names)
    

    return(df)


def get_missing_edges(elist_csv, g):
    df = pd.read_table(elist_csv, sep = " ", 
                       names = colnames) 

    edge_types = g.ep.edge_type.get_2d_array(range(1))[0]
    is_doc = [e_doc == "documented" for e_doc in edge_types]
    g_doc = GraphView(g, efilt = is_doc)
    
    v1s = df.v1.unique()
    v2s = df.v2.unique()

    all_edges = itertools.product(v1s, v2s)

    all_missing = [e for e in all_edges if g_doc.edge(e[0], e[1]) is None]

    return(all_missing)

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
    g.vertex_properties["bipartite"] = vbipartite
    
def prediction_summary(predicted_df):
    held_outs = predicted_df[predicted_df["edge_type"] == 
                                            "held_out"]

    documented = predicted_df[predicted_df["edge_type"] == 
                                                        "documented"]
    
    no_nodes = (predicted_df["edge_type"] == "reconstructed") & \
            ~predicted_df["host"].str.contains("node", na = False)
    reconstructed = predicted_df[no_nodes]
    
    recov_held_outs = sum(held_outs["p"] > 0.5)/len(held_outs["p"])
    print(f"\tRecovered held outs: {recov_held_outs}")
    
    documented_true = sum(documented["p"] > 0.5)/len(documented["p"])
    print(f"\tDocumented correctly maintaned: {documented_true}")
    
    new_reconstructed = sum(reconstructed["p"] > 0.5)
    print(f"\tNew interactions: {new_reconstructed}")
    
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
        
    return(g)

