import networkx as nx
import copy
import random
import sys
import matplotlib.pyplot as plt
from matplotlib import cm
import Colorizer
from datetime import datetime
from pathlib import Path
import pickle
import numpy as np
import csv
from kResiliencePaths import kResiliencePaths
import math
from timeit import default_timer as timer


def kResilienceTrees(s, d, fails, g, version="multiple", file_name=None, draw=False, unranked=True, treeChoice="shortest"):

    sp_length = getTrueShortestPathLength(g.copy(), s, d, fails)

    # we will use a copy of the input graph
    g_copy = copy.deepcopy(g)

    # At the beginning, all the nodes and the edges of the graph are numbered with '0'
    # i.e. they don't belong to any structure yet
    nx.set_edge_attributes(g_copy,"0","attr")
    nx.set_node_attributes(g_copy,"0","attr")
    # The source and the detsination node have their own separate attributes ('s' and 'd' respectively)
    g_copy.nodes[s]["attr"] = "s"
    g_copy.nodes[s]["label"] = "s"
    g_copy.nodes[d]["attr"] = "d"
    g_copy.nodes[d]["label"] = "d"
    try:
        startEDP = timer()

        edge_disjoint_paths = list(nx.edge_disjoint_paths(g_copy,s,d))
    except:
        # if an error is thrown, it means that the source node and the destination node belong
        # to different graph components; the destination node cannot be reached from the source node anymore
        return (True,0,0,[])

    #sort paths long to short
    edge_disjoint_paths.sort(key=lambda x: len(x), reverse=False)

    #give numbers to paths, the higher the shorter
    no_path = 1
    for path in edge_disjoint_paths:
        for i in range(0,len(path)-1):
            if g_copy.nodes[path[i+1]]["attr"] != "d":
                g_copy.nodes[path[i+1]]["attr"] = str(no_path)
            g_copy[path[i]][path[i+1]]["attr"] = str(no_path)

        no_path += 1

    endEDP = timer()

    timeEDP = endEDP - startEDP
    for fail in fails:
        g_copy[fail[0]][fail[1]]["failed"] = True


    if draw:
        Colorizer.colorizeGraph(g_copy.copy(), paths=len(edge_disjoint_paths), file_name=file_name + "-" + version + "-paths")

    startTreeBuilding = timer()

    if version == "multiple":
        makeMultipleTrees(g_copy, edge_disjoint_paths)
    else:
        makeOneTree(g_copy, edge_disjoint_paths, reverse=True) #reverse here, so we make tree from LONGEST path

    endTreeBuilding = timer()
    timeTree = endTreeBuilding - startTreeBuilding

    d_incidents = set()
    # remove incident edges of D from all structures
    for d_edge in g_copy.edges(d):
        d_incidents.add(d_edge[1])
        g_copy[d_edge[0]][d_edge[1]]['attr'] = "-1"

    s_incidents = set()
    for s_edge in g_copy.edges(s):
        s_incidents.add(s_edge[1])

    if draw:
        Colorizer.colorizeGraph(g_copy.copy(), paths=len(edge_disjoint_paths), file_name=file_name + "-" + version + "-noPP-" + version)

    startTreeProcessing = timer()

    trees_changed, overallNodeAdditions = postProcessTree(g_copy, s, d_incidents, s_incidents, edge_disjoint_paths, version=version)

    endTreeProcessing = timer()
    timeProcessing = endTreeProcessing - startTreeProcessing

    overallTime = (timeEDP + timeTree + timeProcessing) * 1000

    """print("Version:",version)
    print("Nodes:",g.number_of_nodes())
    print("Edges:",g.number_of_edges())
    print("Time:",overallTime)
    print("Time EDP:",timeEDP * 1000)
    print("--------------------------------------------")"""

    if draw:
        Colorizer.colorizeGraph(g_copy.copy(), paths=len(edge_disjoint_paths), file_name=file_name + "-" + version + "-PP-" + version)


    #rankings = rankTree(g_copy, s, d, d_incidents, trees_changed)
    rankings = rankTree(g_copy, s, d, d_incidents) #causes crash when trees neighbors S as well

    hops, routed_paths, tree_order_ranked = routeTrees(g_copy.copy(), s, d, d_incidents, fails, rankings, treeChoice=treeChoice)

    if unranked:
        hops_unranked, routed_paths_unranked, tree_order_unranked = routeTrees(g_copy.copy(), s, d, d_incidents, fails, rankings, unranked=True)
    else:
        hops_unranked = -1

    if draw:
        Colorizer.colorizeGraph(g_copy.copy(), paths=len(edge_disjoint_paths), hops=hops, file_name=file_name + "-PProuted-" + version, routed_paths=routed_paths, sp_length=sp_length)

    return hops, sp_length, overallNodeAdditions, hops_unranked, g_copy, routed_paths, rankings, tree_order_ranked, timeEDP * 1000, overallTime


def makeOneTree(g_copy, edge_disjoint_paths, reverse=False):

    if reverse:
        edge_disjoint_paths.reverse()
        no_tree = len(edge_disjoint_paths)
    else:
        no_tree = 1

    for path in edge_disjoint_paths:
        nodes_added = 0
        for i in range(1, len(path) - 1):
            nodes = [path[i]]  # obtain a list with the nodes of the i-th path
            it = 0
            while (it < len(nodes)):
                # obtain a list with all the incident edges of nodes from the i-th path
                list_of_incident_edges = list(g_copy.edges(nodes[it]))
                # obtain a generator of the previous edges, which provides edges with the '0' attribute
                # (meaning that they are not used by any structure yet)
                edge_candidates_gen = (edge for edge in list_of_incident_edges if
                                       g_copy.get_edge_data(edge[0], edge[1]).get("attr") == "0")
                for edge in edge_candidates_gen:
                    if g_copy.nodes[edge[1]]["attr"] == "0":
                        g_copy[edge[0]][edge[1]]["attr"] = str(no_tree)
                        g_copy.nodes[edge[1]]["attr"] = str(no_tree)

                        nodes.append(edge[1])
                        nodes_added += 1

                    # we also give an attribute to the incident edges of the destination node
                    # however, tree leaves of the tree are considered to be the neighbors of the destination node
                    if g_copy.nodes[edge[1]]["attr"] == "d":
                        g_copy[edge[0]][edge[1]]["attr"] = str(no_tree)

                it += 1

        no_tree = no_tree + 1 if not reverse else no_tree - 1
        if DEBUG:
            print("added", nodes_added, "nodes to tree", no_tree)

def makeMultipleTrees(g_copy, edge_disjoint_paths):
    no_tree = 1
    for path in edge_disjoint_paths:
        nodes_added = 0
        for i in range(1, len(path) - 1):
            nodes = [path[i]]  # obtain a list with the nodes of the i-th path
            it = 0
            while (it < len(nodes)):
                list_of_incident_edges = list(g_copy.edges(nodes[it]))

                edge_candidates_gen = (edge for edge in list_of_incident_edges if
                                       g_copy.get_edge_data(edge[0], edge[1]).get("attr") == "0")
                for edge in edge_candidates_gen:
                    # if g_copy.nodes[edge[1]]["attr"] == "0":

                    node_candidate_incident_attrs = [g_copy[e[0]][e[1]]["attr"] for e in g_copy.edges(edge[1])]

                    if str(no_tree) not in node_candidate_incident_attrs and g_copy[edge[0]][edge[1]]["attr"] == "0" and \
                            g_copy.nodes[edge[1]]["attr"] != "s" and g_copy.nodes[edge[1]]["attr"] != "d":
                        g_copy[edge[0]][edge[1]]["attr"] = str(no_tree)
                        g_copy.nodes[edge[1]]["attr"] = str(no_tree)

                        nodes.append(edge[1])
                        nodes_added += 1

                    # we also give an attribute to the incident edges of the destination node
                    # however, tree leaves of the tree are considered to be the neighbors of the destination node
                    if g_copy.nodes[edge[1]]["attr"] == "d":
                        g_copy[edge[0]][edge[1]]["attr"] = str(no_tree)

                it += 1
        if DEBUG:
            print("added", nodes_added, "nodes to tree", no_tree)
        no_tree += 1

def routeTrees(g, s, d, d_incidents, fails, rankings, unranked=False, treeChoice="shortest"):
    hops = 0

    trees_attributes = []
    for node1,node2,data in g.edges(data=True):
        if data['attr'] not in trees_attributes and int(data['attr']) > 0:
            trees_attributes.append(data['attr'])

    #trees_attributes = [str(el2) for el2 in sorted([int(el1) for el1 in trees_attributes], reverse=True)]
    trees_attributes = getTreeOrder(rankings, treeChoice=treeChoice)

    routed_paths = []

    found = False

    for attr in trees_attributes:
        if found:
            break
        T = nx.Graph() # we reconstruct the tree
        for node1,node2,data in g.edges(data=True):
            if data['attr'] == attr:
                T.add_edge(node1,node2)

        if s not in list(T.nodes):
            continue

        dfs_edge_order_list = list(nx.dfs_labeled_edges(T, s))
        for n1,n2,label in dfs_edge_order_list:
            if label == "nontree" or n1 == n2: # we also remove self-loops
                dfs_edge_order_list.remove((n1,n2,label))

        hops_current_tree = 0

        if not unranked and attr in rankings:
            dfs_edge_order_list = rankDfs(T, s, d, dfs_edge_order_list, rankings[attr])

        final_dfs = removeFails(dfs_edge_order_list, fails)
        for n1,n2,label in final_dfs:
            routed_paths.append((n1,n2))
            hops_current_tree += 1

            if n1 in d_incidents and (str(n1),str(d)) not in fails and (str(d),str(n1)) not in fails and (int(n1),int(d)) not in fails and (int(d),int(n1)) not in fails:
                routed_paths.append((n1, d))
                hops_current_tree += 1
                found = True
                break
            elif n2 in d_incidents and (str(n2), str(d)) not in fails and (str(d), str(n2)) not in fails and (int(n2), int(d)) not in fails and (int(d), int(n2)) not in fails:
                routed_paths.append((n2, d))
                hops_current_tree += 1
                found = True
                break

        hops += hops_current_tree
    if DEBUG:
        if unranked:
            print("UNRANKED: TOOK", hops, "HOPS")
        else:
            print("RANKED: TOOK", hops, "HOPS")

    return hops if found else -1, routed_paths, trees_attributes

def postProcessTree(g, s, d_incidents, s_incidents, edge_disjoint_paths, version=None, tree_attr=None):


    # we will test if we can still reach the destination node; for this we will analyze each subgraph (tree) until we will find the destination node
    trees_attributes = []
    trees_changed = set()
    if tree_attr is None:
        for node1, node2, data in g.edges(data=True):
            if data['attr'] not in trees_attributes:
                if int(data['attr']) > 0:
                    trees_attributes.append(data['attr'])
    else:
        trees_attributes.append(tree_attr)

    #DEBUG
    #trees_attributes = ["1"]

    overallNodeAdditions = 0
    for attr in trees_attributes:
        T = nx.Graph()  # we reconstruct the tree

        for node1, node2, data in g.edges(data=True):
            if data['attr'] == attr:
                T.add_edge(node1, node2)
        if s in list(T.nodes):

            dfs_edge_order_list = list(nx.dfs_labeled_edges(T, s))
            for n1,n2,label in dfs_edge_order_list:
                if label == "nontree" or n1 == n2: # we also remove self-loops
                    dfs_edge_order_list.remove((n1,n2,label))
                    #edges_to_remove.add((node1, node2))

            #print(dfs_edge_order_list)

            good_branch_nodes = set()
            visited_nodes = set()
            visited_nodes.add(dfs_edge_order_list[0][0])
            delete_mode = False

            for i in range(len(dfs_edge_order_list)):

                n1,n2,label = dfs_edge_order_list[i]

                if label == "forward":
                    visited_nodes.add(n2)
                elif label == "reverse":
                    visited_nodes.remove(n2)

                #here maybe n1?
                if label == "forward" or n2 in good_branch_nodes:
                    delete_mode = False

                if delete_mode:

                    if DEBUG:
                        print("edge {},{} set to 0".format(n1, n2))

                    g[n1][n2]["attr"] = "0"
                    g.nodes[n2]["attr"] = "0"

                #neighs = list(T.neighbors(n2))
                if i < len(dfs_edge_order_list) - 1:
                    n1_next, n2_next, label_next = dfs_edge_order_list[i+1]
                    if label == "forward" and label_next == "reverse" and str(n2) not in d_incidents and int(n2) not in d_incidents: #and n2 not in s_incidents:
                        delete_mode = True
                    elif str(n2) in d_incidents or int(n2) in d_incidents: #or n2 in s_incidents:
                        [good_branch_nodes.add(el) for el in visited_nodes]
                        if edge_disjoint_paths is not None:
                            cnt = 0
                            for path in edge_disjoint_paths:
                                if n2 not in path:
                                    cnt += 1
                                else:
                                    break #new
                            if cnt == len(edge_disjoint_paths):
                                overallNodeAdditions += 1
                                trees_changed.add(attr)
                                if DEBUG:
                                    print("found node", n2, "that was not originally present in edge disjoint paths")

    return trees_changed, overallNodeAdditions

def getTreeOrder(rankings, treeChoice="shortest"):

    #order trees by their shortest distance from root node to D
    order = []
    for key in rankings:
        if len(rankings[key][0]) == 0:
            order.append((key, 1)) #check again
        else:
            if treeChoice == "shortest":
                tmpval = np.min(list(rankings[key][0].values())[0])
            elif treeChoice == "average":
                tmpval = np.mean(list(rankings[key][0].values())[0])
            elif treeChoice == "edgeCount":
                tmpval = len(list(rankings[key][0].values()))

            order.append((key, tmpval))

    order.sort(key=lambda x: x[1])
    final_order = [str(el[0]) for el in order]

    return final_order

def rankTree(g, s, d, d_incidents, trees_changed=None):

    #if trees_changed not given, just take all attrs in graph (even the one that are not trees, but still paths)
    if trees_changed is None:
        trees_changed = []
        for node1, node2, data in g.edges(data=True):
            if data['attr'] not in trees_changed and int(data['attr']) > 0:
                trees_changed.append(data['attr'])

    trees_ranked = {}
    for attr in trees_changed:
        T = nx.Graph()  # we reconstruct the tree

        for node1, node2, data in g.edges(data=True):
            if data is not None and 'attr' in data and data['attr'] == attr:
                T.add_edge(node1, node2)


        if s not in list(T.nodes):
            continue

        dfs_edge_order_list = list(nx.dfs_labeled_edges(T, s))
        for n1,n2,label in dfs_edge_order_list:
            if label == "nontree" or n1 == n2: # we also remove self-loops
                dfs_edge_order_list.remove((n1,n2,label))

        branching_dict = {}
        direction_dict = {}
        best_dict = {}
        travelled_nodes = set()
        for i in range(len(dfs_edge_order_list)):
            n1, n2, label = dfs_edge_order_list[i]

            #if we can reach d from here and its the best route so far, lock temporary distances
            if label == "forward" and n2 in d_incidents:
                for key in branching_dict:
                    if key not in best_dict:
                        best_dict[key] = []
                        direction_dict[key] = []

                    if key in travelled_nodes:
                        best_dict[key].append(branching_dict[key] + 1) #plus one bc we need one more hop to d
                        direction_dict[key].append(n2)

            #add current node as travelled
            travelled_nodes.add(n2)

            #remove current node from travelled if we move backwards from a node already visited
            if label == "reverse" and n2 in travelled_nodes:
                travelled_nodes.remove(n2)

            if label == "forward":
                for node in travelled_nodes:
                    if node in branching_dict:
                        branching_dict[node] += 1
                    else:
                        branching_dict[node] = 1

            if label == "reverse":
                for node in travelled_nodes:
                    branching_dict[node] -= 1

        trees_ranked[attr] = (best_dict, direction_dict)
    return trees_ranked

def rankDfs(T, s, d, dfs, ranking):
    successor_dict = nx.dfs_successors(T, s)

    ranking_processed = [{}, {}]

    def getAllSuccessors(node, succs, set):
        if node in succs:
            for succ in succs[node]:
                set.add(succ)
                getAllSuccessors(succ, succs, set)

    for node in ranking[0]:
        if node not in successor_dict:
            continue

        #write direct neighbor to ranking dict instead of very last node on path -> enable local routing
        neighs = successor_dict[node]
        for neigh in neighs:
            all_succs = set()
            getAllSuccessors(neigh, successor_dict, all_succs)
            for succ in all_succs:
                if succ in ranking[1][node]:
                    ranking[1][node][ranking[1][node].index(succ)] = neigh

        #if multiple paths go over same direct neighbor, take shortest one
        ranking_processed[0][node] = []
        ranking_processed[1][node] = []
        for distinct in np.unique(ranking[1][node]):
            indices = [i for i, x in enumerate(ranking[1][node]) if x == distinct]

            tmpmin = np.min([ranking[0][node][idx] for idx in indices])
            ranking_processed[0][node].append(tmpmin)
            ranking_processed[1][node].append(distinct)

        #sort neighbor rankings so that shortest route is first
        hops = ranking_processed[0][node].copy()
        hops_sorted = sorted(hops)
        dirs_sorted = []
        for i in range(len(hops_sorted)):
            dirs_sorted.append(ranking_processed[1][node][hops.index(hops_sorted[i])])
            hops[hops.index(hops_sorted[i])] = -1

        ranking_processed[0][node] = hops_sorted
        ranking_processed[1][node] = dirs_sorted

    for key in ranking_processed[0]:
        i = 0
        swaps_idx = [-1] * len(ranking_processed[0][key])
        for n1, n2, label in dfs:
            if n1 == key and label == "forward":
                swaps_idx[ranking_processed[1][key].index(n2)] = i
            if n1 == key and label == "reverse":
                swaps_idx[ranking_processed[1][key].index(n2)] = (swaps_idx[ranking_processed[1][key].index(n2)], i + 1)

            i += 1

        range_start = min(swaps_idx, key=lambda x: x[0])[0]
        range_end = max(swaps_idx, key=lambda x: x[1])[1]
        #dfs_reordered = [None] * len(dfs)
        dfs_reordered = []
        for si in swaps_idx:
            dfs_reordered.extend(dfs[si[0] : si[1]])

        dfs[range_start : range_end] = dfs_reordered


    return dfs

def removeFails(dfs, fails):

    #fails = [(14,16)]
    for fail in fails:
        idx = 0
        start_idx = -1
        end_idx = -1
        search_mode = False
        failed_edge = (-1,-1)
        for n1,n2,label in dfs:
            if ((n1 == fail[0] and n2 == fail[1]) or (n1 == fail[1] and n2 == fail[0])) and label == "forward":
                search_mode = True
                failed_edge = fail
                start_idx = idx

            if search_mode and ((n1 == failed_edge[0] and n2 == failed_edge[1]) or (n1 == failed_edge[1] and n2 == failed_edge[0])) and label == "reverse":
                end_idx = idx
                break #new

            idx += 1
        if start_idx > -1 and end_idx > -1:
            for i in range(end_idx, start_idx-1, -1):
                del dfs[i]

    return dfs

def getTrueShortestPathLength(g, s, d, fails):
    g.remove_edges_from(fails)
    try:
        spLength = nx.shortest_path_length(g,s,d)
    except nx.NetworkXNoPath:
        spLength = -1
    return spLength



def main(rep, p, no_nodes, failureModel="random", treeChoice="shortest", edgeFailFactorPar=None, failurePercentPar=None, failureDropPar=None, draw=False, graphMLpath=None):

    now = datetime.now()
    dt_string = now.strftime("%Y-%m-%d_%H-%M-%S")

    if draw:
        Path("img/" + dt_string).mkdir(parents=True, exist_ok=True)


    if graphMLpath:
        loadedML = nx.Graph(nx.read_graphml(graphMLpath))
        no_nodes = loadedML.number_of_nodes()
        p = graphMLpath

        csv_name = "output/trees-n{}-ff{}.csv".format(no_nodes, edgeFailFactorPar)
    else:
        csv_name = "output/trees-p{}-n{}-ff{}.csv".format(p, no_nodes, failurePercentPar)


    writemode = 'w' if no_nodes == 25 else 'a'
    with open(csv_name, mode='w') as csv_file, open('runtime.csv', mode=writemode) as runtime_file:
        fieldnames = ["rep", "nrNodes", "nrEdges", "edgeProb", "nrPaths", "nrFails", "failureModel", "failurePercent", "hopsRanked",
                      "hopsUnranked", "spLength", "stretch_mult", "nodeAdditions", "hopsPath", "stretch_path", "hopsRankedOne",
                      "hopsUnrankedOne", "stretch_one", "overallNodeAdditionsOne", "treeChoice"]

        fieldnames_runtime = ['nodes', 'edges', 'edp', 'one', 'mult']
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer2 = csv.DictWriter(runtime_file, fieldnames=fieldnames_runtime)

        writer.writeheader()
        writer2.writeheader()
        while rep!=0:
            print("iteration", rep)

            if graphMLpath:
                g = loadedML
            else:
                #g = nx.erdos_renyi_graph(no_nodes, p)
                g = nx.random_regular_graph(p, no_nodes, seed=None)


            #pickle.dump(g, open('graph.txt', 'wb'))
            #g = pickle.load(open('graph.txt', 'rb'))


            edges = list(g.edges)
            csv_dict = {}

            csv_dict["rep"] = rep
            csv_dict["nrNodes"] = no_nodes
            csv_dict["nrEdges"] = len(edges)
            csv_dict["edgeProb"] = p
            csv_dict["treeChoice"] = treeChoice

            nodes = list(g.nodes)
            s = random.choice(nodes)
            #s = 32
            d = random.choice(nodes)
            #d = 6
            while s == d or d in list(g.neighbors(s)):
                s = random.choice(nodes)
                d = random.choice(nodes)
            try:
                edge_disjoint_paths = list(nx.edge_disjoint_paths(g,s,d))
            except Exception as e:
                print(e)
                continue

            no_edge_disjoint_paths = len(edge_disjoint_paths)
            csv_dict["nrPaths"] = no_edge_disjoint_paths

            no_failed_edges = no_edge_disjoint_paths-1

            csv_dict["failureModel"] = failureModel
            csv_dict["failurePercent"] = -1

            if (edges):
                if failureModel == "random":
                    if edgeFailFactorPar is None:
                        edgesToFail = no_failed_edges
                    else:
                        edgesToFail = int(round(no_failed_edges * edgeFailFactorPar))

                    csv_dict["nrFails"] = edgesToFail

                    """failed_edges=[]
                    i=1
                    while i <= edgesToFail:
                        
                        edge = random.choice(edges)
                        if (edge not in failed_edges):
                            failed_edges.append(edge)
                        else:
                            i -= 1
                        i += 1"""
                    try:
                        failed_edges = random.sample(edges, edgesToFail)
                    except ValueError:
                        continue
                    csv_dict["failurePercent"] = -1

                elif failureModel == "adversarial":
                    d_incidents = list(g.edges(d))

                    if failurePercentPar is None:
                        #failurePercent = 0.8
                        failurePercent = random.randint(30, 90) / 100
                    else:
                        failurePercent = failurePercentPar

                    sampleNr = math.floor(failurePercent * len(d_incidents))
                    failed_edges = random.sample(d_incidents, sampleNr)
                    csv_dict["failurePercent"] = failurePercent
                elif failureModel == "clustered":
                    failed_edges = []
                    incidents = list(g.edges(d))

                    if failurePercentPar is None:
                        failureStart = 0.6
                    else:
                        failureStart = failurePercentPar

                    if failureDropPar is None:
                        failureDrop = 0.3
                    else:
                        failureDrop = failureDropPar

                    csv_dict["failurePercent"] = str(failureStart) + ";-" + str(failureDrop)

                    failurePercent = failureStart
                    while failurePercent > 0.0:
                        sampleNr = math.floor(failurePercent * len(incidents))
                        failed_edges.extend(random.sample(incidents, sampleNr))
                        failed_edges = list(dict.fromkeys(failed_edges)) #remove duplicates

                        next_incidents = []
                        for edge in incidents:
                            next_incidents.extend(list(g.edges(edge[1])))

                        next_incidents = list(dict.fromkeys(next_incidents)) #remove duplicates
                        incidents = next_incidents
                        failurePercent -= failureDrop


                #failed_edges = [(32,28)]

                no_failed_edges = len(failed_edges)

                try:
                    hops, sp_length, overallNodeAdditions, hops_unranked, g_res_mult, routed_paths_mult, rankings_mult, tree_order_ranked_mult, timeEDP1, timeMult \
                        = kResilienceTrees(s, d, failed_edges, g, version="multiple", file_name=dt_string+"/run"+str(rep), unranked=True, draw=draw, treeChoice=treeChoice)

                    if DEBUG:
                        print("took", hops, "hops after", no_failed_edges, "edges failed")

                    csv_dict["hopsRanked"] = hops
                    csv_dict["hopsUnranked"] = hops_unranked
                    csv_dict["spLength"] = sp_length
                    csv_dict["stretch_mult"] = hops - sp_length
                    csv_dict["nodeAdditions"] = overallNodeAdditions

                    hops_one, sp_length, overallNodeAdditions_one, hops_unranked_one, g_res_one, routed_paths_one, rankings_one, tree_order_ranked_one, timeEDP2, timeOne \
                        = kResilienceTrees(s, d, failed_edges, g, version="one", file_name=dt_string+"/run"+str(rep), unranked=True, treeChoice=treeChoice, draw=draw)

                    csv_dict["hopsRankedOne"] = hops_one
                    csv_dict["hopsUnrankedOne"] = hops_unranked_one
                    csv_dict["stretch_one"] = hops_one - sp_length

                    csv_dict["overallNodeAdditionsOne"] = overallNodeAdditions_one

                    csv_dict["overallNodeAdditionsOne"] = overallNodeAdditions_one
                    csv_dict["overallNodeAdditionsOne"] = overallNodeAdditions_one

                    csv_dict_runtime = {}
                    csv_dict_runtime["edp"] = np.max([timeEDP1, timeEDP2])
                    csv_dict_runtime["one"] = timeOne
                    csv_dict_runtime["mult"] = timeMult
                    csv_dict_runtime["nodes"] = g.number_of_nodes()
                    csv_dict_runtime["edges"] = g.number_of_edges()


                    hops_path, g_res_path, routed_paths_paths, paths_attributes = kResiliencePaths(s, d, failed_edges, g)
                    if DEBUG:
                        print("ONLY PATHS: TOOK", hops_path, "HOPS")
                    csv_dict["hopsPath"] = hops_path
                    csv_dict["stretch_path"] = hops_path - sp_length

                    writer.writerow(csv_dict)
                    writer2.writerow(csv_dict_runtime)

                    """if hops_path < hops_one and hops_path != -1 and hops_one != -1:

                        Colorizer.colorizeGraph(g_res_one.copy(), paths=len(edge_disjoint_paths), hops=hops_one,
                                                file_name="debug/one" + str(rep),
                                                sp_length=sp_length, fails=failed_edges)
                        Colorizer.colorizeGraph(g_res_mult.copy(), paths=len(edge_disjoint_paths), hops=hops,
                                                file_name="debug/mult" + str(rep),
                                                sp_length=sp_length, fails=failed_edges)
                        Colorizer.colorizeGraph(g_res_path.copy(), paths=len(edge_disjoint_paths), hops=hops_path,
                                                file_name="debug/path" + str(rep),
                                                sp_length=sp_length, fails=failed_edges)

                        Colorizer.colorizeGraph(g_res_one.copy(), paths=len(edge_disjoint_paths), hops=hops_one,
                                                file_name="debug/one-routed" + str(rep), routed_paths=routed_paths_one,
                                                sp_length=sp_length, fails=failed_edges)
                        Colorizer.colorizeGraph(g_res_mult.copy(), paths=len(edge_disjoint_paths), hops=hops,
                                                file_name="debug/mult-routed" + str(rep), routed_paths=routed_paths_mult,
                                                sp_length=sp_length, fails=failed_edges)
                        Colorizer.colorizeGraph(g_res_path.copy(), paths=len(edge_disjoint_paths), hops=hops_path,
                                                file_name="debug/path-routed" + str(rep), routed_paths=routed_paths_paths,
                                                sp_length=sp_length, fails=failed_edges)
                        print("stop")"""

                except nx.NetworkXNoPath as e:
                    print(e)

            rep -= 1


DEBUG = False

if __name__ == "__main__":

    p = 0.15

    no_nodes = 100
    repeats = 200

    graphMLpaths = ["graphml/AttMpls.graphml", "graphml/Cogentco.graphml", "graphml/Deltacom.graphml", "graphml/GtsCe.graphml", "graphml/Interoute.graphml", "graphml/Oteglobe.graphml"]
    #graphMLpath = graphMLpaths[5]
    #graphMLpath = None
    #treeChoices = ["shortest", "average", "edgeCount"]

    """main(repeats, 8, no_nodes, failureModel="random", treeChoice="edgeCount", edgeFailFactorPar=None, failurePercentPar=0,
         failureDropPar=0, draw=False)"""

    """while no_nodes<=105:
        main(repeats, p, no_nodes, failureModel="random", treeChoice="edgeCount", edgeFailFactorPar=0, failurePercentPar=0,
             failureDropPar=0, draw=False)
        no_nodes += 10"""

    failpercent = 0.1
    while failpercent <= 0.9:
        print(failpercent)
        main(repeats, 8, no_nodes, failureModel="clustered", treeChoice="edgeCount", edgeFailFactorPar=0, failurePercentPar=failpercent, failureDropPar=0.3, draw=False)
        failpercent += 0.1
        failpercent = round(failpercent, 1)


