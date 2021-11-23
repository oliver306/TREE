import networkx as nx
import random
import copy
import sys
import Colorizer
import matplotlib.pyplot as plt

DEBUG = False

# undirected graph g
# start node s
# destination node s
# failed edges fails

def kResiliencePathsPrecomp(g):
    return g

def kResiliencePaths(s,d,fails,g):

    # we will use a copy of the input graph
    g_copy = copy.deepcopy(g)
    # At the beginning, all the the edges of the graph are numbered with '0'
    # i.e. they don't belong to any structure yet
    nx.set_edge_attributes(g_copy,"0","attr") 
    try:
        edge_disjoint_paths = list(nx.edge_disjoint_paths(g_copy,s,d))
        edge_disjoint_paths.sort(key=lambda x: len(x), reverse=False)
    except:
        # if an error is thrown, it means that the source node and the detsination node belong 
        # to different graph components; the destination node cannot be reached from the source node anymore
        return (True,0,0,[])

    # Each edge from every edge-disjoint path will receive a corresponding number as an attribute
    # e.g. all the mentioned edges from the first path will receive '1' as an attribute,
    # meaning that they now belong to a structure (namely the first path)
    no_path = 1
    paths_attributes = []
    for path in edge_disjoint_paths:
        for i in range(0,len(path)-1):
            g_copy[path[i]][path[i+1]]["attr"] = str(no_path)
            if str(no_path) not in paths_attributes:
                paths_attributes.append(str(no_path))
        if DEBUG:
            print("Path no.",no_path,": ",path)
        no_path += 1

    d_incidents = set()
    for d_edge in g_copy.edges(d):
        d_incidents.add(d_edge[1])

    #rankings = kRes.rankTree(g_copy, s, d, d_incidents) #causes crash when trees neighbors S as well
    #paths_attributes = kRes.getTreeOrder(rankings, treeChoice="shortest")

    hops, route = routing_kResiliencePaths(s,d,fails,g_copy.copy(),paths_attributes=paths_attributes)
    # (True/False, hops, 0, [])
    return hops, g_copy, route, paths_attributes
    

def routing_kResiliencePaths(s,d,fails,g_copy, paths_attributes=None):
    result = True # we will find out if the graph can guarantee resiliency against failures or not
    hops = 0 # initialize the number of hops/number of traversed edges during the routing process 
    
    # we remove the failed edges and the edges which do not belong to any structure (i.e. with the '0' attribute) from the graph
    for e in list(g_copy.edges):
        if (e[0],e[1]) in fails or (e[1],e[0]) in fails:
            g_copy.remove_edge(e[0],e[1])
        try:
            if g_copy.get_edge_data(e[0], e[1]).get("attr") == "0":
                g_copy.remove_edge(e[0],e[1])
        except:
            pass

    # we will test if we can still reach the destination node; for this we will analyze each subgraph (path) until we will find the destination node
    if paths_attributes is None:
        paths_attributes = []
        for node1,node2,data in g_copy.edges(data=True):
            if data['attr'] not in paths_attributes:
                paths_attributes.append(data['attr'])


    attributes_gen = (attr for attr in paths_attributes if result) # continue generating path attributes until the destination node will be reached
    #print("pa:", paths_attributes)

    route = []
    found = False
    for attr in attributes_gen:
        T = nx.Graph() # we reconstruct the path
        for node1,node2,data in g_copy.edges(data=True):
            if data['attr'] == attr:
                T.add_edge(node1,node2)

        if DEBUG:
            print("(Updated) Path no.",attr,": (after some possible failures)",list(T.edges))
        if s in list(T.nodes):
            # we obtain a DFS order of the edges of the current path
            # the edges are labeled with the attributes ‘forward’, ‘nontree’, and ‘reverse’
            # 'nontree'-labeled edges should be removed from the dfs_edge_order_list (because the edge is not in the DFS tree)
            dfs_edge_order_list = list(nx.dfs_labeled_edges(T, s))
            for n1,n2,label in dfs_edge_order_list:
                if label == "nontree" or n1 == n2: # we also remove self-loops
                    dfs_edge_order_list.remove((n1,n2,label))

            dfs_edges_gen = (edge_dfs for edge_dfs in dfs_edge_order_list if result)
            hops_current_path = 0
            # we check if we can still reach the destiantion node by using the current path
            for edge_dfs in dfs_edges_gen:
                route.append((edge_dfs[0],edge_dfs[1]))
                hops_current_path += 1
                if str(edge_dfs[1]) == str(d) or str(edge_dfs[0]) == str(d): # the destination node has been found
                    result = False
                    found = True

            hops += hops_current_path
            
    # we check if the algorithm functioned correctly, i.e. if the destination node really cannot be reached anymore
    if DEBUG and result == True:
        error = True
        try:
            R = nx.algorithms.shortest_paths.generic.shortest_path(g_copy,s,d)
        except:
            error = False
        assert error == False
        print('There is no path from source node %d to destination node %d anymore.' % (s , d))

    return hops if found else -1, route

# Uncomment in order to execute a test using a randomly generated graph
# A command-line argument (the number of repetitions for the experiment) can be given


if __name__ == "__main__":

    if len(sys.argv) > 1:
        rep = int(sys.argv[1])
    else:
        rep = 1 # default number of repetitions
    
    while rep!=0:
        # input & output files
        fo_name = "kResiliencePaths_output" + str(rep) + ".txt" # the results will be written here
        fi_name = "kResiliencePaths_input" + str(rep) + ".txt" # details about the randomly-generated graph will be written here
        fo = open(fo_name, "w")
        fi = open(fi_name, "w")
        w = "kResiliencePaths output"+ str(rep) + "\n"
        fo.write(w)
        # Test using a randomly generated graph
        p = 0.019 # Probability of edge appearance
        no_nodes = 400
        g = nx.erdos_renyi_graph(no_nodes, p)



        edges = list(g.edges)
        if edges:
            for edge in edges:
                fi.write(str(edge)+"\n")
        fo.write("Number of nodes: " + str(no_nodes) + "\n")
        fo.write("Number of edges: " + str(len(edges)) + "\n")
        fo.write("Probability of edge appearance: " + str(p) + "\n")
        nodes = list(g.nodes)
        s = random.choice(nodes)
        d = random.choice(nodes)
        while s == d:
            s = random.choice(nodes)
            d = random.choice(nodes)
        fi.write(str(s)+"\n"+str(d)+"\n")
        try:
            edge_disjoint_paths = list(nx.edge_disjoint_paths(g,s,d))
            no_edge_disjoint_paths = len(edge_disjoint_paths)
            print(edge_disjoint_paths[0])
            fo.write("Number of edge-disjoint paths: "+str(no_edge_disjoint_paths) + "\n")
            no_failed_edges = no_edge_disjoint_paths-1
            if (edges):
                failed_edges=[]
                i=1
                while i <= no_failed_edges:
                    edge = random.choice(edges)
                    if (edge not in failed_edges):
                        failed_edges.append(edge)
                    else:
                        i -= 1
                    i += 1
                fi.write(str(failed_edges))
                no_failed_edges = len(failed_edges)
                fo.write("Number of failed edges: " + str(no_failed_edges) + "\n")
                result = kResiliencePaths(s, d, failed_edges, g)
                fo.write(str(result))
                # k edge-disjoint paths assure resiliency against k-1 failures
                if no_failed_edges == no_edge_disjoint_paths-1:
                    assert result[0] == False
            else:
                fo.write("The randomly generated graph has no edges.")
        except:
            fo.write("The randomly generated graph is not connected (and the source not and the destination node belong to different graph components)")
        fi.close()
        fo.close()
        rep -= 1

