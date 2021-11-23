import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

def mycmap(cm, i):
    if i > 0:
        return cm.colors[i % len(cm.colors)]
    elif i == 0:
        return "lightgrey"
    elif i == -1:
        return "black"

def colorizeGraph(graph, paths=-1, hops=None, file_name="img", routed_paths=None, sp_length=None, fails=[]):
    pos = nx.kamada_kawai_layout(graph)


    cm = plt.cm.get_cmap('Set1')

    edge_colors = []
    edge_alphas = []
    edge_widths = []
    tree_counter = {}
    for edge in graph.edges:
        data = graph.get_edge_data(*edge)
        attr = data['attr']
        tree_counter[attr] = tree_counter.get(attr, 0) + 1

        edge_colors.append(mycmap(cm, int(attr)))

        if int(attr) > 0:
            edge_alphas.append(1)
        elif int(attr) <= 0:
            edge_alphas.append(0.5)

        if "failed" in data or (edge[1], edge[0]) in fails or (edge[0], edge[1]) in fails:
            if int(attr) == 0:
                edge_colors[-1] = "red"

            edge_widths.append(4)
            edge_alphas[-1] = 0.6
        elif routed_paths and (edge in routed_paths or (edge[1], edge[0]) in routed_paths):
            edge_colors[-1] = "darkgreen"
            edge_widths.append(8)

        else:
            edge_widths.append(0.5)




    node_colors = []
    node_alphas = []
    labels = {}
    for node in graph.nodes:

        if "attr" in graph.nodes[node] and graph.nodes[node]["attr"] != "0":
            node_alphas.append(1)
        else:
            node_alphas.append(0.5)

        if "attr" in graph.nodes[node] and graph.nodes[node]["attr"] == "s":
            node_colors.append("green")
            labels[node] = "S" + str(node)
        elif "attr" in graph.nodes[node] and graph.nodes[node]["attr"] == "d":
            node_colors.append("r")
            labels[node] = "D" + str(node)
        else:
            node_colors.append("lightblue")
            labels[node] = node

    fig = plt.figure()
    for key, val in sorted(tree_counter.items()):
        plt.plot([0], [0], color=mycmap(cm, int(key)), label=val)

    nx.draw_networkx_nodes(graph, pos, node_size=5, node_color=node_colors, alpha=node_alphas)
    nx.draw_networkx_edges(graph, pos, edge_color=edge_colors, width=edge_widths, alpha=edge_alphas)
    nx.draw_networkx_labels(graph, pos, labels, font_size=12)

    nr_nodes = len(graph.nodes)
    nr_edges = len(graph.edges)
    title_str = str(nr_nodes) + " nodes, " + str(nr_edges) + " edges, " + str(paths) + " edge disjoint paths\n"
    if hops:
        title_str += str(hops) + " hops taken"
    if sp_length:
        title_str += ", SP-length:" + str(sp_length)
    if hops and sp_length:
        title_str += ", stretch:" + str(hops-sp_length)

    plt.title(title_str)

    plt.legend(loc="upper right")


    #plt.show()
    fig.savefig('img/' + file_name + '.png', dpi=1000)

    plt.cla()