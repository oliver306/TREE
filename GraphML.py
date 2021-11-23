import networkx as nx
import matplotlib.pyplot as plt

g = nx.read_graphml("graphml/AttMpls.graphml")

nx.draw(g)
plt.draw()  # pyplot draw()
