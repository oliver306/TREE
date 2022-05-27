# TREE (Tree Routing to Extend Edge disjoint paths)

## Overview

The repository contains the implementation of the TREE algorithms as outlined in https://arxiv.org/abs/2111.14123

The algorithms are divided as follows:

1. Edge-disjoint paths from source node and destination node (which are known/provided)
2. One Tree based on the edge-disjoint paths described in point 1
3. Multiple Trees based on the edge-disjoint paths described in point 1

The src folder contains 3 Python files:

1. kResiliencePaths.py - The reference edge-disjoint-path implementation
2. kResilienceTrees.py - The extension of these paths to network spanning TREEs + the required logic to route along them
3. Colorizer.py - Uses to visualize and inspect the tree constructions and routing results

This work builds upon the prior work of Paula-Elena Gheorghe, who laid the groundwork for the presented algorithms.

For experiments and tests one can either generate random graphs (by using, for example, the corresponding methods provided by NetworkX) or use topologies from the topology zoo, which can be found at http://www.topology-zoo.org/dataset.html. In order to run the experiments/tests, one must execute the corresponding Python file (kResiliencePaths.py or kResilienceTrees.py) by calling the corresponding function in main and giving all the necessary arguments. 

Evaluation code and results can be found in the form of Jupyter Notebook files inside the evaluation folder. All CSVs, plots and evaluation code used in the published paper can be found inside this folder.




