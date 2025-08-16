addpath('../source/');

filename = '../data/0-2a.CNG.swc';
plt      = true; 
verbose  = true; 
saveout  = true;

% Show plots on screen, print details, and save logs/PNGs into ../output/
[M,nLst,bLst,brchLst,numNodes,numEdges,meanEdge,maxEdge,minEdge,medEdge] = ...
    getgraphstructure('../data/0-2a.CNG.swc', plt,verbose,saveout);