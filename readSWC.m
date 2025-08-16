function [A,id,pid,coord,r,subset]=readSWC(filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This functin reads in an swc file and outputs
% separate attachements of the cell
% -------------------------------------------------------------------------%
% INPUT: SWC filename
% OUTPUT: A is a matrix with all the information from SWC
%         id is a vector containing all the nodes
%         pid is a vector with the parent id nodes
%         coord is the coordinates [x,y,z] for each node
%         r is the radius at each node
%         subset is the type of node ie soma,dendrite, axon...
% -------------------------------------------------------------------------%
%   Written by James Rosado 09/20/2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is all the attachments in the format
% id, subset, x,y,z, radius, pid
A = dlmread(filename, ' ', 1, 0);
id = A(:,1);

% these are the subsets the id belong
for i=1:length(id)
    if A(i,2)==1
        subset{i}='soma';
    else
        subset{i}='neurite';
    end
end
% these are the coordinates
% usually I don't use this
coord=A(:,3:5);

% these are the radii --> important later 
% for numerical computations, but for now
% I use uniform radii for testing
r=A(:,6);

% these are the parent ids
% can build an edge list with (id, pid)
pid = A(:,7);
end

