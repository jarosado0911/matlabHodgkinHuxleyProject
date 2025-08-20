function A = stencilmaker(n,R,a,C,filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stencilMaker  this generates the stencil matrix for the diffusion solve
% of the operator spliting for the Hodgkin Huxley model equations
%   n = number of vertices in graph geometry
%   R = this is the specific resistance
%   a = vector of radii at each vertex
%   C = this is the membrane capacitance which is fixed
%   fileaname = this is the path/geometry.swc
%   Written by James Rosado 09/20/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long
% Initialize LHS and RHS matrices
A=zeros(n);

% get coordinates of the vertices of geometry
[~,~,~,coord,~,~]=readswc(filename);

% get the neigborlist, and number of nodes
[~,nLst,~,~,n,~,~,~,~,~]=getgraphstructure(filename,false,false,false);

% print some information
fprintf('n = %i nodes\n',n)
fprintf('R = %d  [Ohm.m]\n', R)
fprintf('C = %d  [F/m2]\n',C)

for i=1:n
    % get neighborlist
    nghbr = nLst{i};
    lenN = length(nghbr); %length of neighborlist
    
    % get coordinates of current node
    currCoord = coord(i,:);
    nghbrCoord=zeros(lenN,3);
    edgeLengths = zeros(lenN,1);
    sumRecip = 0;
    
    % get coordinates of neighbor nodes
    for j=1:lenN
        nghbrCoord(j,:) = coord(nghbr(j),:);
        edgeLengths(j) = norm(currCoord - nghbrCoord(j,:))*1e-6;  % compute edge lengths to neighbor
        sumRecip = (1/(edgeLengths(j)*a(i)*(a(nghbr(j))^(-2)+a(i)^(-2))))+sumRecip; % compute sum of reciprocals from finite difference scheme
    end
    aveLengths = mean(edgeLengths);
    
    A(i,i) =  -1*((sumRecip)/(R*C*aveLengths));    
    % Set off diagonal Entries
    for j=1:lenN
        A(i,nghbr(j)) = 1/(R*C*aveLengths*edgeLengths(j)*a(i)*(a(nghbr(j))^(-2)+a(i)^(-2)));   
    end
end

% store as sparse matrix to save memory!
A = sparse(A);
end