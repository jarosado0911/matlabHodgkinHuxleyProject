function [M,nLst,bLst,brchLst,numNodes,numEdges,meanEdge,maxEdge,minEdge,medEdge] = ...
    getgraphstructure(filename, plt, verbose, saveOut)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct graph structure from an SWC file and (optionally) plot & save.
% 
% Inputs:
%   filename  : path to .swc
%   plt       : logical, show plots on screen (default: false)
%   verbose   : logical, print detailed messages (default: false)
%   saveOut   : logical, save printed output and plots to ../output/ (default: false)
%
% Outputs:
%   M         : adjacency matrix (sparse)
%   nLst      : neighbor list per node (cell)
%   bLst      : boundary node list (cell of node ids)
%   brchLst   : branching node list (cell of node ids)
%   numNodes  : number of nodes
%   numEdges  : number of edges
%   meanEdge, maxEdge, minEdge, medEdge : edge length stats
%-------------------------------------------------------------------------%
% Written by James Rosado 09/20/2019
% Updated: add saveOut (log + PNGs) and robust plotting behavior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 2 || isempty(plt),      plt = false;     end
    if nargin < 3 || isempty(verbose),  verbose = false; end
    if nargin < 4 || isempty(saveOut),  saveOut = false; end

    % Resolve output directory (relative to current working directory)
    outDir = fullfile('..','output');
    if saveOut && ~exist(outDir,'dir')
        mkdir(outDir);
    end

    % Start logging if requested
    if saveOut
        logFile = fullfile(outDir,'neuron.log');  % (intended spelling)
        % Ensure previous diary is closed, then open
        diary off;
        diary(logFile);
    end

    % ---- Read SWC (only need id, pid, coord) ----
    [~, id, pid, coord, ~, ~] = readswc(filename);

    % s,t edges (skip the first line which is usually a root with pid=-1)
    s = id(2:end); 
    t = pid(2:end);

    % Construct graph and adjacency
    G = graph(s, t);
    M = adjacency(G);

    % ---- Plotting / Saving figures ----
    makeFigures = plt || saveOut;

    % Visibility: if saving only, keep figures off-screen
    figVisibility = 'on';
    if saveOut && ~plt
        figVisibility = 'off';
    end

    % Derive base name for images
    [~, baseName, ~] = fileparts(filename);
    graphPng = fullfile(outDir, [baseName '_graph.png']);
    sparsityPng = fullfile(outDir, [baseName '_sparsity.png']);

    % ---- Custom 3D plot: nodes and edges colored by distance from origin ----
    if makeFigures
        f1 = figure('Visible', figVisibility);
        hold on
        axis equal
        view(3)
    
        % Compute distance of each node from origin (0,0,0)
        dist = sqrt(sum(coord.^2, 2));
    
        % Normalize distances for colormap mapping
        cmap = jet(256);              % choose a colormap
        distNorm = (dist - min(dist)) / (max(dist) - min(dist) + eps);
        colors = cmap(floor(distNorm*255)+1, :);
    
        % Plot nodes as scatter points
        scatter3(coord(:,1), coord(:,2), coord(:,3), ...
                 36, colors, 'filled');   % 36 = point size, adjust as needed
    
        % Plot edges: match line width to point "diameter"
        lw = 2.5;   % line width roughly similar to scatter point size
        for e = 1:height(G.Edges)
            % get endpoints
            n1 = G.Edges.EndNodes(e,1);
            n2 = G.Edges.EndNodes(e,2);
    
            % average color of the two nodes
            c = mean([colors(n1,:); colors(n2,:)], 1);
    
            % plot edge
            plot3([coord(n1,1), coord(n2,1)], ...
                  [coord(n1,2), coord(n2,2)], ...
                  [coord(n1,3), coord(n2,3)], ...
                  '-', 'Color', c, 'LineWidth', lw);
        end
    
        % Colorbar to show distance scale
        colormap(cmap);
        cb = colorbar;
        ylabel(cb, 'Distance from origin (\mum)');
    
        title(sprintf('Neuron Graph Colored by Distance: %s', baseName), ...
              'Interpreter','none');
        xlabel('x (\mum)'); ylabel('y (\mum)'); zlabel('z (\mum)');
    
        if saveOut
            exportgraphics(f1, graphPng, 'Resolution', 200);
        end
        
        % Sparsity (Hines matrix visualization)
        f2 = figure('Visible', figVisibility);
        spy(M + speye(size(M)));  % show diagonal
        title(sprintf('Hines Sparsity: %s', baseName), 'Interpreter','none');
        if saveOut
            exportgraphics(f2, sparsityPng, 'Resolution', 200);
        end
    end

    % ---- Neighbor, boundary, and branch lists ----
    nLst = {}; bLst = {}; brchLst = {};
    for i = 1:height(G.Nodes)
        nbrs = neighbors(G, i);
        nLst{end+1} = nbrs; %#ok<AGROW>
        if verbose
            fprintf('Node %d has %d neighbor(s)\n', i, numel(nbrs));
        end
        if numel(nbrs) == 1
            bLst{end+1} = i; %#ok<AGROW>
        end
        if numel(nbrs) > 2
            brchLst{end+1} = i; %#ok<AGROW>
        end
    end

    % ---- Edge stats ----
    numEdges = height(G.Edges); 
    numNodes = height(G.Nodes);

    edge_lengths = zeros(1, numEdges);
    for i = 2:numNodes
        % NOTE: assumes id/pid are valid row indices into coord
        d = coord(id(i),:) - coord(pid(i),:);
        edge_lengths(i-1) = sqrt(sum(d.^2));
        if verbose
            fprintf('Edge (%d,%d) length = %.6f microns\n', id(i), pid(i), edge_lengths(i-1));
        end
    end

    meanEdge = mean(edge_lengths);
    maxEdge  = max(edge_lengths);
    minEdge  = min(edge_lengths);
    medEdge  = median(edge_lengths);

    if verbose
        fprintf('\n--- Edge Length Summary ---\n');
        fprintf('Average = %.6f microns\n', meanEdge);
        fprintf('Max     = %.6f microns\n', maxEdge);
        fprintf('Min     = %.6f microns\n', minEdge);
        fprintf('Median  = %.6f microns\n', medEdge);
        if saveOut
            fprintf('Saved plots to:\n  %s\n  %s\n', graphPng, sparsityPng);
            fprintf('Log written to:\n  %s\n', logFile);
        end
    end

    % Always stop diary if we turned it on
    if saveOut
        diary off;
    end
end