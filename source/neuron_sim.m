function neuron_sim(dt,filename)

    % global radius variable
    global a  

    % get radius, subset from readswc function 
    [~,~,~,~,a,subset]=readswc(filename); 
    
    % scale radii to micron
    a = a*1e-6;

    % for now using M adjacency matrix, neighborlist, boundarylist, branchpoint
    % list, number of nodes n, set average edgelength to be dx
    % [M,nLst,bLst,brchLst,numNodes,numEdges, meanEdge,maxEdge,minEdge,medEdge]
    [~,nLst,blst,brchlst,n,nEdges, dx_mean,dx_max,dx_min,dx_med]=getgraphstructure(filename,true,false,true); 

    % this is for recording the voltages and states at branch points
    record_index = cell2mat(brchlst);

    % this outputs some general information abou the cell to the console
    fprintf('\nNeuron Information....')
    fprintf('\n Num nodes:      %s\n',num2str(n));
    fprintf(' Num edges:      %s\n',num2str(nEdges));
    fprintf(' Branch nodes:   %s\n', num2str(record_index));
    fprintf(' Boundary nodes: %s\n', num2str(cell2mat(blst)));
    fprintf('mean dx = %f microns\n', dx_mean);
    fprintf('maxl dx = %f microns\n', dx_max);
    fprintf('minl dx = %f microns\n', dx_min);
    fprintf('medl dx = %f microns\n', dx_med);
    fprintf('time dt = %f seconds\n', dt);

    % I really do not use dx because the stencil has to compute the dx for each
    % edge when making the stencil entries
    dx = dx_mean*1e-6; 
end
