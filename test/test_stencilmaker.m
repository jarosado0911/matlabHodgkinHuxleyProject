
filename = '../data/0-2a.CNG.swc';
[A, id, pid, coord, r, subset] = readswc(filename);

r = r*1e-6;
dt = 1e-5;
P = hh_params(); 
[~,~,~,~,n,~,dx,~,~,~]=getgraphstructure(filename,false,false,false); 
[LHS, RHS] = stencilmaker(n,dt,dx,P.R,r,P.C,filename);
% After you compute LHS and RHS, call:
plotstencilheatmaps(LHS, RHS);