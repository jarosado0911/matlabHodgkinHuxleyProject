recvect = [1,97,153];
recvect = [recvect(2)];
clamp_index = recvect(1);
filename    = '../data/ref1.swc';
outfolder   = '../output';
% make a figure
figure(1)
hold on
sbdf2solve(1e-5,clamp_index,recvect,filename,outfolder);
strangsolve(1e-5,clamp_index,recvect,0,filename,outfolder);
strangsolve(1e-5,clamp_index,recvect,1,filename,outfolder);

recvect = [1,9,125]; % for 2 ref.
recvect = [recvect(2)];
filename    = '../data/ref2.swc';
outfolder   = '../output';
% make a figure
figure(2)
hold on
sbdf2solve(1e-5,clamp_index,recvect,filename,outfolder);
strangsolve(1e-5,clamp_index,recvect,0,filename,outfolder);
strangsolve(1e-5,clamp_index,recvect,1,filename,outfolder);

% 
% recvect = [1,385,609]; % for 3 ref.
% filename    = '../data/ref3.swc';
% outfolder   = '../output';
% sbdf2solve(1e-5,clamp_index,recvect,filename,outfolder);