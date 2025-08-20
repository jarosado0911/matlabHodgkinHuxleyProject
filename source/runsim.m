clamp_index = 1;
recvect = (105); % for 2 ref.
filename    = '../data/ref2.swc';
dt = 0.5e-5;
% make a figure
figure(1)
hold on
outfolder   = '../output/sbdf2_results';
sbdf2solve(dt,clamp_index,recvect,filename,outfolder,'SBDF2');

outfolder   = '../output/strang0TR_results';
strangsolve('tr',dt,clamp_index,recvect,0,filename,outfolder,'STRANGF0-TR');

outfolder   = '../output/strang0HN_results';
strangsolve('hn',dt,clamp_index,recvect,0,filename,outfolder,'STRANGF0-HN');

outfolder   = '../output/strang0BE_results';
strangsolve('be',dt,clamp_index,recvect,0,filename,outfolder,'STRANGF0-BE');

outfolder   = '../output/strang0FE_results';
strangsolve('fe',dt,clamp_index,recvect,0,filename,outfolder,'STRANGF0-FE');

outfolder   = '../output/strang0RK4_results';
strangsolve('rk4',dt,clamp_index,recvect,0,filename,outfolder,'STRANGF0-RK4');