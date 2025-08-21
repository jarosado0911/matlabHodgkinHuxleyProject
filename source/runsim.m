clamp_index = 1;
recvect = (105); % for 2 ref.
filename    = '../data/ref2.swc';
outfolder   = '../output/sbdf2_results';
dt = 10.0e-5;

[~, id, pid, coord, ~, ~] = readswc(filename);
plotneuron3views(coord, id, pid,sprintf('%s/sbdf2neuron.png',outfolder));

% make a figure
figure(2)
hold on
sbdf2solve(dt,clamp_index,recvect,filename,outfolder,'SBDF2',1);

outfolder   = '../output/strang0TRBDF_results';
strangsolve('trbdf',dt,clamp_index,recvect,0,filename,outfolder,'STRANGF0-TRBDF',1);

outfolder   = '../output/strang0MD_results';
strangsolve('mid',dt,clamp_index,recvect,0,filename,outfolder,'STRANGF0-MD',1);

outfolder   = '../output/strang0TR_results';
strangsolve('tr',dt,clamp_index,recvect,0,filename,outfolder,'STRANGF0-TR',1);

outfolder   = '../output/strang0HN_results';
strangsolve('hn',dt,clamp_index,recvect,0,filename,outfolder,'STRANGF0-HN',1);

outfolder   = '../output/strang0BE_results';
strangsolve('be',dt,clamp_index,recvect,0,filename,outfolder,'STRANGF0-BE',1);

outfolder   = '../output/strang0FE_results';
strangsolve('fe',dt,clamp_index,recvect,0,filename,outfolder,'STRANGF0-FE',1);

outfolder   = '../output/strang0RK4_results';
strangsolve('rk4',dt,clamp_index,recvect,0,filename,outfolder,'STRANGF0-RK4',1);

outfolder = '../output/sbdf2_results';
% makematlabmovtraces(outfolder,filename,sprintf('%s/sbdf2video',outfolder),recvect(1));