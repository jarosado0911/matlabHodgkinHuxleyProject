function makematlabmov(dataFolder,filename,outname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makematlabmov this makes a movie of the voltage data from the neuron sim
%  Input:
%  dataFolder: this is the folder that contains the vm voltage files
%  geometryfile: this is the folder/file.swc of the geometry
%  movieName: what do you want to call the movie?
%
%   Written by James Rosado 09/20/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get coordinates of the geometry, the vertices 
[~,~,~,coords,r,~]=readswc(filename);

% Load time data from MatLab simulation
filename=sprintf('%s/time.mat',dataFolder);

% get all the times for showing in plot
load(filename,'t'); t = t';

% this will make the geometry appear realistic by demonstrating the radii
% are not uniform, in particular the cell is thicker near the soma and
% thinner at the ends of dendrites
markerSize = (r./max(r)).*50;

% make a new figure windows
fig=figure('units','normalized','outerposition',[0 0 0.325 1.0]);
% this is for recording the movie
v = VideoWriter(sprintf('%s.mp4',outname),'MPEG-4');
open(v)

% this part is for setting the text values on the color bar, no need to
% modify this except maybe the cmax and cmin if your action potentials have
% lower/higher peak values
yticklabel={};
cmax = 50; cmin = -5;
vals = [cmin:5:cmax];
for i=1:length(vals)
    yticklabel{i}=num2str(vals(i));
end

% if you change the 100 this will affect the length of the movie
for i=0:1:length(t)-1
        % read the voltage data from the .dat files
        u_sol = readmatrix(sprintf('%s/data/vm_t%i.dat',dataFolder,i));
        
        % make a scatter plot
        scatter3(coords(:,1),coords(:,2),coords(:,3),markerSize,'filled','CData',u_sol);
        
        %set labels
        xlabel(sprintf('{\\mu}m'))
        ylabel(sprintf('{\\mu}m'))
        set(gca,'Color', [0.5 0.5 0.5])
        caxis([cmin cmax]*1e-3)
        title(sprintf('MatLab, t = %0.2f [ms]',t(i+1)*1e3))
        colormap('jet')
        colorbar
        
        % set tick labels on colorbar
        c = colorbar;  
        c.Label.String="[mV]";
        c.TickLabels = yticklabel; 
        view(2)
        
        % save the frame to video file
        thisframe=getframe(fig);
        writeVideo(v, thisframe);

        drawnow
        fprintf('frame = %i\n',i)
end
% don't forget to close the video file, if not it will be corrupted!
close(v)