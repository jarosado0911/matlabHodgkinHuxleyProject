function [A, id, pid, coord, r, subset] = readswc(filename)
% READSWC Robust SWC reader tolerant to leading spaces/tabs and long headers.
% Output:
%   A:     [N x 7] double [id, type, x, y, z, r, pid]
%   id:    [N x 1] int32
%   pid:   [N x 1] int32
%   coord: [N x 3] double
%   r:     [N x 1] double
%   subset:{N x 1} labels
%-------------------------------------------------------------------------%
% Written by James Rosado 09/20/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fid = fopen(filename, 'r');
    if fid == -1
        error('readswc:OpenFailed','Cannot open file: %s', filename);
    end

    % textscan will ignore any line starting with '#'
    C = textscan(fid, '%f %f %f %f %f %f %f', ...
        'Delimiter', ' \t', ...            % spaces and tabs
        'MultipleDelimsAsOne', true, ...   % collapse runs of delimiters
        'CommentStyle', '#', ...           % skip header lines
        'CollectOutput', true);            % get a single numeric array
    fclose(fid);

    if isempty(C) || isempty(C{1})
        error('readswc:EmptyData','No numeric SWC data found in "%s".', filename);
    end

    A = C{1};
    if size(A,2) < 7
        error('readswc:BadFormat','Expected 7 columns, found %d in "%s".', size(A,2), filename);
    elseif size(A,2) > 7
        A = A(:,1:7);
    end

    % Extract columns
    id    = int32(round(A(:,1)));
    typ   = A(:,2);
    coord = A(:,3:5);
    r     = A(:,6);
    pid   = int32(round(A(:,7)));

    % Map SWC type codes to labels
    subset = repmat({"unknown"}, numel(id), 1);
    subset(typ==1) = {"soma"};
    subset(typ==2) = {"axon"};
    subset(typ==3) = {"dendrite"};
    subset(typ==4) = {"apical_dendrite"};
    subset(typ==5) = {"custom"};
    subset(typ==6) = {"undefined"};
    subset(typ==7) = {"glia"};
end
