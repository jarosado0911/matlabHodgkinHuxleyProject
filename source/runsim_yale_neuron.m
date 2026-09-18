% RUNSIM_YALE_NEURON
% Cross-validates this project's HH cable model against NEURON (the
% "Yale NEURON" simulator, https://www.neuron.yale.edu/neuron/), using the
% exact same channel kinetics as gates.m/hh_params.m (see source/mod/hhcustom.mod)
% so that only the numerical integration scheme differs from sbdf2solve.m.
%
% Produces a single-trace voltage plot in the same style as
% run_sbdf2_convergence.m (title/labels/ylim/legend/png+fig export), but for
% one morphology and one NEURON run instead of a multi-level convergence sweep.
%
% Usage: set SWCFILE below to the morphology you want to simulate, then run
% this script. It will:
%   1) locate a working NEURON install (nrniv), installing one via pip if
%      none can be found,
%   2) compile source/mod/hhcustom.mod into a NEURON mechanism library
%      (only when the .mod source is newer than the cached binary),
%   3) pick the same "farthest branch point from soma" recording site used
%      by run_sbdf2_convergence.m,
%   4) drive NEURON from a generated .hoc script and read back its voltage
%      trace, and
%   5) plot/save the result.

clear; clc;

% ---- USER: morphology to simulate -------------------------------------
swcfile = '../data/ref2.swc';

% ---- USER: NEURON fixed time step (s) ----------------------------------
dt = 1.0e-6;

outfolder = '../output/yale_neuron_results';
if ~exist(outfolder,'dir'), mkdir(outfolder); end

modDir = fullfile(fileparts(mfilename('fullpath')), 'mod');

% ---- parameters (same sources sbdf2solve.m uses) -----------------------
P = hh_params();
S = sim_params(dt);

% ---- locate/install NEURON, compile the custom mechanism --------------
nrn = find_or_install_neuron();
dllPath = build_mechanism(nrn, modDir);

% ---- recording site: farthest branch point from soma, same rule as
%      run_sbdf2_convergence.m ------------------------------------------
[rec_ind, recCoord, recDist] = find_distal_branch(swcfile);
[~, ~, ~, coord, ~, ~] = readswc(swcfile);
clampCoord = coord(1,:);   % node 1 = soma, matches clamp_index=1 elsewhere in this project

fprintf('%s: recording at branch node %d, coord = (%.3f, %.3f, %.3f), distance from soma = %.3f [um]\n', ...
    swcfile, rec_ind, recCoord(1), recCoord(2), recCoord(3), recDist);

% ---- unit conversions: hh_params.m/sim_params.m (MKS) -> NEURON (cgs-ish) ---
conv = struct( ...
    'Ra'    , P.R   * 100 , ...  % ohm*m   -> ohm*cm
    'cm'    , P.C   * 100 , ...  % F/m^2   -> uF/cm^2
    'gnabar', P.gna * 1e-4, ...  % S/m^2   -> S/cm^2
    'gkbar' , P.gk  * 1e-4, ...
    'gl'    , P.gl  * 1e-4, ...
    'ena'   , P.ena * 1e3 , ...  % V       -> mV
    'ek'    , P.ek  * 1e3 , ...
    'el'    , P.el  * 1e3 , ...
    'dt'      , S.dt      * 1e3, ...  % s -> ms
    'tstop'   , S.endTime * 1e3, ...
    'delay'   , S.delay   * 1e3, ...
    'stop'    , S.stop    * 1e3, ...
    'vStart'  , S.vStart  * 1e3, ...
    'vClamp'  , S.vClamp  * 1e3 ...
);

% ---- generate the NEURON driver and run it -----------------------------
hocfile  = fullfile(outfolder, 'run_yale_neuron.hoc');
tracecsv = fullfile(outfolder, 'trace_yale_neuron.csv');

write_hoc_driver(hocfile, swcfile, tracecsv, conv, clampCoord, recCoord);
run_neuron(nrn, dllPath, hocfile);

% ---- read back the trace and plot (single-graph analogue of
%      run_sbdf2_convergence.m) -----------------------------------------
data = readmatrix(tracecsv);
t_ms = data(:,1);
v_mV = data(:,2);

figure; hold on;
plot(t_ms, v_mV, 'DisplayName', sprintf('NEURON, dt=%.3f [ms]', conv.dt));

ylim([-10, 50]);
title(sprintf('NEURON voltage at distal branch point = (%.3f, %.3f, %.3f)', ...
    recCoord(1), recCoord(2), recCoord(3)));
xlabel('time [ms]');
ylabel('voltage [mV]');
legend('show', 'Location', 'best');

pngPath = fullfile(outfolder, 'yale_neuron_trace.png');
figPath = fullfile(outfolder, 'yale_neuron_trace.fig');
exportgraphics(gcf, pngPath, 'Resolution', 200);
savefig(gcf, figPath);

fprintf('Saved figure to %s\n', pngPath);

% ========================================================================
function [distalId, distalCoord, distalDist] = find_distal_branch(filename)
% Finds the branch node (degree > 2 in the undirected neuron graph) that
% is farthest (Euclidean distance) from the soma (node id 1). Copied from
% run_sbdf2_convergence.m so both scripts pick the same physical site.
    [~,~,~,brchLst,~,~,~,~,~,~] = getgraphstructure(filename, false, false, false);
    [~,~,~,coord,~,~] = readswc(filename);

    somaCoord = coord(1,:);
    nB = numel(brchLst);
    dists = zeros(nB,1);
    for i = 1:nB
        dists(i) = norm(coord(brchLst{i},:) - somaCoord);
    end
    [distalDist, idx] = max(dists);
    distalId = brchLst{idx};
    distalCoord = coord(distalId,:);
end

% ========================================================================
function nrn = find_or_install_neuron()
%FIND_OR_INSTALL_NEURON Locate a working NEURON (nrniv) install, attempting
% a pip install as a fallback when none can be found.
%   nrn.nrniv    - full path to the nrniv executable
%   nrn.home     - NEURON install root (parent of bin/)
%   nrn.sh       - POSIX sh used to run the Windows mod-compiler script
%   nrn.compiler - full path to the mod-compiler script (mknrndll/nrnivmodl)
%   nrn.iswin    - true on Windows

    nrn = struct('nrniv','', 'home','', 'sh','', 'compiler','', 'iswin', ispc);

    if ispc
        candidates = { ...
            'C:\nrn\bin\nrniv.exe', ...
            'C:\Program Files\NEURON\bin\nrniv.exe', ...
            'C:\Program Files (x86)\NEURON\bin\nrniv.exe', ...
            fullfile(getenv('USERPROFILE'), 'AppData', 'Local', 'Programs', 'NEURON', 'bin', 'nrniv.exe') ...
        };
        [status, out] = system('where nrniv');
    else
        candidates = {'/usr/local/nrn/bin/nrniv', '/usr/local/bin/nrniv', '/opt/nrn/bin/nrniv'};
        [status, out] = system('which nrniv');
    end

    if status == 0
        lines = strsplit(strtrim(out), newline);
        candidates = [lines(1), candidates];
    end

    found = '';
    for k = 1:numel(candidates)
        if isfile(candidates{k})
            found = candidates{k};
            break;
        end
    end

    if isempty(found)
        fprintf('NEURON (nrniv) not found on this system - attempting install via pip...\n');
        found = try_install_neuron();
    end

    if isempty(found)
        error('runsim_yale_neuron:NeuronNotFound', ['Could not find or install NEURON automatically.\n' ...
              'Install it manually from https://www.neuron.yale.edu/neuron/download\n' ...
              '(or "pip install neuron" with a Python 3 interpreter), then re-run.']);
    end

    nrn.nrniv = found;
    binDir    = fileparts(found);
    nrn.home  = fileparts(binDir);

    if nrn.iswin
        nrn.compiler = fullfile(binDir, 'mknrndll');
        nrn.sh       = fullfile(nrn.home, 'mingw', 'bin', 'sh.exe');
        if ~isfile(nrn.sh)
            error('runsim_yale_neuron:CompilerNotFound', ...
                'Found nrniv at %s but no bundled sh at %s (needed to run mknrndll).', found, nrn.sh);
        end
    else
        nrn.compiler = fullfile(binDir, 'nrnivmodl');
        if ~isfile(nrn.compiler)
            error('runsim_yale_neuron:CompilerNotFound', ...
                'Found nrniv at %s but no nrnivmodl next to it.', found);
        end
    end

    fprintf('Using NEURON: %s\n', nrn.nrniv);
end

% ========================================================================
function found = try_install_neuron()
% Best-effort install of the "neuron" PyPI package (bundles nrniv) using
% whichever Python 3 interpreter is available. Low blast radius: installs
% into the user's own site-packages, nothing system-wide or hard to undo.
    found = '';
    pyCandidates = {'python3', 'python', 'py -3'};
    for k = 1:numel(pyCandidates)
        py = pyCandidates{k};
        [status, ~] = system(sprintf('%s --version', py));
        if status ~= 0
            continue;
        end

        fprintf('  installing "neuron" package with: %s -m pip install --user neuron\n', py);
        [status, out] = system(sprintf('%s -m pip install --user neuron', py));
        fprintf('%s\n', out);
        if status ~= 0
            continue;
        end

        [status, scriptsDir] = system(sprintf( ...
            '%s -c "import sysconfig; print(sysconfig.get_path(''scripts''))"', py));
        if status ~= 0
            continue;
        end
        scriptsDir = strtrim(scriptsDir);

        candidate = fullfile(scriptsDir, 'nrniv.exe');
        if ~isfile(candidate)
            candidate = fullfile(scriptsDir, 'nrniv');
        end
        if isfile(candidate)
            found = candidate;
            return;
        end
    end
end

% ========================================================================
function dllPath = build_mechanism(nrn, modDir)
%BUILD_MECHANISM Compile source/mod/hhcustom.mod into a NEURON mechanism
% library, rebuilding only when the .mod source is newer than the cached
% binary.
    modFile = fullfile(modDir, 'hhcustom.mod');
    if ~isfile(modFile)
        error('runsim_yale_neuron:ModMissing', 'Cannot find %s', modFile);
    end

    if nrn.iswin
        dllPath = fullfile(modDir, 'nrnmech.dll');
    else
        dllPath = fullfile(modDir, 'x86_64', '.libs', 'libnrnmech.so');
    end

    needBuild = true;
    if isfile(dllPath)
        dMod = dir(modFile);
        dDll = dir(dllPath);
        needBuild = dMod.datenum > dDll.datenum;
    end

    if ~needBuild
        fprintf('Using cached NEURON mechanism: %s\n', dllPath);
        return;
    end

    fprintf('Compiling %s ...\n', modFile);
    if nrn.iswin
        binDir = fullfile(nrn.home, 'bin');
        mingwBin = fullfile(nrn.home, 'mingw', 'bin');
        batPath = fullfile(tempdir, 'build_hhcustom.bat');
        fid = fopen(batPath, 'w');
        fprintf(fid, '@echo off\r\n');
        fprintf(fid, 'set N=%s\r\n', nrn.home);
        fprintf(fid, 'set PATH=%s;%s;%%PATH%%\r\n', binDir, mingwBin);
        fprintf(fid, 'cd /d "%s"\r\n', modDir);
        fprintf(fid, '"%s" "%s" .\r\n', nrn.sh, nrn.compiler);
        fclose(fid);
        [status, out] = system(sprintf('"%s"', batPath));
        delete(batPath);
    else
        cmd = sprintf('cd "%s" && "%s" .', modDir, nrn.compiler);
        [status, out] = system(cmd);
    end
    fprintf('%s\n', out);

    if ~isfile(dllPath)
        error('runsim_yale_neuron:CompileFailed', ...
            'Failed to compile %s (status=%d). See output above.', modFile, status);
    end
    fprintf('Built %s\n', dllPath);
end

% ========================================================================
function write_hoc_driver(hocfile, swcfile, tracecsv, conv, clampCoord, recCoord)
%WRITE_HOC_DRIVER Write a self-contained .hoc script that imports swcfile,
% inserts the hhcustom mechanism with parameters from conv, applies a
% near-ideal voltage clamp at clampCoord (mirrors the Dirichlet clamp in
% sbdf2solve.m), records membrane voltage at recCoord, and writes a
% time,voltage CSV to tracecsv.

    swcPath   = strrep(fullpath_abs(swcfile), '\', '/');
    csvPath   = strrep(fullpath_abs(tracecsv), '\', '/');

    txt = sprintf([ ...
        'load_file("stdlib.hoc")\n' ...
        'objref hoc_sf_\n' ...
        'hoc_sf_ = new StringFunctions()\n' ...
        'load_file("import3d.hoc")\n' ...
        '\n' ...
        'objref swc_, imp_\n' ...
        'swc_ = new Import3d_SWC_read()\n' ...
        'swc_.input("%s")\n' ...
        'imp_ = new Import3d_GUI(swc_, 0)\n' ...
        'imp_.instantiate(nil)\n' ...
        '\n' ...
        '// biophysics converted from hh_params.m (MKS) to NEURON units\n' ...
        'forall {\n' ...
        '    Ra = %.10g\n' ...
        '    cm = %.10g\n' ...
        '    insert hhcustom\n' ...
        '    gnabar_hhcustom = %.10g\n' ...
        '    gkbar_hhcustom  = %.10g\n' ...
        '    gl_hhcustom     = %.10g\n' ...
        '    el_hhcustom     = %.10g\n' ...
        '    ena = %.10g\n' ...
        '    ek  = %.10g\n' ...
        '    nseg = int((L/(0.1*lambda_f(100)))/2)*2 + 1\n' ...
        '}\n' ...
        '\n' ...
        '// nearest section/loc to a target 3D coordinate\n' ...
        'objref nearest_sec\n' ...
        'nearest_loc = 0.5\n' ...
        'proc find_nearest() { local i, d, best\n' ...
        '    best = 1e9\n' ...
        '    forall {\n' ...
        '        for (i=0; i < n3d(); i += 1) {\n' ...
        '            d = sqrt((x3d(i)-$1)*(x3d(i)-$1) + (y3d(i)-$2)*(y3d(i)-$2) + (z3d(i)-$3)*(z3d(i)-$3))\n' ...
        '            if (d < best) {\n' ...
        '                best = d\n' ...
        '                nearest_sec = new SectionRef()\n' ...
        '                if (L > 0) {\n' ...
        '                    nearest_loc = arc3d(i)/L\n' ...
        '                } else {\n' ...
        '                    nearest_loc = 0.5\n' ...
        '                }\n' ...
        '            }\n' ...
        '        }\n' ...
        '    }\n' ...
        '}\n' ...
        '\n' ...
        'find_nearest(%.10g, %.10g, %.10g)\n' ...
        'objref clamp_secref\n' ...
        'clamp_secref = nearest_sec\n' ...
        'clamp_loc = nearest_loc\n' ...
        '\n' ...
        'find_nearest(%.10g, %.10g, %.10g)\n' ...
        'objref rec_secref\n' ...
        'rec_secref = nearest_sec\n' ...
        'rec_loc = nearest_loc\n' ...
        '\n' ...
        '// near-ideal voltage clamp at the soma, mirrors the Dirichlet\n' ...
        '// clamp_index update in sbdf2solve.m: free until delay, held at\n' ...
        '// vClamp until stop, then free again. dur3=0 is deliberate: a\n' ...
        '// nonzero-duration SEClamp phase keeps sinking/sourcing current\n' ...
        '// even with amp3=vStart, which weakly re-pins the soma instead of\n' ...
        '// truly freeing it - dur3=0 fully disconnects the clamp after\n' ...
        '// stop, matching sbdf2solve.m''s unclamped branch exactly.\n' ...
        'objref vclamp_\n' ...
        'clamp_secref.sec { vclamp_ = new SEClamp(clamp_loc) }\n' ...
        'vclamp_.dur1 = %.10g\n' ...
        'vclamp_.amp1 = %.10g\n' ...
        'vclamp_.dur2 = %.10g\n' ...
        'vclamp_.amp2 = %.10g\n' ...
        'vclamp_.dur3 = 0\n' ...
        'vclamp_.amp3 = 0\n' ...
        'vclamp_.rs   = 1e-3\n' ...
        '\n' ...
        'objref tvec_, vvec_\n' ...
        'tvec_ = new Vector()\n' ...
        'vvec_ = new Vector()\n' ...
        'tvec_.record(&t)\n' ...
        'rec_secref.sec { vvec_.record(&v(rec_loc)) }\n' ...
        '\n' ...
        'secondorder = 1\n' ...
        'dt = %.10g\n' ...
        'tstop = %.10g\n' ...
        'finitialize(%.10g)\n' ...
        'while (t < tstop) {\n' ...
        '    fadvance()\n' ...
        '}\n' ...
        '\n' ...
        'objref outfile_\n' ...
        'outfile_ = new File()\n' ...
        'outfile_.wopen("%s")\n' ...
        'for i_ = 0, tvec_.size()-1 {\n' ...
        '    outfile_.printf("%%g,%%g\\n", tvec_.x[i_], vvec_.x[i_])\n' ...
        '}\n' ...
        'outfile_.close()\n' ...
        '\n' ...
        'quit()\n' ...
    ], ...
        swcPath, ...
        conv.Ra, conv.cm, conv.gnabar, conv.gkbar, conv.gl, conv.el, conv.ena, conv.ek, ...
        clampCoord(1), clampCoord(2), clampCoord(3), ...
        recCoord(1), recCoord(2), recCoord(3), ...
        conv.delay, conv.vStart, (conv.stop-conv.delay), conv.vClamp, ...
        conv.dt, conv.tstop, conv.vStart, ...
        csvPath);

    fid = fopen(hocfile, 'w');
    if fid == -1
        error('runsim_yale_neuron:WriteFailed', 'Could not write %s', hocfile);
    end
    fwrite(fid, txt);
    fclose(fid);
end

% ========================================================================
function run_neuron(nrn, dllPath, hocfile)
%RUN_NEURON Invoke nrniv on hocfile with the compiled mechanism preloaded.
    cmd = sprintf('"%s" -nobanner -nogui -dll "%s" "%s"', nrn.nrniv, dllPath, hocfile);
    fprintf('Running NEURON...\n');
    [status, out] = system(cmd);
    fprintf('%s\n', out);
    if status ~= 0
        error('runsim_yale_neuron:NeuronRunFailed', ...
            'nrniv exited with status %d. See output above.', status);
    end
end

% ========================================================================
function p = fullpath_abs(p)
% Resolve p (possibly relative) to an absolute path without requiring the
% file to already exist (needed for the not-yet-written CSV output path).
    if ~(numel(p) >= 2 && (p(2) == ':' || (p(1) == '\' && p(2) == '\')))
        p = fullfile(pwd, p);
    end
end
