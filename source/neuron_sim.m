function neuron_sim(dt,filename,sv)

    % global radius variable
    global a  

    % get radius, subset from readswc function 
    [~,~,~,~,a,subset]=readswc(filename); 
    
    % scale radii to micron
    a = a*1e-6;

    % for now using M adjacency matrix, neighborlist, boundarylist, branchpoint
    % list, number of nodes n, set average edgelength to be dx
    % [M,nLst,bLst,brchLst,numNodes,numEdges, meanEdge,maxEdge,minEdge,medEdge]
    [~,~,blst,brchlst,n,nEdges, dx_mean,dx_max,dx_min,dx_med]=getgraphstructure(filename,true,false,true); 

    % this is for recording the voltages and states at branch points
    record_index = cell2mat(brchlst);

    % this outputs some general information abou the cell to the console
    fprintf('\nNeuron Information....\n')
    fprintf(' Num nodes:        %s\n',num2str(n));
    fprintf(' Num edges:        %s\n',num2str(nEdges));
    fprintf(' Branch nodes:     %s\n', num2str(record_index));
    fprintf(' Boundary nodes:   %s\n', num2str(cell2mat(blst)));
    fprintf(' mean dx =         %f microns\n', dx_mean);
    fprintf(' maxl dx =         %f microns\n', dx_max);
    fprintf(' minl dx =         %f microns\n', dx_min);
    fprintf(' medl dx =         %f microns\n', dx_med);
    fprintf(' time dt =         %f seconds\n', dt);

    % I really do not use dx because the stencil has to compute the dx for each
    % edge when making the stencil entries
    dx = dx_mean*1e-6; 

    % Set the Hodgkinx Huxley Parameters
    global P
    P = hh_params();                    
    print_hh_params(P);
    
    % Set the simulation Parameters
    S = sim_params(dt);         
    print_sim_params(S);

    %------------------------Initialize solution space--------------------%
    u=ones(n,1).*S.vStart;  % this is our voltage solution vector
    somaId=[];
    % find the soma for voltage clamp, it may not always be the ind = 1
    % could be multiple indices in fact, but for our VR it is usually i = 1
    for i=1:n
        if strcmp(subset{i},'soma')
            somaId = [somaId,i];
        end
    end
    clamp = somaId;
    
    % Here I initialize the gating variables m,n, and h
    % these are our state variable solution vectors
    nn=zeros(n,1); nn(:,1)=S.ni;
    mm=zeros(n,1); mm(:,1)=S.mi;
    hh=zeros(n,1); hh(:,1)=S.hi;

    % Make output folder
    dir = sprintf("../output/results/");
    if ~exist(dir), mkdir(dir); end
    if ~exist(dir+"data/") 
        mkdir(dir+"data/"); 
    end

    %---------------------------------------------------------------------%
    % Make sparse stencil matrices
    [LHS, RHS] = stencilmaker(n,dt,dx,P.R,a,P.C,filename);

    % initialize empty recording vectors
    usoma=[]; rec_u =[]; rec_h =[]; rec_m =[]; rec_n =[];

    gN = make_gate_rhs(@an, @bn);   % for n-gate
    gM = make_gate_rhs(@am, @bm);   % for m-gate
    gH = make_gate_rhs(@ah, @bh);   % for h-gate

    fprintf("\nStarting Solve ... \n");
    for i=0:S.nT
        fprintf("  Solving at t = %f\n",i*dt);
        
        % set soma to clamp
        if i*dt >=S.delay
            u(clamp)=S.vClamp;
        end

        % if save is 1 then save the current voltage of entire cell to file
        if sv == 1
            writematrix(u,sprintf('%sdata/vm_t%i.dat',dir,i))
        end

        % we will always record and output the data at the soma and branch
        % points, these locations are important for measuring
        usoma=[usoma;u(clamp)];
        rec_u = [rec_u, u(record_index)];
        rec_h = [rec_h, hh(record_index)];
        rec_m = [rec_m, mm(record_index)];
        rec_n = [rec_n, nn(record_index)];

        % The scheme is a Strang splitting
        % (1) Do a 1/2 time step with the ODEs and Reaction part of OP-split
        nn = RK4(u, nn, dt/2, gN);
        mm = RK4(u, mm, dt/2, gM);
        hh = RK4(u, hh, dt/2, gH);
        u  = RK4_react(u, mm, hh, nn, dt/2, @react,P);
        
        % (2) Do full step of diffusion solve, this is part of OP-split
        u = LHS\(RHS*u);

        % (3) Finish off with other half to time step from ODEs and Reaction
        % part
        u  = RK4_react(u, mm, hh, nn, dt/2, @react,P);
        nn = RK4(u, nn, dt/2, gN);
        mm = RK4(u, mm, dt/2, gM);
        hh = RK4(u, hh, dt/2, gH);

        % set soma to clamp
        if i*dt >=S.delay
            u(clamp)=S.vClamp;
        end  
    end

    % set time values for output
    t=dt*(0:S.nT);
    save(sprintf('%s/time.mat',dir),'t');

    % --- choose output dir (avoid shadowing 'dir') ---
    outDir  = fullfile('..','output','results');
    if ~isfolder(outDir), mkdir(outDir); end
    
    % --- pack useful metadata & signals and save (v7.3 handles big arrays) ---
    save(fullfile(outDir,'trace_data.mat'), ...
         't','usoma','rec_u','rec_h','rec_m','rec_n','-v7.3');

    makematlabmov('../output/results',filename,'../output/testvideo');

    % plottimeseries(t,rec_u,3);
end

function f = make_gate_rhs(aFun, bFun)
%MAKE_GATE_RHS Factory for gating RHS: f(v,x) = a(v)*(1-x) - b(v)*x
%   aFun, bFun are function handles like @an, @bn, @am, @bm, @ah, @bh.
    f = @(v,x) aFun(v).*(1 - x) - bFun(v).*x;
end

% this is the RK4 time stepping scheme set LeVeque 2007 for details
function pout = RK4(v,pp,dt,fun)
    p1 = pp;
    p2 = pp + (0.5).*dt.*fun(v,p1);
    p3 = pp + (0.5).*dt.*fun(v,p2);
    p4 = pp + (1.0).*dt.*fun(v,p3);
    pout = pp + (1/6).*dt.*(fun(v,p1)+2.*fun(v,p2)+2.*fun(v,p3)+fun(v,p4));
end

% I did a separate RK4 function for the reaction term, it involves for
% inputs.
% v = input voltage vector
% mm = input state m for potassium channels
% hh = input state h for potassium channels
% n = input state n for sodium channels
function rout = RK4_react(v,mm,hh,nn,dt,fun,P)
    r1 = v;
    r2 = v + (0.5).*dt.*fun(r1,mm,hh,nn,P);
    r3 = v + (0.5).*dt.*fun(r2,mm,hh,nn,P);
    r4 = v + (1.0).*dt.*fun(r3,mm,hh,nn,P);
    rout = v+(1/6).*dt.*(fun(r1,mm,hh,nn,P)+2.*fun(r2,mm,hh,nn,P)...
                        +2.*fun(r3,mm,hh,nn,P)+fun(r4,mm,hh,nn,P));
end

% this is the reaction term in the HH model
% v  = input voltage vector
% mm = input state m for potassium channels
% hh = input state h for potassium channels
% nn = input state n for sodium channels
function rout = react(v,mm,hh,nn,P)
    rout = (-1/P.C).*(P.gk.*nn.^4.*(v-P.ek)+P.gna.*mm.^3.*hh.*(v-P.ena)+P.gl.*(v-P.el));
end

% These functions are for the gating variables in the Hodgkin-Huxley
% formulism. carefully notice that for this simulation I am using MKS, 
% these gatinG functionS use [mV] and output [ms]^-1 so they need to be 
% properly scaled!!
function out=an(vin)
    vin = vin.*1e3;
    out=(-0.032).*(vin-15)./(exp(-1.*(vin-15)./5)-1);
    out = out*1e3;
end

function out=bn(vin)
    vin = vin.*1e3;
    out=(0.5).*exp(-1.*(vin-10)./40);
    out = out*1e3;
end

function out=am(vin)
    vin = vin.*1e3;
    out=(-0.32).*(vin-13)./(exp(-1.*(vin-13)./4)-1);
    out = out*1e3;
end

function out=bm(vin)
    vin = vin.*1e3;
    out=(0.28).*(vin-40)./(exp((vin-40)./5)-1);
    out = out*1e3;
end

function out=ah(vin)
    vin = vin.*1e3;
    out=(0.128).*exp(-1.*(vin-17)./18);
    out = out*1e3;
end

function out=bh(vin)
    vin = vin.*1e3;
    out=4./(exp((40-vin)./5)+1);
    out = out*1e3;
end