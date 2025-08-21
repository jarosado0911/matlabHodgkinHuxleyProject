function strangsolve(method,dt,clamp_index,rec_ind,force,filename,outputFolder,pname,saveall)

% read radius and subset names from readSWC function
[~,~,~,~,a,~]=readswc(filename);

% need to scale radii to MICRO METERS
a = a*1e-6;
n = length(a);

% this is to make the output folder
dir = sprintf('%s',outputFolder);
mkdir(dir); mkdir(sprintf('%s/data',dir));

% Hodgkin Huxley Parameters
P = hh_params();

% Simulation Parameters
S = sim_params(dt);

u=ones(n,1).*S.vStart;  
nn=zeros(n,1); nn(:,1)=S.ni;
mm=zeros(n,1); mm(:,1)=S.mi;
hh=zeros(n,1); hh(:,1)=S.hi;

% Make sparse stencil matrices
A = stencilmaker(n,P.R,a,P.C,filename);
Id = eye(n);

% Matrix for solving diffusion problem
LHS = Id-dt.*A;
RHS = Id;

dLHS = decomposition(LHS);
dLHS_clamp = [];

% if not force, then we address the clamp condition as a dirichelet bndry
if ~force
    LHS_clamp = LHS;
    LHS_clamp(clamp_index,:) = 0;
    LHS_clamp(clamp_index,clamp_index) = 1;
    dLHS_clamp = decomposition(LHS_clamp);
end

% initialize empty recording vectors
usoma=[]; rec_u =[]; rec_h =[]; rec_m =[]; rec_n =[]; 

% wrappers that match (u,s) -> F(s, a(u), b(u))
Fn = @(u,s) gates.F(s, gates.an(u), gates.bn(u));
Fm = @(u,s) gates.F(s, gates.am(u), gates.bm(u));
Fh = @(u,s) gates.F(s, gates.ah(u), gates.bh(u));

% reaction
rX = @(v, ss, P) gates.X(v, ss, P);

% set ODE time solve
switch(method)
    case 'trbdf'
        ODEstep = @(u, s, dt, fun)    timestep.trbdf2(u, s, dt, fun);
        ODEstepX= @(u, s, dt, fun, P) timestep.trbdf2(u, s, dt, fun, P);
    case 'mid'
        ODEstep = @(u, s, dt, fun)    timestep.mid(u, s, dt, fun);
        ODEstepX= @(u, s, dt, fun, P) timestep.mid(u, s, dt, fun, P);
    case 'tr'
        ODEstep = @(u, s, dt, fun)    timestep.tr(u, s, dt, fun);
        ODEstepX= @(u, s, dt, fun, P) timestep.tr(u, s, dt, fun, P);
    case 'hn'
        ODEstep = @(u, s, dt, fun)    timestep.hn(u, s, dt, fun);
        ODEstepX= @(u, s, dt, fun, P) timestep.hn(u, s, dt, fun, P);
    case 'fe'
        ODEstep = @(u, s, dt, fun)    timestep.fe(u, s, dt, fun);
        ODEstepX= @(u, s, dt, fun, P) timestep.fe(u, s, dt, fun, P);
    case 'be'
        ODEstep = @(u, s, dt, fun)    timestep.be(u, s, dt, fun);
        ODEstepX= @(u, s, dt, fun, P) timestep.be(u, s, dt, fun, P);
    case 'rk4'
        ODEstep = @(u, s, dt, fun)    timestep.rk4(u, s, dt, fun);
        ODEstepX= @(u, s, dt, fun, P) timestep.rk4(u, s, dt, fun, P);
end

for i=0:S.nT
    % set soma to clamp
    if ((i*dt >= S.delay) && (i*dt <= S.stop))
        u(clamp_index)=S.vClamp;
    end

    usoma=[usoma;u(clamp_index)];
    rec_u = [rec_u, u(rec_ind)];
    rec_h = [rec_h, hh(rec_ind)];
    rec_m = [rec_m, mm(rec_ind)];
    rec_n = [rec_n, nn(rec_ind)];

    % The scheme is a Strang splitting
    % (1) Do a 1/2 time step with the ODEs and Reaction part of OP-split
    nn = ODEstep(u, nn, dt/2, Fn);
    mm = ODEstep(u, mm, dt/2, Fm);
    hh = ODEstep(u, hh, dt/2, Fh);
    u  = ODEstepX(u, [mm,nn,hh],dt/2,rX,P);
    
    % (2) Do full step of diffusion solve, this is part of OP-split
    if ((i*dt >= S.delay) && (i*dt <= S.stop) && ~force)
        u = dLHS_clamp\(RHS*u);
    else
        u = dLHS\(RHS*u);
    end

    % (3) Finish off with other half to time step from ODEs and Reaction
    % part
    u  = ODEstepX(u, [mm,nn,hh],dt/2,rX,P);
    nn = ODEstep(u, nn, dt/2, Fn);
    mm = ODEstep(u, mm, dt/2, Fm);
    hh = ODEstep(u, hh, dt/2, Fh);

    % set soma to clamp
    if ((i*dt >= S.delay) && (i*dt <= S.stop))
        u(clamp_index)=S.vClamp;
    end  

    if saveall == 1
        writematrix(u,sprintf('%s/data/vm_t%i.dat',dir,i));
    end

    fprintf('t= %f [s]\n',i*dt)
end

% set time values for output
t=dt*(0:S.nT);
ylim([-10,50])
for i=1:length(rec_ind)
    dispname = sprintf('%s',pname);
    plot(t.*1e3,rec_u(i,:).*1e3,'DisplayName', dispname);
end

% set the figure titles
title('Voltage profiles')
xlabel('time [ms]')
ylabel('voltage [mV]')
legend('show')

save(fullfile(outputFolder,'trace_data.mat'), ...
         't','rec_u','rec_n','rec_h','rec_m','rec_n','-v7.3');

% save soma voltage and time voltage as .mat files
save(sprintf('%s/time.mat',dir),'t')
end