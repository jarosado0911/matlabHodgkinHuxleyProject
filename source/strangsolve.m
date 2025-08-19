function strangsolve(dt,clamp_index,rec_ind,force,filename,outputFolder)

% read radius and subset names from readSWC function
[~,~,~,~,a,~]=readswc(filename);

% need to scale radii to MICRO METERS
a = a*1e-6;
n = length(a);

% this is for make the output folder
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

% Matrix for solving SBDF Method
LHS = Id-dt.*A;
RHS = Id;

if ~force
    LHS_clamp = LHS;
    LHS_clamp(clamp_index,:) =0;
    LHS_clamp(clamp_index,clamp_index) = 1;
end

dLHS = decomposition(LHS);

% this preallocates a matrix for saving recorded voltage
u_rec = zeros(length((0:S.nT)),length(rec_ind));
n_rec = zeros(length((0:S.nT)),length(rec_ind));
m_rec = zeros(length((0:S.nT)),length(rec_ind));
h_rec = zeros(length((0:S.nT)),length(rec_ind));

% wrappers that match (u,s) -> F(s, a(u), b(u))
Fn = @(u,s) gates.F(s, gates.an(u), gates.bn(u));
Fm = @(u,s) gates.F(s, gates.am(u), gates.bm(u));
Fh = @(u,s) gates.F(s, gates.ah(u), gates.bh(u));

% reaction
rX = @(v, ss, P) gates.X(v, ss, P);

for i=0:S.nT
    % set soma to clamp
    if ((i*dt >= S.delay) && (i*dt <= S.stop))
        u(clamp_index)=S.vClamp;
    end

    % The scheme is a Strang splitting
    % (1) Do a 1/2 time step with the ODEs and Reaction part of OP-split
    nn = timestep.rk4(u, nn, dt/2, Fn);
    mm = timestep.rk4(u, mm, dt/2, Fm);
    hh = timestep.rk4(u, hh, dt/2, Fh);
    u  = timestep.rk4(u, [mm,nn,hh],dt/2,rX,P);
    
    % (2) Do full step of diffusion solve, this is part of OP-split
    if ((i*dt >= S.delay) && (i*dt <= S.stop) && ~force)
        u = LHS_clamp\(RHS*u);
    else
        u = LHS\(RHS*u);
    end

    % (3) Finish off with other half to time step from ODEs and Reaction
    % part
    u  = timestep.rk4(u, [mm,nn,hh],dt/2,rX,P);
    nn = timestep.rk4(u, nn, dt/2, Fn);
    mm = timestep.rk4(u, mm, dt/2, Fm);
    hh = timestep.rk4(u, hh, dt/2, Fh);

    % set soma to clamp
    if ((i*dt >= S.delay) && (i*dt <= S.stop))
        u(clamp_index)=S.vClamp;
    end  

    % records the voltage
    for ii=1:length(rec_ind)
        u_rec(i+1,ii) = u(rec_ind(ii));  n_rec(i+1,ii) = nn(rec_ind(ii));
        m_rec(i+1,ii) = mm(rec_ind(ii)); h_rec(i+1,ii) = hh(rec_ind(ii));
    end
    
    fprintf('t= %f [s]\n',i*dt)
end

% set time values for output
t=dt*(0:S.nT);

for i=1:length(rec_ind)
    dispname = sprintf('Index = %i',rec_ind(i));
    plot(t.*1e3,u_rec(:,i)*.1e3,'DisplayName', dispname);
end

% set the figure titles
title('Voltage profiles')
xlabel('time [ms]')
ylabel('voltage [mV]')
legend('show')

% save the output recorded voltage
save(sprintf('%s/u_rec.mat',outputFolder),'u_rec');
save(sprintf('%s/n_rec.mat',outputFolder),'n_rec');
save(sprintf('%s/m_rec.mat',outputFolder),'m_rec');
save(sprintf('%s/h_rec.mat',outputFolder),'h_rec');

% save soma voltage and time voltage as .mat files
save(sprintf('%s/time.mat',dir),'t')
end