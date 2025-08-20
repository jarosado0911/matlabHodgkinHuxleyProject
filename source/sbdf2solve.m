function sbdf2solve(dt,clamp_index,rec_ind,filename,outputFolder,pname,saveall)

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
npre = nn; mpre = mm; 
hpre = hh; upre = u;

% Make sparse stencil matrices
A = stencilmaker(n,P.R,a,P.C,filename);
Id = eye(n);

% Matrix for solving SBDF Method
LHS = Id-(2/3).*dt.*A;
RHS = Id;

dLHS = decomposition(LHS);

% this preallocates a matrix for saving recorded voltage
%u_rec = zeros(length((0:S.nT)),length(rec_ind));
%n_rec = zeros(length((0:S.nT)),length(rec_ind));
%m_rec = zeros(length((0:S.nT)),length(rec_ind));
%h_rec = zeros(length((0:S.nT)),length(rec_ind));

 % initialize empty recording vectors
 usoma=[]; rec_u =[]; rec_h =[]; rec_m =[]; rec_n =[]; 

% wrappers that match (u,s) -> F(s, a(u), b(u))
Fn = @(u,s) gates.F(s, gates.an(u), gates.bn(u));
Fm = @(u,s) gates.F(s, gates.am(u), gates.bm(u));
Fh = @(u,s) gates.F(s, gates.ah(u), gates.bh(u));

for i=0:S.nT

    usoma=[usoma;u(clamp_index)];
    rec_u = [rec_u, u(rec_ind)];
    rec_h = [rec_h, hh(rec_ind)];
    rec_m = [rec_m, mm(rec_ind)];
    rec_n = [rec_n, nn(rec_ind)];

    % update b vector via sbdf2 rule
    b = timestep.sbdf2b(P,dt, nn, mm, hh, u, npre, mpre, hpre, upre);

    % update states using SBDF2
    bN = timestep.sbdf2(dt,u,upre,nn,npre,Fn);
    bM = timestep.sbdf2(dt,u,upre,mm,mpre,Fm);
    bH = timestep.sbdf2(dt,u,upre,hh,hpre,Fh);
        
    % store current states in previous variable   
    upre = u; npre = nn; mpre = mm; hpre = hh;
    % update states
    nn = bN; mm = bM; hh = bH;
    
    if ((i*dt >= S.delay) && (i*dt <= S.stop))
        b(clamp_index) = S.vClamp;
        ej = u.*0; ej(clamp_index)=1;
        rj = LHS(clamp_index,:); rj(clamp_index)=rj(clamp_index)-1;
        
        z = dLHS\(ej);
        y = dLHS\b;
        
        u = y + ((rj*y)/(1-rj*z))*z;
    else
        u = dLHS\(RHS*b);
    end
    
    if saveall == 1
        writematrix(u,sprintf('%s/data/vm_t%i.dat',dir,i));
    end

    % records the voltage
    %for ii=1:length(rec_ind)
    %    u_rec(i+1,ii) = u(rec_ind(ii));  n_rec(i+1,ii) = nn(rec_ind(ii));
    %    m_rec(i+1,ii) = mm(rec_ind(ii)); h_rec(i+1,ii) = hh(rec_ind(ii));
    %end
    
    fprintf('t= %f [s]\n',i*dt)
end

% set time values for output
t=dt*(0:S.nT);

for i=1:length(rec_ind)
    dispname = sprintf('%s',pname);
    plot(t.*1e3,rec_u(i,:)*.1e3,'DisplayName', dispname);
end

% set the figure titles
title('Voltage profiles')
xlabel('time [ms]')
ylabel('voltage [mV]')
legend('show')

% save the output recorded voltage
%save(sprintf('%s/u_rec.mat',outputFolder),'u_rec');
%save(sprintf('%s/n_rec.mat',outputFolder),'n_rec');
%save(sprintf('%s/m_rec.mat',outputFolder),'m_rec');
%save(sprintf('%s/h_rec.mat',outputFolder),'h_rec');
save(fullfile(outputFolder,'trace_data.mat'), ...
         't','rec_u','rec_n','rec_h','rec_m','rec_n','-v7.3');
% save soma voltage and time voltage as .mat files
save(sprintf('%s/time.mat',dir),'t')
end