%% Model specification

cd /projectnb/crc-nak/sherfey/projects/ACC_simulations/ACC_models

Ne=80; % 4,80 (10,20,40,80,100)
Ni=.25*Ne;

% load model parameters fit to single cell intrinsic property data
ipfile='/projectnb/crc-nak/sherfey/projects/ACC_simulations/rat-ACd-results/rat-ACd-model_cell-intrinsic-properties.mat';
load(ipfile,'simdata','expdata');

% collect values for parameters to make heterogeneous in E-cell population
het={'gAHP','gcan','gcat','gh','gkca','gks','gkdr','gnaf','gpas'};
params=[];
for i=1:numel(het)
  params=cat(2,params,simdata.(het{i}));
end
% parameters of the multivariate normal distribution
MU = median(params,1);
SIGMA = cov(params);

% set baseline E-cell parameter values (median) for homogeneous network
parameters={};
for i=1:length(simdata.parameters)
  param=simdata.parameters{i};
  value=median(simdata.(param)); % mean(simdata.(param));
  parameters{end+1}=param;
  parameters{end+1}=value*ones(1,Ne); % homogeneous across population
end
% account for default-value discrepancies b/w implementations (my original Papoutsi [PP13] vs Durstewitz [DS02])
aAHP_scale=1; % 1, 1e6
CAF=600*4.6452; % 600, 600*4.6452
cainf = 5e-5; % .05, 5e-5
parameters{end+1}='aAHP_scale';
parameters{end+1}=aAHP_scale;
parameters{end+1}='CAF';
parameters{end+1}=CAF;
parameters{end+1}='cainf';
parameters{end+1}=cainf;
parameters{end+1}='IC_noise';
% set initial condition noise level (IC_noise*rand offsets each cell's gating variable IC)
IC_noise=.25;
parameters{end+1}=IC_noise;

% E-cell mechanisms
mechanisms={'AHP','cadyn','can','cat','h','iks','kdr','kca','naf','nap','pas'};

% generic inputs and state equations
assembly_size=Ne/2;
Ke1=zeros(1,Ne); Ke1(1:assembly_size)=1;    % input kernel for first E-cell assembly
Ke2=zeros(1,Ne); Ke2(assembly_size+1:Ne)=1; % input kernel for all other E-cell assemblies
E_input_def={'input(V)=-gAMPA.*(s1(k,:)+s2(k,:)).*(X-EAMPA); monitor input; EAMPA=0; gAMPA=0; onset=50; offset=inf;';
     sprintf('s1=getPoissonGating(baseline/2,dcAMPA1,acAMPA1,fAMPA1,phiAMPA1,onset,offset,tauAMPA,T,Npop,%s,kick,ramp_dc_flag1,ramp_ac_flag1);',toString(Ke1));
     sprintf('s2=getPoissonGating(baseline/2,dcAMPA2,acAMPA2,fAMPA2,phiAMPA2,onset,offset,tauAMPA,T,Npop,%s,kick,ramp_dc_flag2,ramp_ac_flag2);',toString(Ke2));
             'baseline=0; dcAMPA1=0; acAMPA1=0; fAMPA1=0; phiAMPA1=0; tauAMPA=2; kick=1; ramp_dc_flag1=0; ramp_ac_flag1=0; dcAMPA2=0; acAMPA2=0; fAMPA2=0; phiAMPA2=0; ramp_dc_flag2=0; ramp_ac_flag2=0;';
            };
I_input_def={'input(V)=-gAMPA.*s1(k,:).*(X-EAMPA); monitor input; EAMPA=0; gAMPA=0; onset=50; offset=inf;';
             's1=getPoissonGating(baseline,dcAMPA1,acAMPA1,fAMPA1,phiAMPA1,onset,offset,tauAMPA,T,Npop,ones(1,Npop),kick,ramp_dc_flag1,ramp_ac_flag1);';
             'baseline=0; dcAMPA1=0; acAMPA1=0; fAMPA1=0; phiAMPA1=0; tauAMPA=2; kick=1; ramp_dc_flag1=0; ramp_ac_flag1=0;';
            };
state_equations='dV/dt=(@current+input(V)+Iapp*(t>onset&t<offset))./Cm; Cm=1; Iapp=0; V(0)=-65;';

tauAMPA=2; 
tauNMDA=95;
tauGABA=13; % ms, inhibition decay time constant

% input parameters
onset=100; % 500 ms, stimulus start time
offset=inf; % ms, stimulus stop time
gAMPAe=.001;
enoise=7500; % Hz
dcAMPA1e=0;     dcAMPA2e=dcAMPA1e;
acAMPA1e=7500;  acAMPA2e=acAMPA1e;
fAMPA1e=0;      fAMPA2e=0;
E_input_parameters={'dcAMPA1',dcAMPA1e,'dcAMPA2',dcAMPA2e,'acAMPA1',acAMPA1e,'acAMPA2',acAMPA2e,'fAMPA1',fAMPA1e,'fAMPA2',fAMPA2e,'gAMPA',gAMPAe,'baseline',enoise,'tauAMPA',tauAMPA,'onset',onset,'offset',offset};
gAMPAi=0;
inoise=0;
dcAMPA1i=0;
acAMPA1i=0;
fAMPA1i=0;
I_input_parameters={'dcAMPA1',dcAMPA1i,'acAMPA1',acAMPA1i,'fAMPA1',fAMPA1i,'gAMPA',gAMPAi,'baseline',inoise,'tauAMPA',tauAMPA,'onset',onset,'offset',offset};

% connectivity parameters and kernels
gAMPAee=0; % E->E (<=.1)
gNMDAee=0; % E->E (.1-.2)
gAMPAie=.1; % I->E (.1)
gAMPAii=1; % I->I
if Ne==80
  gAMPAei=1; % E->I (.1)
else
  gAMPAei=.1; % E->I (.1)
end
Kii=ones(Ni)-eye(Ni);
Kei=ones(Ne,Ni);
Kie=ones(Ni,Ne);

% define assemblies for competition sims:
bwblockee=0; % kernal scaling for connection b/w E-assemblies (i.e., E1->E2 connected with gsyn=gee*bwblockee)
Kee=zeros(Ne,Ne); % E->E, [N_pre x N_post]
% connect E-cells within block
block=ones(assembly_size)-eye(assembly_size);
for i=1:(Ne/assembly_size) % loop over assemblies
  ind=(i-1)*assembly_size+(1:assembly_size);
  Kee(ind,ind)=block;
end
% connect E-cells between blocks
K12=(1-Kee-eye(Ne)); % connections b/w blocks 1 and 2
Kee=Kee+bwblockee*K12;
%Kee=ones(Ne)-eye(Ne);

% normalize kernels by number of presynaptic connections
Kee=Kee./repmat(max(1,sum(Kee,1)),[size(Kee,1) 1]);
Kei=Kei./repmat(max(1,sum(Kei,1)),[size(Kei,1) 1]);
Kie=Kie./repmat(max(1,sum(Kie,1)),[size(Kie,1) 1]);
Kii=Kii./repmat(max(1,sum(Kii,1)),[size(Kii,1) 1]);

% NETWORK SPECIFICATION
spec=[];
% Populations:
% E-cells (heterogeneous parameters)
spec.populations(1).name='E';
spec.populations(1).size=Ne;
spec.populations(1).equations=[state_equations E_input_def{:}];
spec.populations(1).mechanism_list=mechanisms;
spec.populations(1).parameters={parameters{:},E_input_parameters{:}};
% Wang-Buzsaki interneuron model, 1996
spec.populations(2).name='I';
spec.populations(2).size=Ni;
spec.populations(2).equations=[state_equations I_input_def{:}];
spec.populations(2).mechanism_list={'WB96FSiNa','WB96FSiK','WB96FSileak'};
spec.populations(2).parameters={'Cm',1,'Eleak',-65,'gleak',.1,'gNa',35,'gK',9,I_input_parameters{:}};    
% Connections:
% E->E (AMPA,NMDA): 50% assemblies
spec.connections(1).direction='E->E';
spec.connections(1).mechanism_list={'iAMPA','iNMDA'};
spec.connections(1).parameters={'gAMPA',gAMPAee,'netcon',Kee,'tauAMPA',tauAMPA,'tauNMDA',tauNMDA};
% E->I (AMPA,NMDA): all-to-all
spec.connections(2).direction='E->I';
spec.connections(2).mechanism_list={'iAMPA'};
spec.connections(2).parameters={'gAMPA',gAMPAei,'netcon',Kei,'tauAMPA',tauAMPA};
% I->E (GABA): all-to-all
spec.connections(3).direction='I->E';
spec.connections(3).mechanism_list={'iGABA'};
spec.connections(3).parameters={'gGABA',gAMPAie,'netcon',Kie,'tauGABA',tauGABA};
% I->I (GABA): all-to-all
spec.connections(4).direction='I->I';
spec.connections(4).mechanism_list={'iGABA'};
spec.connections(4).parameters={'gGABA',gAMPAii,'netcon',Kii,'tauGABA',tauGABA};

% store baseline model
base=spec;

% limit the number of cores to <= scc limit
maxNumCompThreads(4);

% define random seeds for all simulations (use same for hom and het)
num_seeds=1e3; % <= max_num_realizations
seedfile=sprintf('%gseeds.mat',num_seeds);
if exist(seedfile,'file')
  load(seedfile,'seeds');
else
  seeds=zeros(num_seeds,1);
  for i=1:num_seeds
    rng('shuffle');
    seeds(i)=getfield(rng,'Seed');
    pause(.01);
  end
  % save seeds to file (for reproducibility)
  save(seedfile,'seeds');
end
%fprintf('%i\n',seeds)
% ##############################################

if 1
  
  addpath /projectnb/crc-nak/sherfey/projects/ACC_simulations/ACC_models % add path to getEE()

  f1=[0:7.5:60]; f2=[0:7.5:60]; acAMPA1e=4500; enoise=acAMPA1e; 
  tauGABA=5; % 5, 13
  gAMPAee=[.05]; % E->E (<=.1)
  gNMDAee=[.2]; % (later: 0,.1,.2)
  gAMPAei=1;
  betweenblockee=[0]; % 0,1 kernal scaling for connection b/w E-assemblies (i.e., E1->E2 connected with gsyn=gee*bwblockee)
  withinblockee=1;
  vary={'E->E','gw',withinblockee;'E->E','gb',betweenblockee;'E','fAMPA1',f1;'E','fAMPA2',f2;'E->E','gNMDA',gNMDAee;'E->E','gAMPA',gAMPAee;'E','baseline',enoise;'E','(acAMPA1,acAMPA2)',acAMPA1e;'E','(dcAMPA1,dcAMPA2)',dcAMPA1e;'(I->E,I->I)','tauGABA',tauGABA;'E->I','gAMPA',gAMPAei};

  rep=2; % 1:10
  for rep=1:5

    % simulator options
    tspan=[0 5000]; dt=.01; solver='rk1'; compile_flag=1; downsample_factor=10;
    save_data_flag=0; cluster_flag=1; memory_limit='32G'; sims_per_job=10;
    analysis_functions=@CalcSpikeSync;
    analysis_options={'ROI_pairs',{'E_V',[0 .5],'E_V',[.5 1]}};
    plot_functions={@PlotData,@PlotData};
    plot_options={{'plot_type','rastergram'},{'plot_type','power','xlim',[0 100]}};
    simulator_options={'save_data_flag',save_data_flag,'tspan',tspan,'solver',solver,'dt',dt,'compile_flag',compile_flag,'verbose_flag',1,'cluster_flag',cluster_flag,'plot_functions',plot_functions,'plot_options',plot_options,'memory_limit',memory_limit,'downsample_factor',downsample_factor,'analysis_functions',analysis_functions,'analysis_options',analysis_options,'sims_per_job',sims_per_job};

    spec=base;
    spec.connections(1).mechanism_list={'iAMPAee','iNMDAee'};
    spec.connections(1).parameters={'tauAMPA',tauAMPA,'tauNMDA',tauNMDA};

    study_dir_root=sprintf('studies/sweeps_het-vs-hom_ac%g_EEwithin%g-%g_EEbetween%g-%g_gEEnmda%g-%g_gEEampa%g-%g_gEI%g_%g-%gms',acAMPA1e,withinblockee(1),withinblockee(end),betweenblockee(1),betweenblockee(end),gNMDAee(1),gNMDAee(end),gAMPAee(1),gAMPAee(end),gAMPAei,tspan);
    hom_study_dir=sprintf('%sNe%g_f_%g-%gHz-vs-%g-%gHz_tauI%gms_HOM_%g',study_dir_root,Ne,f1(1),f1(end),f2(1),f2(end),tauGABA,rep);
    [data,hom_studyinfo]=SimulateModel(spec,'random_seed',seeds(rep),'study_dir',hom_study_dir,'vary',vary,simulator_options{:});

    % het sims
    hetdegree=.2; % .2
    rng(seeds(rep)); % set the random seed for drawing parameter values
    values = mvnrnd(MU,hetdegree*SIGMA,Ne/2); % parameter values for one realization of E-cell population
    values = cat(1,values,values); % copy for each assembly
    values = max(0,values);       % turn off currents w/ negative max conductance
    % update model specification with heterogeneous parameter values
    mods={};
    for i=1:numel(het)
      mods=cat(1,mods,{'E',het{i},values(:,i)'});
    end
    spec=ApplyModifications(spec,mods);
    study_dir_root=sprintf('studies/sweeps_het-vs-hom_ac%g_EEwithin%g-%g_EEbetween%g-%g_gEEnmda%g-%g_gEEampa%g-%g_gEI%g_%g-%gms',acAMPA1e,withinblockee(1),withinblockee(end),betweenblockee(1),betweenblockee(end),gNMDAee(1),gNMDAee(end),gAMPAee(1),gAMPAee(end),gAMPAei,tspan);
    het_study_dir=sprintf('%sNe%g_f_%g-%gHz-vs-%g-%gHz_tauI%gms_HET%g_%g',study_dir_root,Ne,f1(1),f1(end),f2(1),f2(end),tauGABA,hetdegree,rep);
    [data,het_studyinfo]=SimulateModel(spec,'random_seed',seeds(rep),'study_dir',het_study_dir,'vary',vary,simulator_options{:});

    % het sims
    hetdegree=.6; % .2
    rng(seeds(rep)); % set the random seed for drawing parameter values
    values = mvnrnd(MU,hetdegree*SIGMA,Ne/2); % parameter values for one realization of E-cell population
    values = cat(1,values,values); % copy for each assembly
    values = max(0,values);       % turn off currents w/ negative max conductance
    % update model specification with heterogeneous parameter values
    mods={};
    for i=1:numel(het)
      mods=cat(1,mods,{'E',het{i},values(:,i)'});
    end
    spec=ApplyModifications(spec,mods);
    study_dir_root=sprintf('studies/sweeps_het-vs-hom_ac%g_EEwithin%g-%g_EEbetween%g-%g_gEEnmda%g-%g_gEEampa%g-%g_gEI%g_%g-%gms',acAMPA1e,withinblockee(1),withinblockee(end),betweenblockee(1),betweenblockee(end),gNMDAee(1),gNMDAee(end),gAMPAee(1),gAMPAee(end),gAMPAei,tspan);
    het_study_dir=sprintf('%sNe%g_f_%g-%gHz-vs-%g-%gHz_tauI%gms_HET%g_%g',study_dir_root,Ne,f1(1),f1(end),f2(1),f2(end),tauGABA,hetdegree,rep);
    [data,het_studyinfo]=SimulateModel(spec,'random_seed',seeds(rep),'study_dir',het_study_dir,'vary',vary,simulator_options{:});
    
    hetdegree=1; % .2
    rng(seeds(rep)); % set the random seed for drawing parameter values
    values = mvnrnd(MU,hetdegree*SIGMA,Ne/2); % parameter values for one realization of E-cell population
    values = cat(1,values,values); % copy for each assembly
    values = max(0,values);       % turn off currents w/ negative max conductance
    % update model specification with heterogeneous parameter values
    mods={};
    for i=1:numel(het)
      mods=cat(1,mods,{'E',het{i},values(:,i)'});
    end
    spec=ApplyModifications(spec,mods);
    study_dir_root=sprintf('studies/sweeps_het-vs-hom_ac%g_EEwithin%g-%g_EEbetween%g-%g_gEEnmda%g-%g_gEEampa%g-%g_gEI%g_%g-%gms',acAMPA1e,withinblockee(1),withinblockee(end),betweenblockee(1),betweenblockee(end),gNMDAee(1),gNMDAee(end),gAMPAee(1),gAMPAee(end),gAMPAei,tspan);
    het_study_dir=sprintf('%sNe%g_f_%g-%gHz-vs-%g-%gHz_tauI%gms_HET%g_%g',study_dir_root,Ne,f1(1),f1(end),f2(1),f2(end),tauGABA,hetdegree,rep);
    [data,het_studyinfo]=SimulateModel(spec,'random_seed',seeds(rep),'study_dir',het_study_dir,'vary',vary,simulator_options{:});
    
  end

  unix(sprintf('cat %s/pbsout/sim_job1.out',hom_studyinfo.simulations(1).batch_dir));
  
  
  
end



% -------------------------------------------------------------------------

addpath /projectnb/crc-nak/sherfey/projects/ACC_simulations/ACC_models % add path to getEE()

f1=[0:7.5:60]; f2=[0:7.5:60]; acAMPA1e=4500; enoise=acAMPA1e; 
tauGABA=5; % 5, 13
gAMPAee=[.05]; % E->E (<=.1)
gNMDAee=[.2]; % (later: 0,.1,.2)
gAMPAei=1;
betweenblockee=[0]; % 0,1 kernal scaling for connection b/w E-assemblies (i.e., E1->E2 connected with gsyn=gee*bwblockee)
withinblockee=1;
vary={'E->E','gw',withinblockee;'E->E','gb',betweenblockee;'E','fAMPA1',f1;'E','fAMPA2',f2;'E->E','gNMDA',gNMDAee;'E->E','gAMPA',gAMPAee;'E','baseline',enoise;'E','(acAMPA1,acAMPA2)',acAMPA1e;'E','(dcAMPA1,dcAMPA2)',dcAMPA1e;'(I->E,I->I)','tauGABA',tauGABA;'E->I','gAMPA',gAMPAei};

rep=2; % 1:10
% for rep=3:10
  
  % simulator options
  tspan=[0 5000]; dt=.01; solver='rk1'; compile_flag=1; downsample_factor=10;
  save_data_flag=1; cluster_flag=1; memory_limit='32G'; sims_per_job=1;
  analysis_functions=@CalcSpikeSync;
  analysis_options={'ROI_pairs',{'E_V',[0 .5],'E_V',[.5 1]}};
  plot_functions={@PlotData,@PlotData};
  plot_options={{'plot_type','rastergram'},{'plot_type','power','xlim',[0 100]}};
  simulator_options={'save_data_flag',save_data_flag,'tspan',tspan,'solver',solver,'dt',dt,'compile_flag',compile_flag,'verbose_flag',1,'cluster_flag',cluster_flag,'plot_functions',plot_functions,'plot_options',plot_options,'memory_limit',memory_limit,'downsample_factor',downsample_factor,'analysis_functions',analysis_functions,'analysis_options',analysis_options,'sims_per_job',sims_per_job};

  spec=base;
  spec.connections(1).mechanism_list={'iAMPAee','iNMDAee'};
  spec.connections(1).parameters={'tauAMPA',tauAMPA,'tauNMDA',tauNMDA};

  study_dir_root=sprintf('studies/sweeps_het-vs-hom_ac%g_EEwithin%g-%g_EEbetween%g-%g_gEEnmda%g-%g_gEEampa%g-%g_gEI%g_%g-%gms',acAMPA1e,withinblockee(1),withinblockee(end),betweenblockee(1),betweenblockee(end),gNMDAee(1),gNMDAee(end),gAMPAee(1),gAMPAee(end),gAMPAei,tspan);
  hom_study_dir=sprintf('%sNe%g_f_%g-%gHz-vs-%g-%gHz_tauI%gms_HOM_%g',study_dir_root,Ne,f1(1),f1(end),f2(1),f2(end),tauGABA,rep);
  [data,hom_studyinfo]=SimulateModel(spec,'random_seed',seeds(rep),'study_dir',hom_study_dir,'vary',vary,simulator_options{:});

  % het sims
  hetdegree=.05; % .2
  rng(seeds(rep)); % set the random seed for drawing parameter values
  values = mvnrnd(MU,hetdegree*SIGMA,Ne/2); % parameter values for one realization of E-cell population
  values = cat(1,values,values); % copy for each assembly
  values = max(0,values);       % turn off currents w/ negative max conductance
  % update model specification with heterogeneous parameter values
  mods={};
  for i=1:numel(het)
    mods=cat(1,mods,{'E',het{i},values(:,i)'});
  end
  spec=ApplyModifications(spec,mods);
  study_dir_root=sprintf('studies/sweeps_het-vs-hom_ac%g_EEwithin%g-%g_EEbetween%g-%g_gEEnmda%g-%g_gEEampa%g-%g_gEI%g_%g-%gms',acAMPA1e,withinblockee(1),withinblockee(end),betweenblockee(1),betweenblockee(end),gNMDAee(1),gNMDAee(end),gAMPAee(1),gAMPAee(end),gAMPAei,tspan);
  het_study_dir=sprintf('%sNe%g_f_%g-%gHz-vs-%g-%gHz_tauI%gms_HET%g_%g',study_dir_root,Ne,f1(1),f1(end),f2(1),f2(end),tauGABA,hetdegree,rep);
  [data,het_studyinfo]=SimulateModel(spec,'random_seed',seeds(rep),'study_dir',het_study_dir,'vary',vary,simulator_options{:});

% end

unix(sprintf('cat %s/pbsout/sim_job1.out',hom_studyinfo.simulations(1).batch_dir));

% analysis
maxlag_time=10; clear homstats1 homstats2 homstats3 hetstats1 hetstats2 hetstats3

homdata=ImportData(hom_study_dir);
data=SelectData(homdata,'time_limits',[500 2500]);
for i=1:length(data)
  homstats1(i)=AnalyzeData(data(i),@CalcSpikeSync,'ROI_pairs',{'E_V',[0 .5],'E_V',[.5 1]},'maxlag_time',maxlag_time);
end
data=SelectData(homdata,'time_limits',[3000 5000]);
for i=1:length(data)
  homstats2(i)=AnalyzeData(data(i),@CalcSpikeSync,'ROI_pairs',{'E_V',[0 .5],'E_V',[.5 1]},'maxlag_time',maxlag_time);
end
data=SelectData(homdata,'time_limits',[500 5000]);
for i=1:length(data)
  homstats3(i)=AnalyzeData(data(i),@CalcSpikeSync,'ROI_pairs',{'E_V',[0 .5],'E_V',[.5 1]},'maxlag_time',maxlag_time);
end

hetdata=ImportData(het_study_dir);
data=SelectData(hetdata,'time_limits',[500 2500]);
for i=1:length(data)
  hetstats1(i)=AnalyzeData(data(i),@CalcSpikeSync,'ROI_pairs',{'E_V',[0 .5],'E_V',[.5 1]},'maxlag_time',maxlag_time);
end
data=SelectData(hetdata,'time_limits',[3000 5000]);
for i=1:length(data)
  hetstats2(i)=AnalyzeData(data(i),@CalcSpikeSync,'ROI_pairs',{'E_V',[0 .5],'E_V',[.5 1]},'maxlag_time',maxlag_time);
end
data=SelectData(hetdata,'time_limits',[500 5000]);
for i=1:length(data)
  hetstats3(i)=AnalyzeData(data(i),@CalcSpikeSync,'ROI_pairs',{'E_V',[0 .5],'E_V',[.5 1]},'maxlag_time',maxlag_time);
end

PlotData(homdata,'plot_type','rastergram');
PlotData(homdata,'plot_type','rastergram','xlim',[500 2500]);
PlotData(homdata,'plot_type','rastergram','xlim',[3000 5000]);

% (uncoupled,hom):
% 4x imagesc(f1,f2): competition and integration for (500-2500ms) & (3000-5000ms)

stats=homstats1; name='hom(500-2500ms)';
% stats=homstats2; name='hom(3000-5000ms)';
stats=homstats3; name='hom(500-5000ms)';
% stats=hetstats1; name='het(500-2500ms)';
% stats=hetstats2; name='het(3000-5000ms)';
% stats=hetstats3; name='het(500-5000ms)';

f1=[stats.E_fAMPA1]; uf1=unique(f1); nf1=length(uf1);
f2=[stats.E_fAMPA2]; uf2=unique(f2); nf2=length(uf2);
n1=zeros(nf1,nf2);
n2=zeros(nf1,nf2);
c12=zeros(nf1,nf2);
for k=1:length(stats)
  i=(f1(k)==uf1);
  j=(f2(k)==uf2);
  n1(i,j)=stats(k).pairs.Nspikes1;
  n2(i,j)=stats(k).pairs.Nspikes2;
  c12(i,j)=stats(k).pairs.xcsum_pops;
end
figure('position',[670 200 560 720]);
subplot(2,1,1); imagesc(uf1,uf2,abs(n1-n2)); title(['|dn| ' name]); axis square xy; colorbar;
%subplot(2,1,1); imagesc(uf1,uf2,abs((n1-n2)./(n1+n2))); title(['|dn*| ' name]); axis square xy; colorbar; caxis([0 1])
subplot(2,1,2); imagesc(uf1,uf2,c12); title(['spkcoh ' name]); axis square xy; colorbar; caxis([0 1])

% resonance frequency
dn=abs(n1-n2); [i,j]=find(dn==min(dn(:))); [uf1(i) uf2(j)]
fr=15;

% stats=ImportResults(hom_study_dir,@CalcSpikeSync);
dnrel=arrayfun(@(x)x.pairs.dNsumN,stats);
c12=arrayfun(@(x)x.pairs.xcsum_pops,stats);
n1=arrayfun(@(x)x.pairs.Nspikes1,stats);
n2=arrayfun(@(x)x.pairs.Nspikes2,stats);
dn=n1-n2;
df=f1-f2;
dfr=abs(f1-fr)-abs(f2-fr);
f0i=(f1==0)|(f2==0);
figure('position',[510 210 740 660]);
subplot(2,2,1); plot(dfr(~f0i),abs(dn(~f0i)),'o',dfr(f0i),abs(dn(f0i)),'x'); xlabel('dfr'); ylabel('|dn|'); line(xlim,[0 0]); line([0 0],ylim);
subplot(2,2,2); plot(df(~f0i),abs(dn(~f0i)),'o',df(f0i),abs(dn(f0i)),'x'); xlabel('df'); ylabel('|dn|'); line(xlim,[0 0]); line([0 0],ylim);
subplot(2,2,3); plot(dfr(~f0i),c12(~f0i),'o',dfr(f0i),c12(f0i),'x'); xlabel('dfr'); ylabel('c12'); ylim([0 1]); line([0 0],ylim);
subplot(2,2,4); plot(df(~f0i),c12(~f0i),'o',df(f0i),c12(f0i),'x'); xlabel('df'); ylabel('c12'); ylim([0 1]); line([0 0],ylim); 
title(name);

figure;
plot(abs(dfr),abs(dnrel),'o'); xlabel('|dfr|'); ylabel('|dnrel|');

% save data with only input and voltage
keepfields={'E_V','time'};
rmfields=setdiff(data(1).labels,keepfields);
data=rmfield(homdata,rmfields);
[data.labels]=deal(keepfields);
%PlotData(data,'plot_type','rastergram');
%stats=AnalyzeData(data(1),@CalcSpikeSync,'ROI_pairs',{'E_V',[0 .5],'E_V',[.5 1]},'maxlag_time',10);

% 4x reduction: remove E_input and set downsample_factor=20
dat1=CalcFR(data);
for i=1:length(data)
  data(i).time=data(i).time(1:2:end);
  data(i).E_V=data(i).E_V(1:2:end,:);
end
dat2=CalcFR(data); % raster on downsample_factor=20
for i=1:length(data)
  data(i).time=data(i).time(1:2:end);
  data(i).E_V=data(i).E_V(1:2:end,:);
end
dat3=CalcFR(data); % raster on downsample_factor=40

% then: redo plots for (uncoupled,het), (coupled,hom), and (coupled,het)

return

data=SelectData(homdata,'time_limits',[500 2500]);
homdN=[]; homco=[]; homnco=[]; homsync=[]; clear homstats
for i=1:length(data)
  homstats(i)=AnalyzeData(data(i),@CalcSpikeSync,'ROI_pairs',{'E_V',[0 .5],'E_V',[.5 1]},'maxlag_time',maxlag_time);
  homdN=[homdN homstats(i).pairs.dNsumN]; homsync=[homsync homstats(i).pairs.xcsum_pops];
  homco=[homco homstats(i).pairs.coactivity]; homnco=[homnco homstats(i).pairs.ncoactivity];
end
[homdN;homsync;homco;homnco]
data=SelectData(hetdata,'time_limits',[1000 3000]);
hetdN=[]; hetco=[]; hetnco=[]; hetsync=[]; clear hetstats
for i=1:length(data)
  hetstats(i)=AnalyzeData(data(i),@CalcSpikeSync,'ROI_pairs',{'E_V',[0 .5],'E_V',[.5 1]},'maxlag_time',maxlag_time);
  hetdN=[hetdN hetstats(i).pairs.dNsumN]; hetsync=[hetsync hetstats(i).pairs.xcsum_pops];
  hetco=[hetco hetstats(i).pairs.coactivity]; hetnco=[hetnco hetstats(i).pairs.ncoactivity];
end
[hetdN;hetsync;hetco;hetnco]

figure;
subplot(2,1,1); plot(betweenblockee,homsync,'bo-',betweenblockee,hetsync,'ro-'); ylabel('sync'); legend('hom','het');
subplot(2,1,2); plot(betweenblockee,abs(homdN),'bo-',betweenblockee,abs(hetdN),'ro-'); ylabel('|dN|'); legend('hom','het');
xlabel('betweenblockee');

xlims=[1000 3000];
PlotData(homdata,'plot_type','rastergram','xlim',xlims);
PlotData(hetdata,'plot_type','rastergram','xlim',xlims);
