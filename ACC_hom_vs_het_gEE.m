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

% next: [hom within only]; hom within & between; then het

% input parameters
f1=[20 42.5]; f2=f1; acAMPA1e=4500; enoise=acAMPA1e; 
f1=20; f2=42.5; acAMPA1e=4500; enoise=acAMPA1e; 
f1=33; f2=37; acAMPA1e=4500; enoise=acAMPA1e; 
f1=20; f2=24; acAMPA1e=4500; enoise=acAMPA1e; 
tauGABA=5; % 5, 13
gAMPAee=[0 .05]; % E->E (<=.1)
gNMDAee=[0 .1]; % (later: 0,.1,.2)
gAMPAei=1;
vary={'E','fAMPA1',f1;'E','fAMPA2',f2;'E->E','gNMDA',gNMDAee;'E->E','gAMPA',gAMPAee;'E','baseline',enoise;'E','(acAMPA1,acAMPA2)',acAMPA1e;'E','(dcAMPA1,dcAMPA2)',dcAMPA1e;'(I->E,I->I)','tauGABA',tauGABA;'E->I','gAMPA',gAMPAei};

% betweenblockee=0; % kernal scaling for connection b/w E-assemblies (i.e., E1->E2 connected with gsyn=gee*bwblockee)
% withinblockee=1;

betweenblockee=1; % kernal scaling for connection b/w E-assemblies (i.e., E1->E2 connected with gsyn=gee*bwblockee)
withinblockee=1;
% connect E-cells within and between blocks
K11=zeros(Ne,Ne); % E->E, [N_pre x N_post]
block=ones(assembly_size)-eye(assembly_size);
for i=1:(Ne/assembly_size) % loop over assemblies
  ind=(i-1)*assembly_size+(1:assembly_size);
  K11(ind,ind)=block;
end
K12=(1-K11-eye(Ne)); % connections b/w blocks 1 and 2
Kee=withinblockee*K11+betweenblockee*K12;
Kee=Kee./repmat(max(1,sum(Kee,1)),[size(Kee,1) 1]);
figure; imagesc(Kee)
spec=ApplyModifications(base,{'E->E','netcon',Kee});
  
% simulator options
tspan=[0 3000]; dt=.01; solver='rk1'; compile_flag=1; downsample_factor=10;
save_data_flag=1; cluster_flag=1; memory_limit='32G'; sims_per_job=1;
analysis_functions=@CalcSpikeSync;
analysis_options={'ROI_pairs',{'E_V',[0 .5],'E_V',[.5 1]}};
plot_functions={@PlotData,@PlotData};
plot_options={{'plot_type','rastergram'},{'plot_type','power','xlim',[0 100]}};
simulator_options={'save_data_flag',save_data_flag,'tspan',tspan,'solver',solver,'dt',dt,'compile_flag',compile_flag,'verbose_flag',1,'cluster_flag',cluster_flag,'plot_functions',plot_functions,'plot_options',plot_options,'memory_limit',memory_limit,'downsample_factor',downsample_factor,'analysis_functions',analysis_functions,'analysis_options',analysis_options,'sims_per_job',sims_per_job};

rep=1;
study_dir_root=sprintf('studies/sweeps_het-vs-hom_ac%g_EEbetween%g_EEwithin%g_gEEnmda%g-%g_gEEampa%g-%g_gEI%g',acAMPA1e,withinblockee,betweenblockee,gNMDAee(1),gNMDAee(end),gAMPAee(1),gAMPAee(end),gAMPAei);
hom_study_dir=sprintf('%sNe%g_f_%g-%gHz-vs-%g-%gHz_tauI%gms_HOM_%g',study_dir_root,Ne,f1(1),f1(end),f2(1),f2(end),tauGABA,rep);
[data,hom_studyinfo]=SimulateModel(spec,'random_seed',seeds(rep),'study_dir',hom_study_dir,'vary',vary,simulator_options{:});

unix(sprintf('cat %s/pbsout/sim_job1.out',hom_studyinfo.simulations(1).batch_dir));

homdata=ImportData(hom_study_dir);
PlotData(homdata,'plot_type','rastergram');

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
study_dir_root=sprintf('studies/sweeps_het-vs-hom_ac%g_EEbetween%g_EEwithin%g_gEEnmda%g-%g_gEEampa%g-%g_gEI%g',acAMPA1e,withinblockee,betweenblockee,gNMDAee(1),gNMDAee(end),gAMPAee(1),gAMPAee(end),gAMPAei);
het_study_dir=sprintf('%sNe%g_f_%g-%gHz-vs-%g-%gHz_tauI%gms_HET%g_%g',study_dir_root,Ne,f1(1),f1(end),f2(1),f2(end),tauGABA,hetdegree,rep);
[data,het_studyinfo]=SimulateModel(spec,'random_seed',seeds(rep),'study_dir',het_study_dir,'vary',vary,simulator_options{:});

unix(sprintf('cat %s/pbsout/sim_job1.out',het_studyinfo.simulations(1).batch_dir));

hetdata=ImportData(het_study_dir);
PlotData(hetdata,'plot_type','rastergram');

% analysis
maxlag_time=10;
data=SelectData(homdata,'time_limits',[1000 3000]);
homdN=[]; homco=[]; homnco=[]; homsync=[];
for i=1:length(data)
  homstats(i)=AnalyzeData(data(i),@CalcSpikeSync,'ROI_pairs',{'E_V',[0 .5],'E_V',[.5 1]},'maxlag_time',maxlag_time);
  homdN=[homdN homstats(i).pairs.dNsumN]; homsync=[homsync homstats(i).pairs.xcsum_pops];
  homco=[homco homstats(i).pairs.coactivity]; homnco=[homnco homstats(i).pairs.ncoactivity];
end
[homdN;homsync;homco;homnco]
data=SelectData(hetdata,'time_limits',[1000 3000]);
hetdN=[]; hetco=[]; hetnco=[]; hetsync=[];
for i=1:length(data)
  hetstats(i)=AnalyzeData(data(i),@CalcSpikeSync,'ROI_pairs',{'E_V',[0 .5],'E_V',[.5 1]},'maxlag_time',maxlag_time);
  hetdN=[hetdN hetstats(i).pairs.dNsumN]; hetsync=[hetsync hetstats(i).pairs.xcsum_pops];
  hetco=[hetco hetstats(i).pairs.coactivity]; hetnco=[hetnco hetstats(i).pairs.ncoactivity];
end
[hetdN;hetsync;hetco;hetnco]

% gAMPAee=[0 .05];
% gNMDAee=[0 .1];
% gAMPA: 0        .05       0         .05 
% gNMDA: 0        0         .1        .1
% f1=20; f2=42.5; acAMPA1e=4500;
% (betweenblockee=0, withinblockee=1): 
% [homdN;homsync;homco;homnco] = 
%     0.5916    0.5285    0.4095    0.2099    homdN
%     0.3753    0.4389    0.4557    0.3237    homsync
%     0.4725    0.5104    0.3838    0.3545    homco
%     0.0112    0.0123    0.0132    0.0061    homnco
% [hetdN;hetsync;hetco;hetnco] = 
%     0.1419    0.0800    0.0119    0.0497    hetdN
%     0.2593    0.2144    0.2086    0.1517    hetsync
%     0.4955    0.3714    0.4143    0.4040    hetco
%     0.0112    0.0082    0.0076    0.0053    hetnco
% (betweenblockee=1, withinblockee=1): 
% [homdN;homsync;homco;homnco] = 
%     0.5916    0.4079    0.2221    0.0749    homdN
%     0.3753    0.4673    0.5306    0.4496    homsync
%     0.4725    0.5421    0.5091    0.4524    homco
%     0.0112    0.0123    0.0141    0.0098    homnco
% [hetdN;hetsync;hetco;hetnco] = 
%     0.1419    0.0959    0.0086    0.0569
%     0.2593    0.2284    0.2117    0.2165
%     0.4955    0.3885    0.5177    0.4000
%     0.0112    0.0087    0.0073    0.0074

% f1=33; f2=37; acAMPA1e=4500;
% (betweenblockee=1, withinblockee=1): 
% [homdN;homsync;homco;homnco] = 
%     0.0425    0.0075   -0.0281    0.0060
%     0.3958    0.4305    0.4917    0.4186
%     0.5306    0.4286    0.5856    0.4891
%     0.0126    0.0137    0.0132    0.0094
% [hetdN;hetsync;hetco;hetnco] = 
%    -0.0022   -0.0223   -0.0142         0
%     0.2856    0.2417    0.2833    0.3128
%     0.4912    0.5079    0.5120    0.5338
%     0.0120    0.0097    0.0108    0.0095

% f1=20; f2=24; acAMPA1e=4500;
% (betweenblockee=1, withinblockee=1): 
% [homdN;homsync;homco;homnco] = 
%     0.1839    0.1033   -0.1737   -0.0898
%     0.3758    0.4379    0.4088    0.4489
%     0.3786    0.4074    0.3398    0.2683
%     0.0084    0.0092    0.0105    0.0080
% [hetdN;hetsync;hetco;hetnco] = 
%    -0.0669   -0.1194   -0.1346   -0.0539
%     0.2620    0.2962    0.2935    0.4513
%     0.4500    0.3866    0.3761    0.4628
%     0.0082    0.0070    0.0073    0.0093

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1=20; f2=24; acAMPA1e=4500; enoise=acAMPA1e; 
tauGABA=5; % 5, 13
gAMPAee=[0 .05 .1 .15 .2]; % E->E (<=.1)
gNMDAee=[0 .1 .2 .3 .4]; % (later: 0,.1,.2)
gAMPAei=1;
vary={'E','fAMPA1',f1;'E','fAMPA2',f2;'E->E','gNMDA',gNMDAee;'E->E','gAMPA',gAMPAee;'E','baseline',enoise;'E','(acAMPA1,acAMPA2)',acAMPA1e;'E','(dcAMPA1,dcAMPA2)',dcAMPA1e;'(I->E,I->I)','tauGABA',tauGABA;'E->I','gAMPA',gAMPAei};

% simulator options
tspan=[0 3000]; dt=.01; solver='rk1'; compile_flag=1; downsample_factor=10;
save_data_flag=1; cluster_flag=1; memory_limit='32G'; sims_per_job=1;
analysis_functions=@CalcSpikeSync;
analysis_options={'ROI_pairs',{'E_V',[0 .5],'E_V',[.5 1]}};
plot_functions={@PlotData,@PlotData};
plot_options={{'plot_type','rastergram'},{'plot_type','power','xlim',[0 100]}};
simulator_options={'save_data_flag',save_data_flag,'tspan',tspan,'solver',solver,'dt',dt,'compile_flag',compile_flag,'verbose_flag',1,'cluster_flag',cluster_flag,'plot_functions',plot_functions,'plot_options',plot_options,'memory_limit',memory_limit,'downsample_factor',downsample_factor,'analysis_functions',analysis_functions,'analysis_options',analysis_options,'sims_per_job',sims_per_job};

betweenblockee=.5; % 1, .5 kernal scaling for connection b/w E-assemblies (i.e., E1->E2 connected with gsyn=gee*bwblockee)
withinblockee=1;
% connect E-cells within and between blocks
K11=zeros(Ne,Ne); % E->E, [N_pre x N_post]
block=ones(assembly_size)-eye(assembly_size);
for i=1:(Ne/assembly_size) % loop over assemblies
  ind=(i-1)*assembly_size+(1:assembly_size);
  K11(ind,ind)=block;
end
K12=(1-K11-eye(Ne)); % connections b/w blocks 1 and 2
Kee=withinblockee*K11+betweenblockee*K12;
Kee=Kee./repmat(max(1,sum(Kee,1)),[size(Kee,1) 1]);
figure; imagesc(Kee)
spec=ApplyModifications(base,{'E->E','netcon',Kee});
  
rep=1;
study_dir_root=sprintf('studies/sweeps_het-vs-hom_ac%g_EEbetween%g_EEwithin%g_gEEnmda%g-%g_gEEampa%g-%g_gEI%g',acAMPA1e,withinblockee,betweenblockee,gNMDAee(1),gNMDAee(end),gAMPAee(1),gAMPAee(end),gAMPAei);
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
study_dir_root=sprintf('studies/sweeps_het-vs-hom_ac%g_EEbetween%g_EEwithin%g_gEEnmda%g-%g_gEEampa%g-%g_gEI%g',acAMPA1e,withinblockee,betweenblockee,gNMDAee(1),gNMDAee(end),gAMPAee(1),gAMPAee(end),gAMPAei);
het_study_dir=sprintf('%sNe%g_f_%g-%gHz-vs-%g-%gHz_tauI%gms_HET%g_%g',study_dir_root,Ne,f1(1),f1(end),f2(1),f2(end),tauGABA,hetdegree,rep);
[data,het_studyinfo]=SimulateModel(spec,'random_seed',seeds(rep),'study_dir',het_study_dir,'vary',vary,simulator_options{:});


homdata=ImportData(hom_study_dir);
hetdata=ImportData(het_study_dir);

% analysis
maxlag_time=10;
data=SelectData(homdata,'time_limits',[1000 3000]);
homdN=[]; homco=[]; homnco=[]; homsync=[];
for i=1:length(data)
  homstats(i)=AnalyzeData(data(i),@CalcSpikeSync,'ROI_pairs',{'E_V',[0 .5],'E_V',[.5 1]},'maxlag_time',maxlag_time);
  homdN=[homdN homstats(i).pairs.dNsumN]; homsync=[homsync homstats(i).pairs.xcsum_pops];
  homco=[homco homstats(i).pairs.coactivity]; homnco=[homnco homstats(i).pairs.ncoactivity];
end
[homdN;homsync;homco;homnco]
data=SelectData(hetdata,'time_limits',[1000 3000]);
hetdN=[]; hetco=[]; hetnco=[]; hetsync=[];
for i=1:length(data)
  hetstats(i)=AnalyzeData(data(i),@CalcSpikeSync,'ROI_pairs',{'E_V',[0 .5],'E_V',[.5 1]},'maxlag_time',maxlag_time);
  hetdN=[hetdN hetstats(i).pairs.dNsumN]; hetsync=[hetsync hetstats(i).pairs.xcsum_pops];
  hetco=[hetco hetstats(i).pairs.coactivity]; hetnco=[hetnco hetstats(i).pairs.ncoactivity];
end
[hetdN;hetsync;hetco;hetnco]

%{
(gb=1, gw=1):
HOM:
    0.1839    0.1033    0.0871    0.0129         0
    0.3758    0.4379    0.6090    1.3043    2.2545
    0.3786    0.4074    0.3839    0.3469    0.4143
    0.0084    0.0092    0.0087    0.0139    0.0271

   -0.1737   -0.0898   -0.0256   -0.0222    0.0033
    0.4088    0.4489    0.3799    0.9840    2.7463
    0.3398    0.2683    0.3817    0.4348    0.4200
    0.0105    0.0080    0.0053    0.0107    0.0460

   -0.0290   -0.0110   -0.0139    0.0159    0.0131
    0.3990    0.3040    0.4148    1.6223    2.2410
    0.3932    0.3381    0.3971    0.3534    0.3649
    0.0077    0.0051    0.0045    0.0170    0.0316

   -0.0222   -0.0141    0.0101    0.0120    0.0219
    0.3134    0.3189    0.3450    0.4520    1.4103
    0.3355    0.4494    0.4145    0.4733    0.4048
    0.0050    0.0051    0.0042    0.0049    0.0140

   -0.0093    0.0031    0.0192         0         0
    0.2453    0.3212    0.3191    3.0910    3.2645
    0.3765    0.3812    0.5115    0.6500    0.5000
    0.0040    0.0047    0.0044    0.0070    0.0269

HET:
   -0.0669   -0.1194   -0.0693   -0.0298   -0.0058
    0.2620    0.2962    0.4701    0.5857    0.9545
    0.4500    0.3866    0.3806    0.4496    0.4773
    0.0082    0.0070    0.0085    0.0084    0.0098

   -0.1346   -0.0539   -0.0131         0    0.0166
    0.2935    0.4513    0.3478    0.4009    1.3130
    0.3761    0.4628    0.3931    0.4551    0.5678
    0.0073    0.0093    0.0056    0.0059    0.0129

   -0.0409   -0.0010    0.0130    0.0189    0.0103
    0.3623    0.3233    0.3265    0.4957    0.9335
    0.3504    0.3818    0.4331    0.4172    0.5455
    0.0075    0.0061    0.0050    0.0057    0.0106

    0.0146    0.0118    0.0076    0.0126   -0.0006
    0.3153    0.3339    0.2981    0.5054    0.7692
    0.3608    0.4620    0.5145    0.4407    0.4277
    0.0053    0.0055    0.0044    0.0053    0.0080

    0.0025    0.0183    0.0145   -0.0159   -0.0241
    0.3201    0.2358    0.2474    0.4180    2.3345
    0.5325    0.4310    0.5238    0.5163    0.4571
    0.0056    0.0042    0.0041    0.0045    0.0256

%}


%{
(gb=.5, gw=1):
HOM:
    0.1839    0.0295   -0.0052    0.0395   -0.0154
    0.3758    0.3650    0.4225    0.0686    0.9631
    0.3786    0.3458    0.3017    0.1849    0.1667
    0.0084    0.0080    0.0080    0.0009    0.0052

   -0.1185   -0.1159   -0.0336   -0.0230   -0.0248
    0.4882    0.4232    0.1054    0.0603    1.0556
    0.3267    0.3455    0.1923    0.0857    0.1882
    0.0099    0.0091    0.0017    0.0002    0.0033

   -0.0766   -0.0530   -0.0325    0.0657    0.0145
    0.3113    0.2573    0.1223    0.1224    1.1402
    0.2984    0.3382    0.1429    0.1020    0.1848
    0.0058    0.0039    0.0013    0.0007    0.0059

   -0.0302   -0.0018   -0.0031    0.0225    0.0397
    0.2644    0.3356    0.1874    0.2095    1.2152
    0.3944    0.4416    0.1628    0.1739    0.2444
    0.0040    0.0052    0.0014    0.0018    0.0100

   -0.0197    0.0146    0.0259   -0.0231         0
    0.2893    0.3091    0.2876    0.2801    3.0512
    0.4107    0.4364    0.3862    0.0845    0.2821
    0.0049    0.0045    0.0036    0.0001    0.0038

HET:
   -0.0669   -0.0790   -0.0670   -0.0178    0.0124
    0.2620    0.3200    0.3593    0.3152    0.7786
    0.4500    0.3250    0.3840    0.3985    0.3985
    0.0082    0.0081    0.0071    0.0053    0.0076

   -0.1323   -0.0679    0.0252   -0.0063    0.0214
    0.3171    0.3949    0.3715    0.3211    0.6681
    0.4310    0.4375    0.3793    0.3907    0.3551
    0.0086    0.0088    0.0060    0.0051    0.0072

   -0.0377         0    0.0260    0.0368   -0.0108
    0.3053    0.3555    0.3374    0.3778    1.1255
    0.3237    0.4082    0.4500    0.4114    0.4348
    0.0061    0.0068    0.0050    0.0058    0.0137

   -0.0030   -0.0025    0.0060    0.0100   -0.0111
    0.2777    0.3397    0.3481    0.3942    0.2359
    0.3721    0.4878    0.5549    0.5032    0.3143
    0.0052    0.0057    0.0052    0.0058    0.0028

   -0.0042    0.0126    0.0080   -0.0092   -0.0218
    0.2912    0.2510    0.2504    0.3112    1.3369
    0.4971    0.5167    0.5598    0.5114    0.4485
    0.0050    0.0045    0.0044    0.0042    0.0121
%}

% -------------------------------------------------------------------------

f1=20; f2=24; acAMPA1e=4500; enoise=acAMPA1e; 
tauGABA=5; % 5, 13
gAMPAee=[0 .05 .1 .15 .2]; % E->E (<=.1)
gNMDAee=[0 .1 .2 .3 .4]; % (later: 0,.1,.2)
gAMPAei=1;
vary={'E','fAMPA1',f1;'E','fAMPA2',f2;'E->E','gNMDA',gNMDAee;'E->E','gAMPA',gAMPAee;'E','baseline',enoise;'E','(acAMPA1,acAMPA2)',acAMPA1e;'E','(dcAMPA1,dcAMPA2)',dcAMPA1e;'(I->E,I->I)','tauGABA',tauGABA;'E->I','gAMPA',gAMPAei};

% simulator options
tspan=[0 3000]; dt=.01; solver='rk1'; compile_flag=1; downsample_factor=10;
save_data_flag=1; cluster_flag=1; memory_limit='32G'; sims_per_job=1;
analysis_functions=@CalcSpikeSync;
analysis_options={'ROI_pairs',{'E_V',[0 .5],'E_V',[.5 1]}};
plot_functions={@PlotData,@PlotData};
plot_options={{'plot_type','rastergram'},{'plot_type','power','xlim',[0 100]}};
simulator_options={'save_data_flag',save_data_flag,'tspan',tspan,'solver',solver,'dt',dt,'compile_flag',compile_flag,'verbose_flag',1,'cluster_flag',cluster_flag,'plot_functions',plot_functions,'plot_options',plot_options,'memory_limit',memory_limit,'downsample_factor',downsample_factor,'analysis_functions',analysis_functions,'analysis_options',analysis_options,'sims_per_job',sims_per_job};

addpath /projectnb/crc-nak/sherfey/projects/ACC_simulations/ACC_models
spec=base;
spec.connections(1).mechanism_list={'iAMPAee','iNMDAee'};
spec.connections(1).parameters={'tauAMPA',tauAMPA,'tauNMDA',tauNMDA};
m=GenerateModel(spec);
d=SimulateModel(m,'compile_flag',1,'verbose_flag',1,'solver','rk1');
figure; imagesc(d.model.fixed_variables.E_E_iAMPAee_netcon)

vary={'E->E','gw',[0 1];'E->E','gb',[0 1]};
d=SimulateModel(spec,'vary',vary,'compile_flag',1,'verbose_flag',1,'solver','rk1');
figure; 
subplot(2,2,1); imagesc(d(1).model.fixed_variables.E_E_iAMPAee_netcon)
subplot(2,2,2); imagesc(d(2).model.fixed_variables.E_E_iAMPAee_netcon)
subplot(2,2,3); imagesc(d(3).model.fixed_variables.E_E_iAMPAee_netcon)
subplot(2,2,4); imagesc(d(4).model.fixed_variables.E_E_iAMPAee_netcon)

% -------------------------------------------------------------------------

addpath /projectnb/crc-nak/sherfey/projects/ACC_simulations/ACC_models % add path to getEE()

f1=20; f2=24; acAMPA1e=4500; enoise=acAMPA1e; 
tauGABA=5; % 5, 13
gAMPAee=[.05]; % E->E (<=.1)
gNMDAee=[.2]; % (later: 0,.1,.2)
gAMPAei=1;
betweenblockee=[0 .25 .5 .75 1]; % 1, .5 kernal scaling for connection b/w E-assemblies (i.e., E1->E2 connected with gsyn=gee*bwblockee)
withinblockee=1;
vary={'E->E','gw',withinblockee;'E->E','gb',betweenblockee;'E','fAMPA1',f1;'E','fAMPA2',f2;'E->E','gNMDA',gNMDAee;'E->E','gAMPA',gAMPAee;'E','baseline',enoise;'E','(acAMPA1,acAMPA2)',acAMPA1e;'E','(dcAMPA1,dcAMPA2)',dcAMPA1e;'(I->E,I->I)','tauGABA',tauGABA;'E->I','gAMPA',gAMPAei};

rep=2; % 1:10
% for rep=3:10
  
  % simulator options
  tspan=[0 3000]; dt=.01; solver='rk1'; compile_flag=1; downsample_factor=10;
  save_data_flag=1; cluster_flag=1; memory_limit='32G'; sims_per_job=1;
  analysis_functions=@CalcSpikeSync;
  analysis_options={'ROI_pairs',{'E_V',[0 .5],'E_V',[.5 1]}};
  plot_functions={@PlotData,@PlotData};
  plot_options={{'plot_type','rastergram'},{'plot_type','power','xlim',[0 100]}};
  simulator_options={'save_data_flag',save_data_flag,'tspan',tspan,'solver',solver,'dt',dt,'compile_flag',compile_flag,'verbose_flag',1,'cluster_flag',cluster_flag,'plot_functions',plot_functions,'plot_options',plot_options,'memory_limit',memory_limit,'downsample_factor',downsample_factor,'analysis_functions',analysis_functions,'analysis_options',analysis_options,'sims_per_job',sims_per_job};

  spec=base;
  spec.connections(1).mechanism_list={'iAMPAee','iNMDAee'};
  spec.connections(1).parameters={'tauAMPA',tauAMPA,'tauNMDA',tauNMDA};

  study_dir_root=sprintf('studies/sweeps_het-vs-hom_ac%g_EEwithin%g-%g_EEbetween%g-%g_gEEnmda%g-%g_gEEampa%g-%g_gEI%g',acAMPA1e,withinblockee(1),withinblockee(end),betweenblockee(1),betweenblockee(end),gNMDAee(1),gNMDAee(end),gAMPAee(1),gAMPAee(end),gAMPAei);
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
  study_dir_root=sprintf('studies/sweeps_het-vs-hom_ac%g_EEwithin%g-%g_EEbetween%g-%g_gEEnmda%g-%g_gEEampa%g-%g_gEI%g',acAMPA1e,withinblockee(1),withinblockee(end),betweenblockee(1),betweenblockee(end),gNMDAee(1),gNMDAee(end),gAMPAee(1),gAMPAee(end),gAMPAei);
  het_study_dir=sprintf('%sNe%g_f_%g-%gHz-vs-%g-%gHz_tauI%gms_HET%g_%g',study_dir_root,Ne,f1(1),f1(end),f2(1),f2(end),tauGABA,hetdegree,rep);
  [data,het_studyinfo]=SimulateModel(spec,'random_seed',seeds(rep),'study_dir',het_study_dir,'vary',vary,simulator_options{:});

% end

homdata=ImportData(hom_study_dir);
hetdata=ImportData(het_study_dir);

% analysis
maxlag_time=10;
data=SelectData(homdata,'time_limits',[1000 3000]);
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

% ---------------------------------------------------------
% competition test (f1=20Hz, f2=42.5Hz) --> [RESULT]

% model modifications
spec=base;
spec.connections(1).mechanism_list={'iAMPAee','iNMDAee'};
spec.connections(1).parameters={'tauAMPA',tauAMPA,'tauNMDA',tauNMDA};

% parameter modifications
f1=20; f2=42.5; acAMPA1e=4500; enoise=acAMPA1e; 
tauGABA=5; gAMPAee=[.05]; gNMDAee=[.2]; gAMPAei=1;
betweenblockee=0; withinblockee=0; % only varying ratio makes a diffence
rep=2;
vary={'E->E','gw',withinblockee;'E->E','gb',betweenblockee;'E','fAMPA1',f1;'E','fAMPA2',f2;'E->E','gNMDA',gNMDAee;'E->E','gAMPA',gAMPAee;'E','baseline',enoise;'E','(acAMPA1,acAMPA2)',acAMPA1e;'E','(dcAMPA1,dcAMPA2)',dcAMPA1e;'(I->E,I->I)','tauGABA',tauGABA;'E->I','gAMPA',gAMPAei};

% simulation controls
tspan=[0 6000]; dt=.01; solver='rk1'; compile_flag=1; downsample_factor=10;
save_data_flag=0; cluster_flag=0; memory_limit='32G'; sims_per_job=1;
analysis_functions=[]; analysis_options=[];
plot_functions=[]; plot_options=[];
simulator_options={'save_data_flag',save_data_flag,'tspan',tspan,'solver',solver,'dt',dt,'compile_flag',compile_flag,'verbose_flag',1,'cluster_flag',cluster_flag,'plot_functions',plot_functions,'plot_options',plot_options,'memory_limit',memory_limit,'downsample_factor',downsample_factor,'analysis_functions',analysis_functions,'analysis_options',analysis_options,'sims_per_job',sims_per_job};

data=SimulateModel(spec,'random_seed',seeds(rep),'vary',vary,simulator_options{:});
PlotData(data,'plot_type','rastergram');
PlotData(data,'variable','E_nap_h');

dat=SelectData(hetdata,'time_limits',[5000 8000]);
stats=AnalyzeData(dat,@CalcSpikeSync,'ROI_pairs',{'E_V',[0 .5],'E_V',[.5 1]},'maxlag_time',maxlag_time);
stats.pairs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------

rep=1:5;
hetdegree=[.05 .1 .25 1];

% Batch 1: how sync benefit depends on (f1,f2) (f1-f2 and proximity to fr)
% - relate to dependence of competition on (f1,f2)
f1=10:2.5:60; 
f2=f1; 
betweenblockee=[0];
withinblockee=[0 1];

% Batch 2: how sync benefit depends on connectivity within & b/w assemblies
% - interpret wrt learning and assembly integration
f1=[20 30 40 50]; 
f2=[20 30 40 50]; 
betweenblockee=[0 .25 .5 1]; % 1, .5 kernal scaling for connection b/w E-assemblies (i.e., E1->E2 connected with gsyn=gee*bwblockee)
withinblockee=[0 .25 .5 1];

acAMPA1e=4500; 
enoise=acAMPA1e; 
tauGABA=5;
gAMPAee=.05;
gNMDAee=.2;
gAMPAei=1;
vary={'E->E','gw',withinblockee;'E->E','gb',betweenblockee;'E','fAMPA1',f1;'E','fAMPA2',f2;'E->E','gNMDA',gNMDAee;'E->E','gAMPA',gAMPAee;'E','baseline',enoise;'E','(acAMPA1,acAMPA2)',acAMPA1e;'E','(dcAMPA1,dcAMPA2)',dcAMPA1e;'(I->E,I->I)','tauGABA',tauGABA;'E->I','gAMPA',gAMPAei};



% mods={'E->E','netcon','Kee./repmat(max(1,sum(Kee,1)),[size(Kee,1) 1])';
%       'E->E','Kee','withinblockee*K11+betweenblockee*K12';
%       'E->E','K11',K12;'E->E','K12',K12;
%       'E->E','withinblockee',1;'E->E','betweenblockee',1};
% spec=ApplyModifications(base,mods);
% mods={'E->E','netcon','getEE(80,gw,1)';'E->E','gw',1};
% spec=ApplyModifications(base,mods);


% betweenblockee=1; % 1, .5 kernal scaling for connection b/w E-assemblies (i.e., E1->E2 connected with gsyn=gee*bwblockee)
% withinblockee=1;
% % connect E-cells within and between blocks
% K11=zeros(Ne,Ne); % E->E, [N_pre x N_post]
% block=ones(assembly_size)-eye(assembly_size);
% for i=1:(Ne/assembly_size) % loop over assemblies
%   ind=(i-1)*assembly_size+(1:assembly_size);
%   K11(ind,ind)=block;
% end
% K12=(1-K11-eye(Ne)); % connections b/w blocks 1 and 2
% Kee=withinblockee*K11+betweenblockee*K12;
% Kee=Kee./repmat(max(1,sum(Kee,1)),[size(Kee,1) 1]);
% figure; imagesc(Kee)
