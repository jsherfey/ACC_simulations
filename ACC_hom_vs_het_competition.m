% Hom vs Het net with poisson inputs (24-Apr-2016)
% ------------------------------------------------
% Hom/Het Net Simulations:
% 1. Do competition sims for het nets (multivariate normal distribution)
%    note: save enough info for the sims to be reproducible
% 2. Generate rastergrams and box plots comparing hom vs het net responses
% ------------------------------------------------

% examine params from models with experimentally-constrained IPs:
ipfile='/projectnb/crc-nak/sherfey/projects/ACC_simulations/rat-ACd-results/rat-ACd-model_cell-intrinsic-properties.mat';
load(ipfile,'simdata','expdata');
x=1:length(simdata.Cm); p=simdata.parameters;
h1=figure; h2=figure; nbin=100;
for i=1:length(p)
  figure(h1); subplot(3,5,i); plot(x,simdata.(p{i}),'o'); ylabel(p{i}); xlabel('cell index');
  figure(h2); subplot(3,5,i); hist(simdata.(p{i}),nbin); xlabel(p{i}); ylabel('count');
end

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

%% Simulations: example het vs hom: (f1,f2)=(5Hz,20Hz) (n=10)

% limit the number of cores to <= scc limit
maxNumCompThreads(4);

% number realizations of hom or het network
num_realizations=2; enoise=4500; acAMPA1e=4500; hetdegree=.2;

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

% input parameters
f1=35; f2=20; tauGABA=13; % (35,5), (35,20), (35,40)
vary={'E','fAMPA1',f1;'E','fAMPA2',f2;'E','baseline',enoise;'E','(acAMPA1,acAMPA2)',acAMPA1e;'E','(dcAMPA1,dcAMPA2)',dcAMPA1e;'(I->E,I->I)','tauGABA',tauGABA};

% simulator options
tspan=[0 3000]; dt=.01; solver='rk1'; compile_flag=1; downsample_factor=10;
save_data_flag=1; study_dir_root=sprintf('studies/example_het-vs-hom_ac%g_hetdeg%g',acAMPA1e,hetdegree);
cluster_flag=0; memory_limit='64G';
plot_functions={@PlotData,@PlotData};
plot_options={{'plot_type','rastergram'},{'plot_type','power','xlim',[0 100]}};
simulator_options={'save_data_flag',save_data_flag,'tspan',tspan,'solver',solver,'dt',dt,'compile_flag',compile_flag,'verbose_flag',1,'cluster_flag',cluster_flag,'plot_functions',plot_functions,'plot_options',plot_options,'memory_limit',memory_limit,'downsample_factor',downsample_factor};

% homogeneous network simulations
clear homdata
for rep=1:num_realizations
  study_dir=sprintf('%s/Ne%g_f_%gHz-vs-%gHz_tauI%gms_HOM_%g',study_dir_root,Ne,f1,f2,tauGABA,rep);
  [data,studyinfo]=SimulateModel(base,'random_seed',seeds(rep),'study_dir',study_dir,'vary',vary,simulator_options{:});
  if rep>1, data=orderfields(data,homdata); end
  homdata(rep)=data;
  clear data
end

% heterogeneous network simulations
spec=base;
clear hetdata
for rep=1:num_realizations % loop over network realizations
  rng(seeds(rep)); % set the random seed for drawing parameter values
  % draw heterogeneous parameters from multivariate distribution preserving
  % correlations b/w parameters around their median values:
  values = mvnrnd(MU,hetdegree*SIGMA,Ne/2); % parameter values for one realization of E-cell population
  values = cat(1,values,values); % copy for each assembly
  values = max(0,values);       % turn off currents w/ negative max conductance
  % update model specification with heterogeneous parameter values
  mods={};
  for i=1:numel(het)
    mods=cat(1,mods,{'E',het{i},values(:,i)'});
  end
  spec=ApplyModifications(spec,mods);
  study_dir=sprintf('%s/Ne%g_f_%gHz-vs-%gHz_tauI%gms_HET_%g',study_dir_root,Ne,f1,f2,tauGABA,rep);
  [data,studyinfo]=SimulateModel(spec,'random_seed',seeds(rep),'study_dir',study_dir,'vary',vary,simulator_options{:});
  if rep>1, data=orderfields(data,hetdata); end
  hetdata(rep)=data;
  clear data
end

PlotData(homdata,'plot_type','rastergram')
PlotData(hetdata,'plot_type','rastergram')

% analysis
maxlag_time=10;
data=SelectData(homdata,'time_limits',[1000 3000]);
homdN=[]; homco=[]; homnco=[]; homsync=[];
for i=1:length(data)
  homstats(i)=AnalyzeData(data(i),@CalcSpikeSync,'ROI_pairs',{'E_V',[0 .5],'E_V',[.5 1]},'maxlag_time',maxlag_time);
  homdN=[homdN homstats(i).pairs.dNsumN]; homsync=[homsync homstats(i).pairs.xcsum_pops];
  homco=[homco homstats(i).pairs.coactivity]; homnco=[homnco homstats(i).pairs.ncoactivity];
end
data=SelectData(hetdata,'time_limits',[1000 3000]);
hetdN=[]; hetco=[]; hetnco=[]; hetsync=[];
for i=1:length(data)
  hetstats(i)=AnalyzeData(data(i),@CalcSpikeSync,'ROI_pairs',{'E_V',[0 .5],'E_V',[.5 1]},'maxlag_time',maxlag_time);
  hetdN=[hetdN hetstats(i).pairs.dNsumN]; hetsync=[hetsync hetstats(i).pairs.xcsum_pops];
  hetco=[hetco hetstats(i).pairs.coactivity]; hetnco=[hetnco hetstats(i).pairs.ncoactivity];
end
[homdN;homsync;homco;homnco]
[hetdN;hetsync;hetco;hetnco]



%% experiment

% ##############################################
% new coactivity measure: eliminate outliers
% justification: rare coactive periods occur on cycles where the precessing
% inputs transiently synchronize. this transient spike in coactivity is not
% present in the het networks.
i=1;
stats=homstats(i);
edges=stats(i).pairs.bin_edges;
n1=stats(i).pairs.binned_spikes1; homn1=n1;
n2=stats(i).pairs.binned_spikes2; homn2=n2;
stats=hetstats(i);
n1=stats(i).pairs.binned_spikes1; hetn1=n1;
n2=stats(i).pairs.binned_spikes2; hetn2=n2;
th=99; % percentile threshold
tmp=homn1.*homn2; rm=tmp>prctile(tmp,th); tmp(rm)=0; ncoactive=sum(tmp)/(sum(homn1(~rm))*sum(homn2(~rm)))
tmp=hetn1.*hetn2; rm=tmp>prctile(tmp,th); tmp(rm)=0; ncoactive=sum(tmp)/(sum(hetn1(~rm))*sum(hetn2(~rm)))
tmp=homn1.*homn2; rm=tmp>(mean(tmp(:))+5*std(tmp(:))); tmp(rm)=0; ncoactive=sum(tmp)/(sum(homn1(~rm))*sum(homn2(~rm)))
tmp=hetn1.*hetn2; rm=tmp>(mean(tmp(:))+5*std(tmp(:))); tmp(rm)=0; ncoactive=sum(tmp)/(sum(hetn1(~rm))*sum(hetn2(~rm)))

data=SelectData(homdata,'time_limits',[1000 3000]);
s1=AnalyzeData(data(1),@CalcSpikeSync,'ROI_pairs',{'E_V',[0 .5],'E_V',[.5 1]},'maxlag_time',maxlag_time);
data=SelectData(hetdata,'time_limits',[1000 3000]);
s2=AnalyzeData(data(1),@CalcSpikeSync,'ROI_pairs',{'E_V',[0 .5],'E_V',[.5 1]},'maxlag_time',maxlag_time);
s1.pairs.ncoactivity
s2.pairs.ncoactivity

% ##############################################

% input parameters
f1=5:7.5:80; f2=f1; acAMPA1e=4500; enoise=acAMPA1e; 
tauGABA=5; % 5, 13
gNMDAee=0; % (later: 0,.1,.2)
gAMPAei=1;
vary={'E','fAMPA1',f1;'E','fAMPA2',f2;'E','baseline',enoise;'E','(acAMPA1,acAMPA2)',acAMPA1e;'E','(dcAMPA1,dcAMPA2)',dcAMPA1e;'(I->E,I->I)','tauGABA',tauGABA;'E->E','gNMDA',gNMDAee;'E->I','gAMPA',gAMPAei};

% simulator options
tspan=[0 3000]; dt=.01; solver='rk1'; compile_flag=1; downsample_factor=10;
save_data_flag=1; cluster_flag=1; memory_limit='32G'; sims_per_job=10;
analysis_functions=@CalcSpikeSync;
analysis_options={'ROI_pairs',{'E_V',[0 .5],'E_V',[.5 1]}};
plot_functions={@PlotData,@PlotData};
plot_options={{'plot_type','rastergram'},{'plot_type','power','xlim',[0 100]}};
simulator_options={'save_data_flag',save_data_flag,'tspan',tspan,'solver',solver,'dt',dt,'compile_flag',compile_flag,'verbose_flag',1,'cluster_flag',cluster_flag,'plot_functions',plot_functions,'plot_options',plot_options,'memory_limit',memory_limit,'downsample_factor',downsample_factor,'analysis_functions',analysis_functions,'analysis_options',analysis_options,'sims_per_job',sims_per_job};

rep=1;
study_dir_root=sprintf('studies/sweeps_het-vs-hom_ac%g_gEEnmda%g_gEI%g',acAMPA1e,gNMDAee,gAMPAei);
hom_study_dir=sprintf('%sNe%g_f_%g-%gHz-vs-%g-%gHz_tauI%gms_HOM_%g',study_dir_root,Ne,f1(1),f1(end),f2(1),f2(end),tauGABA,rep);
[data,hom_studyinfo]=SimulateModel(base,'random_seed',seeds(rep),'study_dir',hom_study_dir,'vary',vary,simulator_options{:});

unix(sprintf('cat %s/pbsout/sim_job1.out',hom_studyinfo.simulations(1).batch_dir));

% het sims
hetdegree=.2;
rng(seeds(rep)); % set the random seed for drawing parameter values
values = mvnrnd(MU,hetdegree*SIGMA,Ne/2); % parameter values for one realization of E-cell population
values = cat(1,values,values); % copy for each assembly
values = max(0,values);       % turn off currents w/ negative max conductance
% update model specification with heterogeneous parameter values
mods={};
for i=1:numel(het)
  mods=cat(1,mods,{'E',het{i},values(:,i)'});
end
spec=ApplyModifications(base,mods);
het_study_dir=sprintf('%sNe%g_f_%g-%gHz-vs-%g-%gHz_tauI%gms_HET%g_%g',study_dir_root,Ne,f1(1),f1(end),f2(1),f2(end),tauGABA,hetdegree,rep);
[data,het_studyinfo]=SimulateModel(spec,'random_seed',seeds(rep),'study_dir',het_study_dir,'vary',vary,simulator_options{:});

unix(sprintf('cat %s/pbsout/sim_job1.out',het_studyinfo.simulations(1).batch_dir));

% compare hom vs het
homstats=ImportResults(hom_studyinfo,analysis_functions); 
f1=[homstats.E_fAMPA1]; uf1=unique(f1); nf1=length(uf1);
f2=[homstats.E_fAMPA2]; uf2=unique(f2); nf2=length(uf2);
homDN=zeros(nf1,nf2);
homCO=zeros(nf1,nf2);
for i=1:nf1
  for j=1:nf2
    stats=homstats(f1==uf1(i)&f2==uf2(j));
    homDN(i,j)=stats.pairs.dNsumN;
    homCO(i,j)=stats.pairs.ncoactivity;
  end
end
hetstats=ImportResults(het_studyinfo,analysis_functions); 
hetDN=zeros(nf1,nf2);
hetCO=zeros(nf1,nf2);
for i=1:nf1
  for j=1:nf2
    stats=hetstats(f1==uf1(i)&f2==uf2(j));
    hetDN(i,j)=stats.pairs.dNsumN;
    hetCO(i,j)=stats.pairs.ncoactivity;
  end
end

[homCO(3,5) homDN(3,5)]
[hetCO(3,5) hetDN(3,5)]

ind=find(f1==35 & f2==20);
homdat=load(hom_studyinfo.simulations(ind).data_file);
hetdat=load(het_studyinfo.simulations(ind).data_file);
dat=SelectData(homdat,'time_limits',[1000 3000]);
s1=AnalyzeData(dat,@CalcSpikeSync,'ROI_pairs',{'E_V',[0 .5],'E_V',[.5 1]},'maxlag_time',10);
dat=SelectData(hetdat,'time_limits',[1000 3000]);
s2=AnalyzeData(dat,@CalcSpikeSync,'ROI_pairs',{'E_V',[0 .5],'E_V',[.5 1]},'maxlag_time',10);
s1.pairs.ncoactivity
s2.pairs.ncoactivity

% manually load and re-analyze competition
hetDN=zeros(nf1,nf2);
hetCO=zeros(nf1,nf2);
hetCOn=zeros(nf1,nf2);
tstart=tic;
for i=1:nf1
  for j=1:nf2
    ind=find(f1==uf1(i)&f2==uf2(j));
    homdat=load(hom_studyinfo.simulations(ind).data_file);
    hetdat=load(het_studyinfo.simulations(ind).data_file);
    dat=SelectData(homdat,'time_limits',[1000 3000]);
    s1=AnalyzeData(dat,@CalcSpikeSync,'ROI_pairs',{'E_V',[0 .5],'E_V',[.5 1]},'maxlag_time',10);
    dat=SelectData(hetdat,'time_limits',[1000 3000]);
    s2=AnalyzeData(dat,@CalcSpikeSync,'ROI_pairs',{'E_V',[0 .5],'E_V',[.5 1]},'maxlag_time',10);
    homDN(i,j)=s1.pairs.dNsumN;
    homCO(i,j)=s1.pairs.coactivity;
    homCOn(i,j)=s1.pairs.ncoactivity;
    hetDN(i,j)=s2.pairs.dNsumN;
    hetCO(i,j)=s2.pairs.coactivity;
    hetCOn(i,j)=s2.pairs.ncoactivity;
  end
  fprintf('%g of %g (%gmin)\n',i,nf1,toc(tstart)/60);
end
figure
subplot(3,3,1); imagesc(f1,f2,homDN); title('hom dN'); caxis([-1 1]); axis xy
subplot(3,3,2); imagesc(f1,f2,homCO); title('hom co'); axis xy; colorbar%caxis([0 .04])
subplot(3,3,3); imagesc(f1,f2,homCOn); title('hom co-n'); axis xy; caxis([0 .04])
subplot(3,3,4); imagesc(f1,f2,hetDN); title('het dN'); caxis([-1 1]); axis xy
subplot(3,3,5); imagesc(f1,f2,hetCO); title('het co'); axis xy; colorbar%caxis([0 .04])
subplot(3,3,6); imagesc(f1,f2,hetCOn); title('het co-n'); axis xy; caxis([0 .04])
subplot(3,3,7); imagesc(f1,f2,abs(hetDN)<abs(homDN)); title('|dN| (het<hom)'); axis xy
subplot(3,3,8); imagesc(f1,f2,homCO<hetCO); title('co (het>hom)'); axis xy
subplot(3,3,9); imagesc(f1,f2,homCOn<hetCOn); title('co-n (het>hom)'); axis xy


% todo: add time_limits to dN and sync calcs in CalcSpikeSync

%% Simulations: frequency sweeps: (f1,f2) (n=10)

% number realizations of hom or het network
num_realizations=5; 

% load random seeds for all simulations (use same for hom and het)
num_seeds=1e3; % <= max_num_realizations
seedfile=sprintf('%gseeds.mat',num_seeds);
load(seedfile,'seeds');

% input parameters
f1=40; f2=5:60; tauGABA=13;
vary={'E','fAMPA1',f1;'E','fAMPA2',f2;'E','baseline',enoise;'E','(acAMPA1,acAMPA2)',acAMPA1e;'E','(dcAMPA1,dcAMPA2)',dcAMPA1e;'(I->E,I->I)','tauGABA',tauGABA};

% simulator options
tspan=[0 3000]; dt=.01; solver='rk1'; compile_flag=1; downsample_factor=10;
save_data_flag=1; study_dir_root='studies/het-vs-hom';
cluster_flag=1; memory_limit='32G'; sims_per_job=10;
analysis_functions=@CalcSpikeSync;
analysis_options={'ROI_pairs',{'E_V',[0 .5],'E_V',[.5 1]}};
plot_functions={@PlotData,@PlotData};
plot_options={{'plot_type','rastergram'},{'plot_type','power','xlim',[0 100]}};
simulator_options={'save_data_flag',save_data_flag,'tspan',tspan,'solver',solver,'dt',dt,'compile_flag',compile_flag,'verbose_flag',1,'cluster_flag',cluster_flag,'plot_functions',plot_functions,'plot_options',plot_options,'memory_limit',memory_limit,'downsample_factor',downsample_factor,'analysis_functions',analysis_functions,'analysis_options',analysis_options,'sims_per_job',sims_per_job};

% homogeneous network simulations
clear homdata
for rep=1:num_realizations
  study_dir=sprintf('%s/Ne%g_f_%g-%gHz-vs-%g-%gHz_tauI%gms_HOM_%g',study_dir_root,Ne,f1(1),f1(end),f2(1),f2(end),tauGABA,rep);
  [data,studyinfo]=SimulateModel(base,'random_seed',seeds(rep),'study_dir',study_dir,'vary',vary,simulator_options{:});
  if ~isempty(data)
    if rep>1, data=orderfields(data,homdata); end
    homdata(rep)=data;
  end
  clear data
end

% heterogeneous network simulations
spec=base;
clear hetdata
for rep=1:num_realizations % loop over network realizations
  rng(seeds(rep)); % set the random seed for drawing parameter values
  % draw heterogeneous parameters from multivariate distribution preserving
  % correlations b/w parameters around their median values:
  values = mvnrnd(MU,hetdegree*SIGMA,Ne/2); % parameter values for one realization of E-cell population
  values = cat(1,values,values); % copy for each assembly
  values = max(0,values);       % turn off currents w/ negative max conductance
  % update model specification with heterogeneous parameter values
  mods={};
  for i=1:numel(het)
    mods=cat(1,mods,{'E',het{i},values(:,i)'});
  end
  spec=ApplyModifications(spec,mods);
  study_dir=sprintf('%s/Ne%g_f_%g-%gHz-vs-%g-%gHz_tauI%gms_HET_%g',study_dir_root,Ne,f1(1),f1(end),f2(1),f2(end),tauGABA,rep);
  [data,studyinfo]=SimulateModel(spec,'random_seed',seeds(rep),'study_dir',study_dir,'vary',vary,simulator_options{:});
  if ~isempty(data)
    if rep>1, data=orderfields(data,hetdata); end
    hetdata(rep)=data;
  end
  clear data
end

if cluster_flag
  !qstat -u sherfey
  unix(sprintf('cat %s/pbsout/sim_job1.out',studyinfo.simulations(1).batch_dir));
end

% load analysis results
for rep=1:num_realizations
  study_dir=sprintf('%s/Ne%g_f_%g-%gHz-vs-%g-%gHz_tauI%gms_HOM_%g',study_dir_root,Ne,f1(1),f1(end),f2(1),f2(end),tauGABA,rep);
  homstats{rep}=ImportResults(study_dir,analysis_functions); 
  study_dir=sprintf('%s/Ne%g_f_%g-%gHz-vs-%g-%gHz_tauI%gms_HET_%g',study_dir_root,Ne,f1(1),f1(end),f2(1),f2(end),tauGABA,rep);
  hetstats{rep}=ImportResults(study_dir,analysis_functions); 
end

% figure; 
% r1=nanmean(homstats(1).pairs.r1,2);
% r2=nanmean(homstats(1).pairs.r2,2);
% rms=sqrt(mean((r1-r2).^2));
% rms/mean([r1;r2]) % 3.7669
% corr([r1 r2]) % .3597
% sum(r1.*r2) % .5566
% sum(r1.*r2)/(sum(r1)*sum(r2)) % .0026
% subplot(2,1,1); plot(1:2001,r1,'b',1:2001,r2+max(r1),'r'); title('hom'); axis tight
% r1=nanmean(hetstats(1).pairs.r1,2);
% r2=nanmean(hetstats(1).pairs.r2,2);
% rms=sqrt(mean((r1-r2).^2));
% rms/mean([r1;r2]) % 2.5966
% corr([r1 r2]) % .1597
% sum(r1.*r2) % .1969
% sum(r1.*r2)/(sum(r1)*sum(r2)) % .
% subplot(2,1,2); plot(1:2001,r1,'b',1:2001,r2+max(r1),'r'); title('het'); axis tight


% figure; 
% i=1;
% stats=homstats(i);
% edges=stats(i).pairs.bin_edges;
% n1=stats(i).pairs.binned_spikes1; homn1=n1;
% n2=stats(i).pairs.binned_spikes2; homn2=n2;
% subplot(2,1,1); plot(edges,n1,'bo-',edges,n2,'ro-'); title('hom');
% stats=hetstats(i);
% n1=stats(i).pairs.binned_spikes1; hetn1=n1;
% n2=stats(i).pairs.binned_spikes2; hetn2=n2;
% subplot(2,1,2); plot(edges,n1,'bo-',edges,n2,'ro-'); title('het');
% 
% n1=homn1; n2=homn2;
% ncoactive=sum(n1.*n2)/(sum(n1)*sum(n2))
% ncoactive=sum(homn1.*homn2)/(sum(homn1)*sum(homn2))
% ncoactive=sum(hetn1.*hetn2)/(sum(hetn1)*sum(hetn2))

% ncoactive=sum(sqrt(homn1.*homn2))/(sum(homn1)*sum(homn2))
% ncoactive=sum(sqrt(hetn1.*hetn2))/(sum(hetn1)*sum(hetn2))
% 
% i1=homn1>0;i2=homn2>0; ncoactive=median(homn1(i1&i2).*homn2(i1&i2))/(median(homn1(i1))*median(homn2(i2)))
% i1=hetn1>0;i2=hetn2>0; ncoactive=median(hetn1(i1&i2).*hetn2(i1&i2))/(median(hetn1(i1))*median(hetn2(i2)))
% 
% stats=homstats(1);
% stats=hetstats(1);
% raster1=stats.pairs.raster1;
% raster2=stats.pairs.raster2;
% winsize=5;
% cnt=0;
% for i=1:size(raster1,1)
%   if any((raster2(:,1)>=raster1(i,1)-winsize)&(raster2(:,1)<=raster1(i,1)+winsize))
%     cnt=cnt+1;
%   end
% end
% p12=cnt/size(raster1,1);
% raster1=stats.pairs.raster2;
% raster2=stats.pairs.raster1;
% winsize=5;
% cnt=0;
% for i=1:size(raster1,1)
%   if any((raster2(:,1)>=raster1(i,1)-winsize)&(raster2(:,1)<=raster1(i,1)+winsize))
%     cnt=cnt+1;
%   end
% end
% p21=cnt/size(raster1,1);
% [p12 p21]
