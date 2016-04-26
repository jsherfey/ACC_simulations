Hom net (mean param values) with poisson inputs
------------------------------------------------
NEW SCRIPT


Compare cells/nets w/ original vs new mech implementations
homogeneous poisson input: Determine (gsyn,DC) s.t. LFP = (gamma|tauI=5ms) and (beta2|tauI=13ms) & FR=low
rhythmic poisson input (f1,AC=DC): Determine resonant frequency
competition sims (f1,f2): Choose (f1=fr, f2<fr) s.t. |relative activity|>>0

no additive gaussian noise (only noise from poisson)

% define random seeds for all simulations
nsets=5; % one for each parameter, one for hom net, one het varying all params
nreps=5; % repetitions of the same network
nsims=nsets*nreps;
seeds=zeros(nreps,nsets);
for i=1:nsims
  rng('shuffle');
  seeds(i)=getfield(rng,'Seed');
  pause(.01);
end
% save seeds to file (for reproducibility)
seedfile='rand_seeds.par';
fid=fopen(seedfile,'wt');
fprintf(fid,'%10i\n',seeds);
fclose(fid);


% before drawing het params for network:
rng(seeds(i));
% draw parameters for hetergeneous network (uniform random (min,max))
% ...
% run simulation
SimulateModel(model,'random_seed',seeds(i));


% multivariate normal parameters with experimentally-constrained covariance
load('rat-ACd-results/rat-ACd-model_cell-intrinsic-properties.mat','simdata');
collect={'gAHP','gcan','gcat','gh','gkca','gks','gkdr','gnaf','gpas'};
params=[];
for i=1:numel(collect)
  params=cat(2,params,simdata.(collect{i}));
end
MU = mean(params,1);
SIGMA = cov(params);
values = mvnrnd(MU,SIGMA);  % parameter values for one cell in one simulation
values = max(0,values);     % turn off currents w/ negative max conductance
% values=max(0,mvnrnd(MU,SIGMA,ncells));




%{
Model param ranges
------------------------------
From ACd_genbatches.m: (~/models/dnsim/aro/ACd_hetero-vs-homo)

% mean max conductance for homogeneous networks
gAHP_mu = mean([.1 0]);
taurca_mu = mean([28.5714,57.1429]);
gcan_mu = mean([0,.01]);
gcat_mu = mean([0,.01]);
gh_mu = mean([0,2]);
eh_mu = mean([-30,-10]);
gks_mu = mean([0,1]); % .1
gkdr_mu = mean([3,9]);
gkca_mu = mean([0,2]); % .1
gnaf_mu = 75;%mean([25,100]); % 75
gnap_mu = 0.0005;
epas_mu = mean([-90,-60]);
gpas_mu = mean([.04,.1]); % .06
noise_mu = mean([0,2]);
% extra mechanism parameters
gNa_DS = 80;
% gK_DS = 36;

default_params = [gAHP_mu,taurca_mu,gcan_mu,gcat_mu,gh_mu,eh_mu,gks_mu,...
  gkdr_mu,gkca_mu,gnaf_mu,gnap_mu,epas_mu,gpas_mu,noise_mu,0];%,gK_DS];

variable_params = {[0 gAHP_mu],[0 taurca_mu],[0 gcan_mu],[0 gcat_mu],...
  [0 gh_mu],[0 eh_mu],[0 gks_mu],[4 gkdr_mu],[0 gkca_mu],[0 gnaf_mu],...
  [0 gnap_mu],[-60 epas_mu],[0 gpas_mu],[0 noise_mu],[0 gNa_DS]};%,[0 gK_DS]

param_labels = {'gAHP','taurca','gcan','gcat','gh','eh','gks',...
                'gkdr','gkca','gnaf','gnap','epas','gpas','noise',...
                'gNa'};%,'gK'};

params_to_vary={'gAHP','gcan','gcat','gh','gks','gkca','gkdr','epas'};%,'gnaf','gNa'

hetero_mechanisms = {...
'AHPdist','cadyndist','candist','catdist','hdist','iksdist','kdrdist',...
'kcadist','nafdist','nap','pasdist','randndist','InputGenerator',...
'iNa_DS'};%,'iK_DS'};

homo_mechanisms = hetero_mechanisms; % <-- same model for homo vs hetero (only parameter values differ)

num_params = length(default_params);              
[var_i,jnk] = match_str(param_labels,params_to_vary);
const_i = setdiff(1:num_params,var_i);

% finalize list of parameters that are constant and that are
% (1) variable across nets in homo sims as well as
% (2) variable across cells in hetero sims
variable_params = variable_params(var_i);
variable_params_labels = param_labels(var_i)
num_variable_params = length(variable_params);
constant_params = default_params(const_i);
constant_params_labels = param_labels(const_i)
num_constant_params = length(constant_params);

% homo: all mechs effectively w/o "dist" and params=[mean or 0 per param]
% hetero: all mechs w/ "dist" and params=[grid per params]

% generate parameter space spanned by population
v=variable_params;
n=numel(v);
x=cell(n,1);
[x{1:n,1}]=ndgrid(v{end:-1:1});
p=reshape(cat(n+1,x{:}),[],n);
paramspace=p(:,end:-1:1);
% disp(p); size(p)

% population sizes
nE = size(paramspace,1) % number of E-cells
nI = ceil(.25*nE)

% prepare static parameter key-value pairs
constant_keyval={};
for i=1:num_constant_params
  constant_keyval{end+1}=constant_params_labels{i};
  constant_keyval{end+1}=repmat(constant_params(i),[nE 1]);
end

hetero_keyval={};
for i=1:num_params
  hetero_keyval{end+1}=param_labels{i};
  if ismember(param_labels{i},variable_params_labels)
    % vary this parameter across the heterogeneous population
    index=find(ismember(variable_params_labels,param_labels{i}));
    hetero_keyval{end+1}=paramspace(:,index);
  else
    % use the default value for all cells in the population
    index=find(ismember(constant_params_labels,param_labels{i}));
    hetero_keyval{end+1}=repmat(constant_params(index),[nE 1]);
  end
end

cd /home/jason/models/dnsim/aro/ACd_hetero-vs-homo/mechanisms

scc_parentdir='/projectnb/crc-nak/sherfey/projects/heterogeneity';
PROJECT_DIR = fullfile(scc_parentdir,'ACd_hetero-vs-homo');

local_parentdir='/home/jason/models/dnsim/aro';
PROJECT_DIR = fullfile(local_parentdir,'ACd_hetero-vs-homo');

mechpath = fullfile(PROJECT_DIR,'mechanisms');
studydriver=fullfile(PROJECT_DIR,'ACd_simstudy_driver.m'); % 'ACd_net_simstudy.m');
cd(PROJECT_DIR);

%}