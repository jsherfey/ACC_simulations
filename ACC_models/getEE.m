function Kee=getEE(Ne,gw,gb)
if nargin<3, gb=1; end
if nargin<2, gw=1; end

betweenblockee=gb; % 1, .5 kernal scaling for connection b/w E-assemblies (i.e., E1->E2 connected with gsyn=gee*bwblockee)
withinblockee=gw;
% connect E-cells within and between blocks
K11=zeros(Ne,Ne); % E->E, [N_pre x N_post]
assembly_size=Ne/2;
block=ones(assembly_size)-eye(assembly_size);
for i=1:(Ne/assembly_size) % loop over assemblies
  ind=(i-1)*assembly_size+(1:assembly_size);
  K11(ind,ind)=block;
end
K12=(1-K11-eye(Ne)); % connections b/w blocks 1 and 2
Kee=withinblockee*K11+betweenblockee*K12;
Kee=Kee./repmat(max(1,sum(Kee,1)),[size(Kee,1) 1]);
