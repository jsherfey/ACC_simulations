% Clustering IPs (simulated data)
% -----------------------------------

addpath ardid-cluster-method
load('rat-ACd-results/rat-ACd-model_cell-intrinsic-properties.mat','simdata','expdata');

% select IPs to cluster
sel=1:length(simdata.RMP);%1:100;
IP=[simdata.AHPtime2trough(sel) simdata.spikewidth(sel) simdata.threshrate(sel) simdata.RMP(sel) simdata.steprate(sel)];
% sel=1:84;
% IP=[expdata.AHPtime2trough(sel) expdata.spikewidth(sel) expdata.threshrate(sel) expdata.RMP(sel) expdata.steprate(sel)];
IPlabels=expdata.properties;

% ---- CLUSTERING ---- %
[ncell,nmet]=size(IP);
resultsdir='sim_clusters';
mkdir(resultsdir);

% replace NaN's by zeros
IP(isnan(IP))=0;

% scale each IP to [0,1]
for i=1:size(IP,2)
  IP(:,i)=IP(:,i)-min(IP(:,i));
  IP(:,i)=IP(:,i)/max(IP(:,i));
end

type = 'Kmeans';
datanorm4Cluster=IP;

numClusters=6;%2:6;
visible='on';
tic
[labels,clustFilt,percElements] = function_pairingOfClusterElements(type,datanorm4Cluster,numClusters,visible,resultsdir);
toc

disp('num clusters and percent cells clustered');
numClusters,percElements

for nc=numClusters%4:6 % number of clusters

celltype=repmat({'broad'},[1 ncell]);
visible='on';
ClusterIndex = find(numClusters==nc);
main_dendrogram_clustersCentroids;

varIP = var(IP);
[sortedVar_datanormAll,indSortedVar_datanormAll] = sort(varIP,2,'descend');																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																				
explainedVar_datanormAll = 100*sortedVar_datanormAll/sum(sortedVar_datanormAll);
cumExplainedVar_datanormAll = cumsum(explainedVar_datanormAll);

% cutoff 80% --> exclude IP #s ...
% resort IPs and plot

datanorm4Cluster = IP(:,indSortedVar_datanormAll);
IPlabels = IPlabels(indSortedVar_datanormAll);
ClusterIndex = find(numClusters==nc);
main_clustering_heatMaps;

figure
imagesc(1:numMeasures,1:numNeurons,sortedData);
colormap(1-gray); axis tight; %axis off; 
try for i=1:length(boundaries), hline(boundaries(i),'color','k','linewidth',3); end; end
title(sprintf('boundaries=[%s] (%g%% clustered)',num2str(boundaries'),percElements(ClusterIndex)))
set(gca,'xtick',1:length(IPlabels),'xticklabel',IPlabels); colorbar; axis xy
file = [resultsdir,'/custom_',type,num2str(numClust),'ClusteringResults_sortedHeatMap'];
% print(gcf,'-djpeg',file);
% print(gcf,'-depsc',file);

disp('cumulative explained variance and sorted IP list');
IPlabels
cumExplainedVar_datanormAll

% figure; imagesc(sortedData)
binds = [1 boundaries' size(sortedData,1)+1];
avgIPs = zeros(nc,size(sortedData,2));
for i=1:nc
  dat=sortedData(binds(i):binds(i+1)-1,:);
  avgIPs(i,:)=mean(dat,1);
end

disp('class boundaries and avg normalized IPs');
binds,avgIPs

% cellIDs=IPs_id;
% fprintf('cluster_cell_ids={};\n');
% for i=1:nc
%   fprintf('cluster_cell_ids{%g} = [',i); 
%   ids=cellIDs(clustFilt{ClusterIndex,clustOrder(i)});
%   for j=1:length(ids)
%     fprintf('%g ',ids(j));
%   end
%   fprintf('];\n');
% end
% 
% file = [resultsdir,'/',type,'_ClusteringResults_cluster_cell_indices'];
% fid=fopen(file,'a+');
% fprintf(fid,'%g clusters (%g%% cells clustered): IPs (%s) (%g cells)\n',numClust,percElements(ClusterIndex),ipstr,ncell);
% fprintf(fid,'cluster_cell_ids={};\n');
% for i=1:nc
%   fprintf(fid,'cluster_cell_ids{%g} = [',i); 
%   ids=cellIDs(clustFilt{ClusterIndex,clustOrder(i)});
%   for j=1:length(ids)
%     fprintf(fid,'%g ',ids(j));
%   end
%   fprintf(fid,'];\n');
% end
% fprintf(fid,'\n');
% fclose(fid);
% save([resultsdir,'/',type,num2str(numClust),'ClusteringResults.mat']);

end

% save simdata_6clusters.mat