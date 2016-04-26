% savesvg=0; if savesvg==0, visible='on'; else visible='off'; end

% clear;
% close all;
% clc;

% type = 'Kmeans';
% resultsdir = 'clusteringResults';
% matfile = [resultsdir,'/',type,'ClusteringResults'];
% load(matfile);
numClust=numClusters(ClusterIndex);
matfile = [resultsdir,'/',type,num2str(numClust),'ClusteringResults'];
if ~exist(resultsdir,'dir'), mkdir(resultsdir); end

% calculating centers of k-means clusters
% numClusters = 7;
% offset = 4;
% ClusterIndex = numClusters-offset;

numOfDim=size(datanorm4Cluster,2);

nC=length(clustFilt(1,:));
k=0; centers=[]; centersSD=[]; numElem=[];
tmp=1:size(datanorm4Cluster,1);
for i=1:nC
    if(~isempty([clustFilt{ClusterIndex,i}]))
        k = k+1;
        centers(k,:) = mean(datanorm4Cluster([clustFilt{ClusterIndex,i}],:),1);
        centersSD(k,:) = std(datanorm4Cluster([clustFilt{ClusterIndex,i}],:),0,1);
        numElem(k) = length([clustFilt{ClusterIndex,i}]);
        tmp = setdiff(tmp,[clustFilt{ClusterIndex,i}]);
    end
end

% remaining elements if any
if ~isempty(tmp)
  k = k+1;
  clustFilt{ClusterIndex,k} = tmp;
  numElem(k) = size(datanorm4Cluster,1)-sum(numElem);
  centers(k,:) = mean(datanorm4Cluster([clustFilt{ClusterIndex,k}],:),1);
  centersSD(k,:) = std(datanorm4Cluster([clustFilt{ClusterIndex,k}],:),0,1);
end

Ap = [];
for i = 1:k
  % I add a term only to keep the proper structure of the dendrogram according to the neurons in the clusters
  Ap = [Ap;(ones(numElem(i),1)+0.01*((1:numElem(i))'/numElem(i)-0.5))*centers(i,:)];
end
Y = pdist(Ap);
Z = linkage(Y,'average');
c = cluster(Z,'maxclust',k);

% threshold color
t = sort(Z(:,3));
Cth = t(size(Z,1)+2-k);

figure('color','none','visible',visible);
set(gca,'layer','top','color','none')
[h,~,outperm] = dendrogram(Z,0,'ColorThreshold',Cth,'orientation','left');
hold on
set(h,'LineWidth',1)
axis off
ylabel('Neurons','fontsize',16);
xlabel('Distance','fontsize',16);
set(gca,'fontSize',16,'LineWidth',1,'TickDir','out','Box','off','XTick',0:.2:1,'YTick',[])
% plot2svg([resultsdir,'/dendrogram_',type,'_centroids.svg'])
file = [resultsdir,'/',type,num2str(numClust),'ClusteringResults_dendrogram_centroids'];
print(gcf,'-djpeg',file);
print(gcf,'-depsc',file);

% matching the dendrogram with the results of the kmeans to determine the cluster ordering
clustOrder=unique(c(outperm),'stable')';
% clustersNumOfElements = find(diff(c(outperm))~=0);
% clustersNumOfElements = [clustersNumOfElements(1);diff([clustersNumOfElements;length(c)])];
% clustOrder=zeros(1,length(clustersNumOfElements));
% for i = 1:length(clustersNumOfElements)
%   clustOrder(i) = find(clustersNumOfElements(i)==numElem,1,'first');
% end

% save(matfile,'clustFilt','clustOrder','-append');

if isequal(visible,'off'), close; end
