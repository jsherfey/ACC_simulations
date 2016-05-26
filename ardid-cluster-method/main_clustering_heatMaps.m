% savesvg=0; if savesvg==0, visible='on'; else visible='on'; end

% clear;
% close all;
% clc;

numClust=numClusters(ClusterIndex);
% type = 'Kmeans';
% resultsdir = 'clusteringResults';
matfile = [resultsdir,'/',type,num2str(numClust),'ClusteringResults'];
% load(matfile);
if ~exist(resultsdir,'dir'), mkdir(resultsdir); end
[numNeurons,numMeasures] = size(datanorm4Cluster);

% figure('color','none','visible',visible)
% hold on
% set(gca,'layer','top','color','none')
% p=mypcolor(1:numMeasures,1:numNeurons,datanorm4Cluster);
% rgb_red=[1,0.35,0.35]; % red
% rgb_blue=[0.35,0.35,1]; % blue
% map = diverging_map(0:1/255:1,rgb_red,rgb_blue);
% colormap(map)
% set(p,'EdgeColor','interp');
% axis tight
% axis off
% filename=[resultsdir,'/rawHeatMap'];
% svgFile=[filename,'.svg'];
% plot2svg(svgFile);
figure('visible',visible)
hold on
imagesc(1:numMeasures,1:numNeurons,datanorm4Cluster);
rgb_red=[1,0.35,0.35]; % red
rgb_blue=[0.35,0.35,1]; % blue
map = diverging_map(0:1/255:1,rgb_red,rgb_blue);
colormap(map); axis tight; %axis off; 
file = [resultsdir,'/',type,num2str(numClust),'ClusteringResults_rawHeatMap'];
print(gcf,'-djpeg',file);
print(gcf,'-depsc',file);

% resorting data according to the dendrogram of the k-means clustering
% numClusters=7;
% offset=6;%4; % it starts with numClusters=5
% ClusterIndex=numClusters-offset;

k=1;
boundaries=[];
sortedData=zeros(size(datanorm4Cluster));
numCellType = nan(size(celltype));
numCellType(strcmp('narrow',celltype)) = 1;
numCellType(strcmp('broad',celltype)) = 2;
numCellType(strcmp('fuzzy',celltype)) = 1.5;
sortedCelltype=zeros(size(numCellType));
for i=clustOrder
  these_ids=[clustFilt{ClusterIndex,i}];
    if(~isempty(these_ids))
        tmp = datanorm4Cluster(these_ids,:);
        [tmp2,ind]=sort(numCellType(these_ids));
        sortedCellTypes(k:k+length(these_ids)-1)=tmp2;
        sortedData(k:k+length(these_ids)-1,:)=tmp(ind,:);
        k=k+length(these_ids);
        boundaries=[boundaries;k];
    end
end
boundaries(end)=[];

if 0
  figure('color','none','visible',visible)
  hold on
  set(gca,'layer','top','color','none')
  p=mypcolor(1:numMeasures,1:numNeurons,sortedData);
  colormap(map);
  colorbar
  set(p,'EdgeColor','interp');
  linesx=ones(length(clustOrder)-1,1)*[1 numMeasures+1];
  linesy=[boundaries boundaries];
  plot(linesx',linesy','k-','LineWidth',2)
  axis tight
  axis off
  % filename=[resultsdir,'/colorMap'];
  % svgFile=[filename,'.svg'];
  % plot2svg(svgFile);
  file = [resultsdir,'/',type,num2str(numClust),'ClusteringResults_colorMap'];
  print(gcf,'-djpeg',file);
  print(gcf,'-depsc',file);
end

% figure('visible',visible)
% hold on
% set(gca,'layer','top','color','none')
% p=mypcolor(1:numMeasures,1:numNeurons,sortedData);
% colormap(map);
% set(p,'EdgeColor','interp');
% linesx=ones(length(clustOrder)-1,1)*[1 numMeasures+1];
% linesy=[boundaries boundaries];
% plot(linesx',linesy','k:','LineWidth',1)
% axis tight
% axis off
% filename=[resultsdir,'/sortedHeatMap'];
% svgFile=[filename,'.svg'];
% plot2svg(svgFile);

figure('visible',visible)
hold on
imagesc(1:numMeasures,1:numNeurons,sortedData);
colormap(map); axis tight; %axis off; 
try for i=1:length(boundaries), hline(boundaries(i),'color','k','linewidth',3); end; end
title(['boundaries=[' num2str(boundaries') ']'])
set(gca,'xticklabel',IPlabels); colorbar
file = [resultsdir,'/',type,num2str(numClust),'ClusteringResults_sortedHeatMap'];
print(gcf,'-djpeg',file);
print(gcf,'-depsc',file);

return

figure('color','none','visible',visible)
hexcolors={'1D0091','007BFF','03FC5E','FFFF00','FFAE00','FF4000','000000','A8100D'};
colororder=zeros(length(hexcolors),3);
for i=1:length(hexcolors)
    colororder(i,:) = rgbconv(hexcolors{i});
end
set(gca, 'ColorOrder', colororder);
hold on
set(gca,'layer','top','color','none')
linesx=ones(length(clustOrder),2);
linesy=[1 boundaries' ;boundaries' length(numCellType)]';
plot(linesx',linesy','LineWidth',5)
axis tight
axis off
% filename=[resultsdir,'/clusterColors'];
% svgFile=[filename,'.svg'];
% plot2svg(svgFile);
file = [resultsdir,'/',type,num2str(numClust),'ClusteringResults_clusterColors'];
print(gcf,'-djpeg',file);
print(gcf,'-depsc',file);

figure('color','none','visible',visible)
hold on
set(gca,'layer','top','color','none')
p=mypcolor(1,1:length(numCellType),sortedCellTypes');
map=[
  1 0 0
  .5 .5 .5
  0 0 1
];
colormap(map);
axis tight
xlim([0 7.5]);
axis off
% filename=[resultsdir,'/cellTypeColors'];
% svgFile=[filename,'.svg'];
% plot2svg(svgFile);
file = [resultsdir,'/',type,num2str(numClust),'ClusteringResults_cellTypeColors'];
print(gcf,'-djpeg',file);
print(gcf,'-depsc',file);

save(matfile,'numNeurons','numMeasures','numCellType','clustOrder','hexcolors','-append');

if isequal(visible,'off'), close all; end
