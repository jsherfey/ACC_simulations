clear all
cd /projectnb/crc-nak/sherfey/projects/ACC_simulations/ACC_models

Ne=80; Ni=.25*Ne;
f1=[0:7.5:60]; f2=[0:7.5:60]; acAMPA1e=4500; enoise=acAMPA1e; 
tauGABA=5;
gAMPAee=[.05]; % E->E
gNMDAee=[.2];
gAMPAei=1;
tspan=[0 5000]; dt=.01; solver='rk1'; compile_flag=1; downsample_factor=10;

acAMPA1e=4500; tauGABA=5;
% acAMPA1e=5500; tauGABA=5; % (0,1) only
% acAMPA1e=7500; tauGABA=5; % (0,1) only
% acAMPA1e=6500; tauGABA=5; % (0,1) only
% acAMPA1e=6500; tauGABA=13; % (0,1) only

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gee={[0 0],[0 1],[1 1]};
  % betweenblockee=0; withinblockee=0;
  % betweenblockee=0; withinblockee=1;
  % betweenblockee=1; withinblockee=1;
for g=1:length(gee)
  betweenblockee=gee{g}(1);
  withinblockee=gee{g}(2);

  hetdegrees=[0 .05 .1 .2 .6 1]; % add: [.4 .8] then plot: [0 .2 .4 .6 .8 1]
  reps=1:5;

  for i=1:length(hetdegrees)
    hetdegree=hetdegrees(i);
    for j=1:length(reps)
      fprintf('g %g/%g, i %g/%g, j %g/%g\n',g,length(gee),i,length(hetdegrees),j,length(reps));
      rep=reps(j);
      %rep=2; hetdegree=.05;
      study_dir_root=sprintf('studies/sweeps_het-vs-hom_ac%g_EEwithin%g-%g_EEbetween%g-%g_gEEnmda%g-%g_gEEampa%g-%g_gEI%g_%g-%gms',acAMPA1e,withinblockee(1),withinblockee(end),betweenblockee(1),betweenblockee(end),gNMDAee(1),gNMDAee(end),gAMPAee(1),gAMPAee(end),gAMPAei,tspan);
      if hetdegree==0
        study_dir=sprintf('%sNe%g_f_%g-%gHz-vs-%g-%gHz_tauI%gms_HOM_%g',study_dir_root,Ne,f1(1),f1(end),f2(1),f2(end),tauGABA,rep);
      else
        study_dir=sprintf('%sNe%g_f_%g-%gHz-vs-%g-%gHz_tauI%gms_HET%g_%g',study_dir_root,Ne,f1(1),f1(end),f2(1),f2(end),tauGABA,hetdegree,rep);
      end
      statfile=[study_dir '_stats.mat'];
      if ~exist(study_dir,'dir') || exist(statfile,'file')
        continue;
      end
%       try
        stats=ImportResults(study_dir,@CalcSpikeSync);
        % reduce size of structure
        for s=1:length(stats)
          stats(s).pairs.Power_SUA1=[]; % note: Power_SUA* is 95% of memory
          stats(s).pairs.Power_SUA2=[];
        end
        % size results for this batch
        save(statfile);
%       end
    end
  end

end

return

data=ImportData(study_dir);
PlotData(data,'plot_type','rastergram');
file=[study_dir '_rastergram']; set(gcf,'PaperPositionMode','auto'); print(gcf,[file '.jpg'],'-djpeg'); print(gcf,[file '.eps'],'-depsc');

% plot results for a single batch
fr=15;
f1=[stats.E_fAMPA1]; uf1=unique(f1); nf1=length(uf1);
f2=[stats.E_fAMPA2]; uf2=unique(f2); nf2=length(uf2);
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

% plot image(f1,f2)
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
