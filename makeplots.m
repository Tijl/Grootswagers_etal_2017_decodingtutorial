%% This script was used to compare several options along the analysis pipeline for Grootswagers, Wardle, & Carlson
% It uses the massive RESULTS.mat file, created by compare.m
% Some of these plots require 
% fdr.m (http://www-personal.umich.edu/~nichols/FDR/)
% shadederrorbar.m (http://au.mathworks.com/matlabcentral/fileexchange/26311-shadederrorbar)
% the RSA toolbox (http://www.mrc-cbu.cam.ac.uk/methods-and-resources/toolboxes/)
% 
% These scripts are supplied for reference only.
% Tijl Grootswagers, June 2015


clear all
load RESULTS.mat

%% plot all

PLT=struct('field',{'default','dimredmethod','windowsize','samplerate','averaging','classifier','classifiernopca','cvmethod'},...
    ...%legend ordered as in the result struct
    'legend',{  {'Animacy decoding (default pipeline)','Standard error accross subjects','Decoding higher than chance (p<0.05)','Chance','Stimulus on'},...
                {'None','PCA','Anova'},...
                {'1ms','5ms','10ms','20ms'},...
                {'50Hz','100Hz','200Hz','1000Hz'},...
                {'No averaging','Averaging 4 trials','Averaging 8 trials','Averaging 16 trials','Averaging 32 trials'},...
                {'Linear Discriminant (LDA)','Gaussian Naive Bayes (GNB)','Linear SVM','Pearson''s Correlation','Spearman''s Rank Correlation'},...
                {'Linear Discriminant (LDA)','Gaussian Naive Bayes (GNB)','Linear SVM','Pearson''s Correlation','Spearman''s Rank Correlation'},...
                {'No cross-validation','2-fold','10-fold','Leave-one-exemplar-out','Leave-one-trial-out'}},...
    ...%change drawing order for some plots for simplicity
    'order',{[1 2 3 4 5],[2 3 1],[4 3 2 1],[],[5 4 3 2 1],[1 2 3 5 4],[1 2 3 5 4],[1 5 3 2 4]},...    
    ...%some plots need an alternative ylim
    'ylim',{[],[],[],[],[43 80],[44 70],[44 70],[43 80]},...  
    ...%some plots need a title
    'title',{[],[],'B','A',[],'A','B',[]});

for p=1:length(PLT)
    f=figure(p);clf;hold on;a=gca;a.FontSize=16;f.PaperPosition = [0 0 8 6]*1.3;a.LineWidth=2;
    if p==1
        res = [RESULTS.averaging{:,2}]';
    else
        rres = getfield(RESULTS,PLT(p).field);
        res = struct('balancedpcorr',{},'timevect',{});
        for i=1:size(rres,2);
            ss=~cellfun(@isempty,rres(:,i));
            r = [rres{ss,i}];
            res(ss,i)=struct('balancedpcorr',{r.balancedpcorr},'timevect',{r.timevect});
        end
    end
    if isempty(PLT(p).order);order = 1:size(res,2);else order = PLT(p).order;end;
    legendvalues = PLT(p).legend;
    
    if isempty(PLT(p).ylim);ylim([45 70]);else ylim(PLT(p).ylim);end
    hh=[];hp=[];
    for i = 1:size(res,2)
        data = 100*vertcat(res(:,i).balancedpcorr);
        timevect = res(1,order(i)).timevect;
        hh{i} = shadedErrorBar(timevect,data,{@mean,@standarderror},{'color',a.ColorOrder(i,:),'LineWidth',2},1);
        hp{i} = plotstarvect(timevect,data,50,48-order(i)*range(a.YLim)/50,{'color',a.ColorOrder(i,:)});
    end
    xlim(minmax(timevect));
    hc=plot(timevect,50+0*timevect,'k--');hc.LineWidth=1;
    hs=fill([0 66 66 0],[a.YLim(1) a.YLim(1) repmat(a.YLim(1)+range(a.YLim)/50,1,2)],[.5 .5 .5],'EdgeAlpha',1);
    hs.LineWidth=2;
    if p==1;h=[hh{1}.mainLine, hh{1}.patch, hp{1}, hc, hs];else h = [hh{:}];h=[h.mainLine];end;
    
    leg = legend(h(order),legendvalues(order),'Location','NW','box','off');
    
    if ~isempty(PLT(p).title)
        title(PLT(p).title,'HorizontalAlignment','Left','FontSize',40,...
            'Units','Normalized','Position',[-0.25 1.95 0])
    end
    
    xlabel('time (ms)');
    ylabel('classifier accuracy (% correct)');
    drawnow;
    saveas(gcf, sprintf('Figures/%i_%s',p,PLT(p).field),'png')
    saveas(gcf, sprintf('Figures/%i_%s',p,PLT(p).field),'tif')
    saveas(gcf, sprintf('Figures/%i_%s',p,PLT(p).field),'fig')
    
end

%% 3.1 crossdecoding
f=figure(1);clf;a=gca;cm = colormap('jet');colormap(cm);f.PaperPosition = [0 0 8 6]*1.3;
res = [RESULTS.crossdecoding{:}];
data = 100*cat(3,res.pcorr);
p=zeros(size(data,1),size(data,2));
for i=1:size(p,1);
    for j=1:size(p,2);
        p(i,j) = signrank(squeeze(data(i,j,:)),50,'tail','both');
    end
end
% reorder

meandata = squeeze(mean(data,3));

meandata_sig_fdr = nan(size(meandata));
h2 = squeeze(p<fdr(p(:),0.01));
meandata_sig_fdr(h2) = meandata(h2);

h=imagesc(flipud(meandata));axis square
a.FontSize=16;
c=colorbar;c.Label.String='Classifier Accuracy (%)';
xlabel('Training time (ms)')
ylabel('Test time (ms)')
a.XTick=1:20:141;
a.YTick=a.XTick;
a.XTickLabel = res(1).timevect(a.XTick);
a.YTickLabel = res(1).timevect(fliplr(a.YTick));
title('A','HorizontalAlignment','Left','FontSize',40,'Position',[-22 10 0])
saveas(gcf, sprintf('Figures/%i_%s',8,'crossdecoding'),'png')
saveas(gcf, sprintf('Figures/%i_%s',8,'crossdecoding'),'tif')
saveas(gcf, sprintf('Figures/%i_%s',8,'crossdecoding'),'fig')


f=figure(2);clf;a=gca;cm = colormap('jet');colormap([0 1 0;.5 0 0]);f.PaperPosition = [0 0 8 6]*1.3;

h=imagesc(flipud(h2),[.5 1]);axis square
a.FontSize=16;
c=colorbar;c.Label.String='';c.Visible='off';
xlabel('Training time (ms)')
ylabel('Test time (ms)')
a.XTick=1:20:141;
a.YTick=a.XTick;
a.XTickLabel = res(1).timevect(a.XTick);
a.YTickLabel = res(1).timevect(fliplr(a.YTick));
title('B','HorizontalAlignment','Left','FontSize',40,'Position',[-22 10 0])
saveas(gcf, sprintf('Figures/%i_%s',8,'crossdecoding_sig'),'png')
saveas(gcf, sprintf('Figures/%i_%s',8,'crossdecoding_sig'),'tif')
c.Visible='on'; %or matlab crashes :S
saveas(gcf, sprintf('Figures/%i_%s',8,'crossdecoding_sig'),'fig')




%% 2.1 rsa crossdecoding
figure(1);clf;a=gca;cm = [0 0 0; colormap('jet')];colormap(cm);
res = [RESULTS.dsms{:}];
data = cat(3,res.DSM);
groupdata = [];
for s=1:size(data,3)
    groupdata(s,:,:) = 1-pdist2(squeeze(data(:,:,s)),squeeze(mean(data(:,:,1:size(data,3)~=s),3)),'spearman');
end
[h,p] = ttest(groupdata(:,:,:));
meandata = squeeze(mean(groupdata));
meandata_sig = nan(size(meandata));
h1 = squeeze(p<0.01);
meandata_sig(h1) = meandata(h1);
meandata_sig_fdr = nan(size(meandata));
h2 = squeeze(p<fdr(p(:),0.01));
meandata_sig_fdr(h2) = meandata(h2);

hh=imagesc(flipud(meandata_sig_fdr),[0 .3]);axis square
c=colorbar;c.Label.String='Similarity (Spearman''s \rho)';c.TickLabels{1}='N.S.';
xlabel('time (ms)')
ylabel('time (ms)')
a.XTick=1:20:141;
a.YTick=a.XTick;
a.XTickLabel = RESULTS.crossdecoding{1}.timevect(a.XTick);
a.YTickLabel = RESULTS.crossdecoding{1}.timevect(fliplr(a.YTick));
saveas(gcf, sprintf('Figures/%i_%s',8,'crossrsa'),'png')
saveas(gcf, sprintf('Figures/%i_%s',8,'crossrsa'),'fig')

%% 2.1 rsa model testing
animacymodel=vectorizeRDM(categoricalRDM(1:48<25,[],0));
naturalmodel=vectorizeRDM(categoricalRDM([ ones(1,24) 1 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 0],[],0));

stimnames1 = arrayfun(@(x) ['../Experiment1/Pictures/Animate/' x.name],dir('../Experiment1/Pictures/Animate/*.png'), 'UniformOutput', false)';
stimnames2 = arrayfun(@(x) ['../Experiment1/Pictures/Inanimate/' x.name],dir('../Experiment1/Pictures/Inanimate/*.png'), 'UniformOutput', false)';

%load all stimuli
silhouette = zeros(48,512,512);
imagelabels = {};
imagelabels.blackdisks=0;
imagelabels.images = [];
imagelabels.nRows=4;
for i=1:48
    if i<=24
        [im map a]= imread(stimnames1{i});
    else
        [im map a]= imread(stimnames2{i-24});
    end
    silhouette(i,:,:) = a;
    imagex = {};
    imagex.image = im;
    imagex.alpha = a;
    imagelabels.images = [imagex imagelabels.images];
end

jaccardmodel = pdist(silhouette(:,:),'jaccard');

res = [RESULTS.dsms{:}];
data = cat(3,res.DSM);

timevect = RESULTS.crossdecoding{1}.timevect;
groupdata_animacy = [];groupdata_mammal=[];groupdata_jaccard=[];
type = 'Kendall';
for s=1:size(data,3)
    groupdata_animacy(s,:) = corr(animacymodel',data(:,:,s)','type',type);
    groupdata_natural(s,:) = corr(naturalmodel',data(:,:,s)','type',type);
    groupdata_jaccard(s,:) = corr(jaccardmodel',data(:,:,s)','type',type);
end

%ranktransform all data for noise ceilings
for s=1:size(data,3)
    for t=1:size(data,1)
        data(t,:,s) = rankTransform(squeeze(data(t,:,s)));
    end
end

noise_ceiling = [];
for t=1:size(data,1);
    for s=1:size(data,3)
        noise_ceiling(1,t,s) = corr(mean(data(t,:,:),3)',data(t,:,s)','type',type);
        noise_ceiling(2,t,s) = corr(mean(data(t,:,~(s==(1:size(data,3)))),3)',data(t,:,s)','type',type);
    end
end
noise_ceiling=mean(noise_ceiling,3);

%pooled correlations
%bootci(5, @(x) corr(animacymodel',mean(data(:,:,x),3)','type',type), 1:20)
% [groupdata_animacy,groupdata_animacy_p] = corr(animacymodel',mean(data(:,:,:),3)','type',type);
% [groupdata_natural,groupdata_natural_p] = corr(naturalmodel',mean(data(:,:,:),3)','type',type);
% [groupdata_jaccard,groupdata_jaccard_p] = corr(jaccardmodel',mean(data(:,:,:),3)','type',type);

%% plot some RDM's (add niko's toolbox!)
addpath(genpath('~/PHD/Matlab/rsatoolbox'))
timepoints = [-50 100 250 400];
tidx = find(ismember(timevect,timepoints));
res = [RESULTS.dsms{:}];
data = cat(3,res.DSM);
mrdm = mean(data(:,:,:),3);

r = [];
for t=tidx;
    rd={};
    rd.RDM=mrdm(t,:);
    rd.name=sprintf('MEG RDM at %ims',timevect(t));
    r = [r rd];
end

for i=1:4
    close all;
    if i==3
        figure(1);clf;showRDMs(r(i),figure(1),1,[],0,1,imagelabels);
    else
        figure(1);clf;showRDMs(r(i),figure(1),1,[],0,1);
    end
    a=gca;a.Position=[0 0 .95 .95];axis equal;title('')
    f=gcf;f.Renderer='painters';
    set(f,'PaperPosition',2*f.PaperPosition)
    set(f,'PaperSize',2*f.PaperSize)
    saveas(gcf, sprintf('Figures/rdms%i',i),'png')
    saveas(gcf, sprintf('Figures/rdms%i',i),'fig')
end
%% plot the models

r=[];
rd.RDM=jaccardmodel;
rd.name='Silhouette';
r=[r rd];
rd.RDM=animacymodel;
rd.name='Animate - Inanimate';
r=[r rd];
rd.RDM=naturalmodel;
rd.name='Natural - Man-made';
r=[r rd];

for i=1:3
    close all;
    figure(1);clf;showRDMs(r(i),figure(1),1,[0 1],0,1);
    a=gca;a.Position=[0 0 .95 .95];axis equal;title('')
    f=gcf;f.Renderer='painters';
    set(f,'PaperPosition',2*f.PaperPosition)
    set(f,'PaperSize',2*f.PaperSize)
    saveas(gcf, sprintf('Figures/rdm_%s',r(i).name),'png')
    saveas(gcf, sprintf('Figures/rdm_%s',r(i).name),'fig')
end

%% one colourbar
close all;
figure(1);clf;colormap(RDMcolormap());axis off
c=colorbar
c.Ticks=[];
c.Location='North';
c.Label.String = {['\bfDissimilarity'], '[percentile of 1-r]'};
c.Label.Rotation=0;
c.Label.HorizontalAlignment='Left';
c.Label.VerticalAlignment='middle';
c.Label.FontSize=12;
c.Label.Position(1)=2.5;
saveas(gcf, 'Figures/rdm_colourbar','png')
saveas(gcf, 'Figures/rdm_colourbar','fig')

%% plot the rsa model testing
f=figure(1);clf;hold on;a=gca;cm = [0 0 0; colormap('jet')];colormap(cm);a.FontSize=16;f.PaperPosition = [0 0 8 6]*1.3;
h=[];

if 0
    f=patch([timevect,fliplr(timevect)],[noise_ceiling(2,:) fliplr(noise_ceiling(1,:))],[.7 .7 .7]);
    f.FaceAlpha=.3;f.EdgeAlpha=0;f.EdgeColor=f.FaceColor;
else
    f=plot(timevect,noise_ceiling(2,:),':','Color',[.5 .5 .5],'LineWidth',2);
end
    
h1 = shadedErrorBar(timevect,groupdata_animacy,{@mean,@standarderror},{'LineWidth',2,'Color',a.ColorOrder(1,:)},1);
h2 = shadedErrorBar(timevect,groupdata_natural,{@mean,@standarderror},{'LineWidth',2,'Color',a.ColorOrder(2,:)},1);
h3 = shadedErrorBar(timevect,groupdata_jaccard,{@mean,@standarderror},{'LineWidth',2,'Color',a.ColorOrder(3,:)},1);

p1 = plotstarvect(timevect,groupdata_animacy,0,-.018,{'Color',a.ColorOrder(1,:)},'right');
p2 = plotstarvect(timevect,groupdata_natural,0,-.024,{'Color',a.ColorOrder(2,:)},'right');
p3 = plotstarvect(timevect,groupdata_jaccard,0,-.026,{'Color',a.ColorOrder(3,:)},'right');

% h1 = plot(timevect,groupdata_animacy,'LineWidth',2,'Color',a.ColorOrder(1,:))
% h2 = plot(timevect,groupdata_natural,'LineWidth',2,'Color',a.ColorOrder(2,:))
% h3 = plot(timevect,groupdata_jaccard,'LineWidth',2,'Color',a.ColorOrder(3,:))


ylabel({'Mean RDM correlation (Kendall''s tau)'});
xlabel('time (ms)')
xlim([-100 600]);ylim([-.02 .1])
legend([h1.mainLine h2.mainLine h3.mainLine f],{'Animacy model' 'Natural model','Silhouette model','Noise ceiling (lower bound)'},'location','ne','box','off')
plot(timevect,0*timevect,'k--')
saveas(gcf, sprintf('Figures/%i_%s',8,'animacyrsa'),'png')
saveas(gcf, sprintf('Figures/%i_%s',8,'animacyrsa'),'fig')