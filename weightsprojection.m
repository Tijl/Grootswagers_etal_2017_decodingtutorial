subdir = '/Users/43546021/DATA/MEGMethodsPaper/';
files = dir([subdir '*_1000Hz_none.mat']);
nsubjects = length(files);
for i=1:nsubjects
    subjects{i} = files(i).name(19:20);
end

%%
res = [];
for s=1:nsubjects
    %classify and add the activation projections
    [data,B] = loaddata(subjects{s},200,'none');
    [avdata,avlabels] = averagetrials(data.class_dat,B.exemplar,4);
    r = timeseriesdecoding(avdata,avlabels>24,'cvfolds',0,'pcavariance',99,'timevect',data.timevect,'classifier','linear','verbose',2,'weights',1);
    
    r.activationpatterns = r.weights;
    r.correctedactivationpatterns = r.correctedweights;
        
    res{s} = r;
end

%% use fieldtrip for plotting

r = [res{:}];

timevect=r(1).timevect;
timepoints = [-100:100:600]; % timepoints we want to have weights for
addpath('~/fieldtrip-20150712');ft_defaults;

A = cat(3,r.activationpatterns);
A = A(:,:,:);
A = mean(A,3);
A(:) = zscore(A(:));

AC = cat(3,r.correctedactivationpatterns);
AC = AC(:,:,:);
AC = mean(AC,3);
%AC = AC./repmat(max(AC,[],2),1,160);
%AC = AC/max(AC(:));
AC(:) = zscore(AC(:));

ft1 = s1tofieldtrip(A',timevect,[],'timelock');
ft2 = s1tofieldtrip(AC',timevect,[],'timelock');

cfg=ft_read_header('~/DATA/Experiment1/delaytest.con');
ft1.label=cfg.label(1:160);
ft1.grad = cfg.grad;
ft2.label=cfg.label(1:160);
ft2.grad = cfg.grad;
cfg.layout = ft_prepare_layout(cfg);

%
% optional: estimate planar gradients
% cfg.method          = 'triangulation';
% cfg.feedback        = 'no';
% cfg.neighbours      = ft_prepare_neighbours(cfg, ft1);
% cfg.planarmethod    = 'sincos';
% ft1planar        = ft_megplanar(cfg, ft1);
% ft2planar        = ft_megplanar(cfg, ft2);
% 
% %combine planars
% ft1planarcomb = ft_combineplanar(cfg,ft1planar);
% ft2planarcomb = ft_combineplanar(cfg,ft2planar);

%prepare plots
cfg.marker = 'off';
cfg.style = 'straight';   
cfg.colormap = 'jet';
cfg.comment = 'no';
cfg.zlim = [-4 4];

%all in one plot
figure(3);clf;
cfg.colorbar = 'no';
for i=1:length(timepoints)
    
    cfg.xlim = [timepoints(i) timepoints(i)];
    
    subplot(2,length(timepoints),i);
    ft_topoplotER(cfg,ft1);
    title(sprintf('%ims',timepoints(i)))
    
    subplot(2,length(timepoints),length(timepoints)+i);
    ft_topoplotER(cfg,ft2);
    title(sprintf('%ims corrected',timepoints(i)))
    
    drawnow;
end

saveas(figure(3),'Figures/ft_weightprojections_all','fig')
saveas(figure(3),'Figures/ft_weightprojections_all','png')

%individual plots
cfg.colorbar = 'yes';
for i=1:length(timepoints)
    
    figure(1);clf;
    cfg.xlim = [timepoints(i) timepoints(i)];
    ft_topoplotER(cfg,ft1);
    
    figure(2);clf;
    cfg.xlim = [timepoints(i) timepoints(i)];
    ft_topoplotER(cfg,ft2);
    
    for f=1:2
        figure(f);a=gca;a.Position=[0.2 0 .8 1]
        
        c=colorbar;c.Location='westoutside';
        c.Label.String = 'Reconstructed Pattern (A.U.)';
        c.Label.FontSize = 26;
        c.Ticks = [];%[-1:.5:1];
        c.FontSize = 16;
        
    end
    saveas(figure(1),sprintf('Figures/ft_weightprojections_%ims',timepoints(i)),'fig')
    saveas(figure(1),sprintf('Figures/ft_weightprojections_%ims',timepoints(i)),'png')
    saveas(figure(2),sprintf('Figures/ft_corrected_weightprojections_%ims',timepoints(i)),'fig')
    saveas(figure(2),sprintf('Figures/ft_corrected_weightprojections_%ims',timepoints(i)),'png')
end

