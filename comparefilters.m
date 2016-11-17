%% signal filter tests


clear
close all

addpath('~/fieldtrip-20150712');ft_defaults;

rng(1292493237)
Ntrials = 500;
SampRate = 1000;
NoiseLevel = 1;
TimeVect = -50:150;
LegendString = {'Raw 1000Hz Signal' '200Hz low-pass filter' '100Hz low-pass filter' '30Hz low-pass filter'};

Sig = zeros(Ntrials,length(TimeVect));
Sig(:,101:end) = 1;

figure(1);clf;
Sig1 = Sig  + NoiseLevel*randn(size(Sig));
n=4;t='but';dir='twopass';
Sig2 = ft_preproc_lowpassfilter(Sig1,1000,100,n,t,dir);
Sig3 = ft_preproc_lowpassfilter(Sig1,1000,50,n,t,dir);
Sig5 = ft_preproc_lowpassfilter(Sig1,1000,30,n,t,dir);


figure(1);clf;a=gca;hold on

Sig1mu = mean(Sig1,1);
Sig2mu = mean(Sig2,1);
Sig3mu = mean(Sig3,1);
Sig5mu = mean(Sig5,1);

fastplotMEG(TimeVect,[Sig1mu; Sig2mu; Sig3mu; Sig5mu],a.ColorOrder(1:5,:),LegendString);hold on
plot(TimeVect,0*TimeVect,'k--')
a=gca;
a.XLim=[-50,150];a.YLim=[-.25 1.25];
%title('Underlying signals (same onset)')
legend(LegendString,'Location','NW')

%%
figure(1);clf;
a=gca;
[star_vec] = MakeAShineyPlot(Sig1,TimeVect,a.ColorOrder(1,:),'-',a.ColorOrder(1,:),[-.25 1.25],0,-.1);
[star_vec] = MakeAShineyPlot(Sig2,TimeVect,a.ColorOrder(2,:),'-',a.ColorOrder(2,:),[-.25 1.25],0,-.15);
[star_vec] = MakeAShineyPlot(Sig3,TimeVect,a.ColorOrder(3,:),'-',a.ColorOrder(3,:),[-.25 1.25],0,-.2);
[star_vec] = MakeAShineyPlot(Sig5,TimeVect,a.ColorOrder(5,:),'-',a.ColorOrder(5,:),[-.25 1.25],0,-.25);
a.XLim=[0 100];a.YLim=[-.35 1.5];hold on
plot(TimeVect,0*TimeVect,'k--')
xlabel('time (ms)')
ylabel('Mean amplitude')
%title('Simulated data with significance testing (500 samples, \sigma_{noise} = 1)')
legend(a.Children(end-2:-3:2),LegendString,'Location','NW')

saveas(gcf,'Figures/filtering','png');
saveas(gcf,'Figures/filtering','tif');
saveas(gcf,'Figures/filtering','fig');





