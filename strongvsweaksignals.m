
clear
close all
rng(123456)
Ntrials = 500;
SampRate = 199;
NoiseLevel = 1;
TimeVect = -50:149;
LegendString = {'Strong signal' 'Medium signal' 'Weak signal'};

Sig = sin(-pi:2*pi/SampRate:pi);
Sig(1:100) = 0;

plot(Sig);

Sig1 = Sig;
Sig2 = .5*Sig;
Sig3 = .25*Sig;

subplot(2,1,1);
h = plot(TimeVect,[Sig1; Sig2; Sig3]','LineWidth',2);hold on
plot(TimeVect,0*TimeVect,'k--')
xlabel('time(ms)')
ylabel('Signal strength')
a=gca;a.Box='off'
a.XLim=[0 150];a.YLim=[-.25 1.25];
t=title('A:  Underlying Signal');t.Position=[0 1.3 0];t.HorizontalAlignment='left';
l=legend(LegendString)
l.Location='northwest'

Sig1 = repmat(Sig1,[Ntrials 1]) + NoiseLevel*randn([Ntrials length(Sig1)]);
Sig2 = repmat(Sig2,[Ntrials 1]) + NoiseLevel*randn([Ntrials length(Sig2)]);
Sig3 = repmat(Sig3,[Ntrials 1]) + NoiseLevel*randn([Ntrials length(Sig3)]);

Sig1mu = mean(Sig1);
Sig2mu = mean(Sig2);
Sig3mu = mean(Sig3);
% 
% 
% subplot(2,2,2)
% plot(TimeVect,[Sig1mu; Sig2mu; Sig3mu]');
% axis([ -50 150 -.2 1.25])
% title('Average signal (500 samples) sigma = 0.2)')
% legend(LegendString)

%%
subplot(2,1,2)
fastplotMEG(TimeVect,[Sig1mu; Sig2mu; Sig3mu],[h(1).Color;h(2).Color;h(3).Color],LegendString)

%%

subplot(2,1,2);a=gca();
[star_vec] = MakeAShineyPlot(Sig1,TimeVect,h(1).Color,'-',h(1).Color,[-.25 1.25],0,-.1);
[star_vec] = MakeAShineyPlot(Sig2,TimeVect,h(2).Color,'-',h(2).Color,[-.25 1.25],0,-.15);
[star_vec] = MakeAShineyPlot(Sig3,TimeVect,h(3).Color,'-',h(3).Color,[-.25 1.25],0,-.2);
a.XLim=[0 150];a.YLim=[-.25 1.25];hold on
plot(TimeVect,0*TimeVect,'k--');
xlabel('time(ms)')
ylabel('Decoding accuracy')
t=title('B:  Simulated decoding performance');t.Position=[0 1.3 0];t.HorizontalAlignment='left';
l=legend();l.Location='northwest'
saveas(gcf,'Figures/decodingonsets','png');
saveas(gcf,'Figures/decodingonsets','tif');
saveas(gcf,'Figures/decodingonsets','fig');




