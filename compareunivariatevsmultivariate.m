%% make a figure showing separability by combining channels

rng(4)
ntrials = 500;
classA = randn(ntrials,1);
classA = [.5+classA -.5+classA] + .5*randn(ntrials,2);
classB = randn(ntrials,1);
classB = [-.5+classB .5+classB] + .5*randn(ntrials,2);

%% first plot
f=figure(1);f.Position = [100 100 800 800];clf;hold on;
a=gca;a.Position=[0.05 0.05 0.8 0.8];
sa=scatter(classA(:,1),classA(:,2),'o','filled');sa.MarkerFaceColor=sa.CData;
sb=scatter(classB(:,1),classB(:,2),'s','filled');sb.MarkerFaceColor=sb.CData;

legend([sa sb],{'class A','class B'},'Location','NW')

a.XTick=[];a.YTick=[];

xlabel('Voltage channel 1');
ylabel('Voltage channel 2');


xlim([-5 5])
ylim([-5 5])
h=line([-4 4],[-4 4]);h.LineStyle='--';h.Color='k';h.LineWidth=1;

xx=-4:0.1:4;

xa=axes('position',[a.Position(1) a.Position(2) a.Position(3) a.Position(4)/7]);
hold on;axis off
xlim([-5 5])
ylim([0 .6])
plot(xx, normpdf(xx,0.5,1),'LineWidth',2)
plot(xx, normpdf(xx,-0.5,1),'LineWidth',2)

xb=axes('position',[a.Position(1) a.Position(2) a.Position(3)/7 a.Position(4)]);
hold on;axis off
plot(xx, normpdf(xx,0.5,1),'LineWidth',2)
plot(xx, normpdf(xx,-0.5,1),'LineWidth',2)
xlim([-5 5])
ylim([0 .6])
xb.View = [90 90];

xc=axes('position',[a.Position(1)+.65*a.Position(3) a.Position(2)+.65*a.Position(4) a.Position(3) a.Position(4)   ]);
hold on;axis off
plot(xx, normpdf(xx,1,.5),'LineWidth',2)
plot(xx, normpdf(xx,-1,.5),'LineWidth',2)
plot(xc.XLim, 0.*xc.XLim,'k')
xlim([-5 5])
ylim([0 4])
xc.View = [45 90];

set(gcf, 'PaperPosition', [.2 .2 8 8]);
set(gcf, 'Papersize', [11 11]);
saveas(gcf, 'Figures/univmulti','png')
saveas(gcf, 'Figures/univmulti','fig')

%% plot erp components
figure(2)

x = -100:600;
classA=2*sin(x/50+4)+abs(1.05*sin(x/50+4));
classB=2.2*sin(x/50+4)+abs(1.04*sin(x/50+4));

classA(1:225)=.05*sin(x(1:225)/10-4);
classB(1:225)=.05*sin(x(1:225)/10-6);
classA(545:end)=.25*sin(x(545:end)/20+4.9);
classB(545:end)=.25*sin(x(545:end)/20+5.1);
windowSize = 50;
b = (1/windowSize)*ones(1,windowSize);
classA=filter(b,1,classA);
classB=filter(b,1,classB);

ntrials = 500;
classA=repmat(.5*classA,ntrials,1)+.1*randn(ntrials,length(classA));
classB=repmat(.5*classB,ntrials,1)+.1*randn(ntrials,length(classB));

figure(2);clf
subplot(1,2,1);title('ERP observed at channel 1')
hold on;a=gca;xlim(minmax(x));ylim([-.7 1.8]);a.XTick=[];a.YTick=[];a.Box='on';xlabel('time');ylabel('Voltage')
plot(x,0*x,'k--')
ha=shadedErrorBar(x,classB,{@mean, @std},{'color',a.ColorOrder(1,:),'linewidth',2},1)
hb=shadedErrorBar(x,classA,{@mean, @std},{'color',a.ColorOrder(2,:),'linewidth',2},1)
legend([ha.mainLine,hb.mainLine],{'Class A','Class B'},'Location','NW')
p=patch([210,220,220,210],[1.1,1.1,1.75,1.75],[.5 .5 .5]);p.FaceAlpha=.75;p.EdgeAlpha=0;

% subplot(1,2,2);title('ERP observed at channel 2')
% hold on;a=gca;xlim(minmax(x));ylim([-1.8 .8]);a.XTick=[];a.YTick=[];a.Box='on';xlabel('time');ylabel('Voltage')
% plot(x,0*x,'k--')
% hb=shadedErrorBar(x,-classB,{@mean, @std},{'color',a.ColorOrder(1,:),'linewidth',2},1)
% ha=shadedErrorBar(x,-classA,{@mean, @std},{'color',a.ColorOrder(2,:),'linewidth',2},1)
% legend([hb.mainLine,ha.mainLine],{'Class A','Class B'},'Location','NW')
% p=patch([210,220,220,210],[-1.75,-1.75,-1.1,-1.1],[.5 .5 .5]);p.FaceAlpha=.75;p.EdgeAlpha=0;

subplot(1,2,2);title('ERP observed at channel 2')
hold on;a=gca;xlim(minmax(x));ylim([-.7 1.8]);a.XTick=[];a.YTick=[];a.Box='on';xlabel('time');ylabel('Voltage')
plot(x,0*x,'k--')
hb=shadedErrorBar(x,classA,{@mean, @std},{'color',a.ColorOrder(1,:),'linewidth',2},1)
ha=shadedErrorBar(x,classB,{@mean, @std},{'color',a.ColorOrder(2,:),'linewidth',2},1)
legend([hb.mainLine,ha.mainLine],{'Class A','Class B'},'Location','NW')
p=patch([210,220,220,210],[1.1,1.1,1.75,1.75],[.5 .5 .5]);p.FaceAlpha=.75;p.EdgeAlpha=0;

set(gcf, 'PaperPosition', [1 1 10 3]);
set(gcf, 'Papersize', [11 11]);
saveas(gcf, 'Figures/univmulti2','png')
saveas(gcf, 'Figures/univmulti2','fig')

%% new way; plot all in one

rng(4)
ntrials = 500;
classA = randn(ntrials,1);
classA = [.5+classA -.5+classA] + .5*randn(ntrials,2);
classB = randn(ntrials,1);
classB = [-.5+classB .5+classB] + .5*randn(ntrials,2);

% first plot
f=figure(1);clf;f.Position = [-1000 -100 700 1000];
a=axes('Position',[0.05 0.05 0.8571428 0.60]);hold on;a.FontSize=16;
t=title('Combining the activation patterns at one time point');t.HorizontalAlignment='left';
t.Position = [-5 5.05];
text(-5.5,5,'B','FontSize',30,'FontWeight','bold','VerticalAlignment','bottom');

sa=scatter(classA(:,1),classA(:,2),'o','filled');sa.MarkerFaceColor=sa.CData;
sb=scatter(classB(:,1),classB(:,2),'s','filled');sb.MarkerFaceColor=sb.CData;

legend([sa sb],{'Class A','Class B'},'Location','NW')

a.XTick=[];a.YTick=[];

xlabel('Voltage channel 1');
ylabel('Voltage channel 2');

a.XLim=[-5 5];
a.YLim=[-5 5];
h=line([-4 4],[-4 4]);h.LineStyle='--';h.Color='k';h.LineWidth=2;

xx=-4:0.1:4;

xa=axes('position',[a.Position(1) a.Position(2) a.Position(3) a.Position(4)/7]);
hold on;axis off
xlim([-5 5])
ylim([0 .6])
plot(xx, normpdf(xx,0.5,1),'LineWidth',2)
plot(xx, normpdf(xx,-0.5,1),'LineWidth',2)

xb=axes('position',[a.Position(1) a.Position(2) a.Position(3)/7 a.Position(4)]);
hold on;axis off
plot(xx, normpdf(xx,0.5,1),'LineWidth',2)
plot(xx, normpdf(xx,-0.5,1),'LineWidth',2)
xlim([-5 5])
ylim([0 .6])
xb.View = [90 90];

xc=axes('position',[a.Position(1)+.65*a.Position(3) a.Position(2)+.65*a.Position(4) a.Position(3) a.Position(4)   ]);
hold on;axis off
plot(xx, normpdf(xx,1,.5),'LineWidth',2)
plot(xx, normpdf(xx,-1,.5),'LineWidth',2)
plot(xc.XLim, 0.*xc.XLim,'k')
xlim([-5 5])
ylim([0 4])
xc.View = [45 90];

% plot erp components

x = -100:600;
classA=2*sin(x/50+4)+abs(1.05*sin(x/50+4));
classB=2.2*sin(x/50+4)+abs(1.04*sin(x/50+4));

classA(1:225)=.05*sin(x(1:225)/10-4);
classB(1:225)=.05*sin(x(1:225)/10-6);
classA(545:end)=.25*sin(x(545:end)/20+4.9);
classB(545:end)=.25*sin(x(545:end)/20+5.1);
windowSize = 50;
b = (1/windowSize)*ones(1,windowSize);
classA=filter(b,1,classA);
classB=filter(b,1,classB);

ntrials = 500;
classA=repmat(.5*classA,ntrials,1)+.1*randn(ntrials,length(classA));
classB=repmat(.5*classB,ntrials,1)+.1*randn(ntrials,length(classB));

a=axes('Position',[0.05 0.75 0.43 0.2]);a.FontSize=16;
title('ERP observed at channel 1')
hold on;a=gca;xlim(minmax(x));ylim([-.7 1.8]);a.XTick=[];a.YTick=[];a.Box='on';xlabel('time');ylabel('Voltage')
plot(x,0*x,'k--')
ha=shadedErrorBar(x,classB,{@mean, @std},{'color',a.ColorOrder(1,:),'linewidth',2},1);
hb=shadedErrorBar(x,classA,{@mean, @std},{'color',a.ColorOrder(2,:),'linewidth',2},1);
legend([ha.mainLine,hb.mainLine],{'Class A','Class B'},'Location','NW')
p=patch([210,220,220,210],[1.1,1.1,1.75,1.75],[.5 .5 .5]);p.FaceAlpha=.75;p.EdgeAlpha=0;

text(-165,1.8,'A','FontSize',30,'FontWeight','bold','VerticalAlignment','bottom');

a=axes('Position',[0.55 0.75 0.43 0.2]);a.FontSize=16;
title('ERP observed at channel 2')
hold on;xlim(minmax(x));ylim([-.7 1.8]);a.XTick=[];a.YTick=[];a.Box='on';xlabel('time');ylabel('Voltage')
plot(x,0*x,'k--')
hb=shadedErrorBar(x,classA,{@mean, @std},{'color',a.ColorOrder(1,:),'linewidth',2},1);
ha=shadedErrorBar(x,classB,{@mean, @std},{'color',a.ColorOrder(2,:),'linewidth',2},1);
legend([hb.mainLine,ha.mainLine],{'Class A','Class B'},'Location','NW')
p=patch([210,220,220,210],[1.1,1.1,1.75,1.75],[.5 .5 .5]);p.FaceAlpha=.75;p.EdgeAlpha=0;

f.PaperPositionMode='auto';
saveas(gcf,'Figures/univmulti','png')
saveas(gcf,'Figures/univmulti','tif')
saveas(gcf,'Figures/univmulti','fig')




