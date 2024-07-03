format short
format compact

%Fixed parameters for Figure 1 of paper
Chains=20;
N=10000;
T=5;
dt=0.05;
p=6;

%Obtain HMC and CHMC sampling data for d=400
d=400;
[qJ0Ch4,qLFCh4]=ChiSampler(d,p,N,Chains,T,dt);
qJ04=reshape(qJ0Ch4,[],1);
qLF4=reshape(qLFCh4,[],1);
%Obtain HMC and CHMC sampling data for d=800
d=800;
[qJ0Ch8,qLFCh8]=ChiSampler(d,p,N,Chains,T,dt);
qJ08=reshape(qJ0Ch8,[],1);
qLF8=reshape(qLFCh8,[],1);
%Obtain HMC and CHMC sampling data for d=1200
d=1200;
[qJ0Ch12,qLFCh12]=ChiSampler(d,p,N,Chains,T,dt);
qJ012=reshape(qJ0Ch12,[],1);
qLF12=reshape(qLFCh12,[],1);

%Plotting
xvalues = {'$1/2$','$1/4$','$1/8$','$1/16$'};
yvalues = {'$1$','$2$','$3$','$4$','$5$'};

alphaLevel = 0.1;

if p==2
    width=3;
end
if p==4
    width=.75;
end
if p>=6
    width=.3;
end

figure(1)
tcl = tiledlayout(3,3,TileSpacing="tight");
histPlots = gobjects(9,1); 


d=400;
colorLF = [0, 0.4470, 0.7410];
colorJ0 = [0.4660, 0.6740, 0.1880];

nexttile(1)
xVal = linspace(0,4,200);
args = (d-1)*log(xVal)-abs(xVal).^p/p - (d/p-1)*log(p)-log(gamma(d/p));
yVal = exp(args);
hold on
plot(xVal,yVal,'linewidth',3,'color','r')
hold off
xlim([0 4])
ylim([0 12])
title("$d=400$",'interpreter','latex');
ylabel('Probability density function','interpreter','latex');

nexttile(4)
ChiSam=AGSamChi(100000,d,p);
xp=(d-1)^(1/p);
dx=(2*width)/100;
edges=xp-width:dx:xp+width;
TF=histfit(ChiSam);
CurveX = TF(2).XData;
CurveY = TF(2).YData;
In=trapz(CurveX,CurveY);

DX=(max(CurveX)-min(CurveX))/2;

histogram(qJ04,edges,'Normalization','pdf','FaceColor',colorJ0)
hold on
histogram(qLF4,edges,'Normalization','pdf','FaceColor',colorLF)
plot(CurveX,CurveY/In,'linewidth',3,'color','r')
hold off
set(gca,'XTick',[])
xlim([xp-0.25 xp+0.25])
ylim([0 12])
ylabel('Histogram of combined samples','interpreter','latex');

nexttile(7)
hold on
for Ch=1:Chains
    samples = qJ0Ch4(Ch,(isfinite(qJ0Ch4(Ch,:))));
    Violin({samples'},Ch,'QuartileStyle','shadow','DataStyle', 'none','ShowMedian', true,'ShowMean', true, 'ViolinColor',{colorJ0});
    samples = qLFCh4(Ch,(isfinite(qLFCh4(Ch,:))));
    if max(samples)-min(samples)>10*eps
        Violin({samples'},10+Ch,'QuartileStyle','shadow','DataStyle', 'none','ShowMedian', true,'ShowMean', true, 'ViolinColor',{colorLF});
    end
end
set(gca,'XTick',[])
hold off
ylim([min(CurveX)-DX/2 max(CurveX)+DX/2])
xlim([0 21])
ylabel({'$\xi$'},'Interpreter','latex')
xlabel('Violin plot of chain samples','interpreter','latex')
view([90,90])

d=800;
colorLF = [0, 0.4470, 0.7410];
colorJ0 = [0.4660, 0.6740, 0.1880];

nexttile(2)
xVal = linspace(0,4,200);
args = (d-1)*log(xVal)-abs(xVal).^p/p - (d/p-1)*log(p)-log(gamma(d/p));
yVal = exp(args);
hold on
plot(xVal,yVal,'linewidth',3,'color','r')
hold off
xlim([0 4])
ylim([0 12])
title("$d=800$",'interpreter','latex');

nexttile(5)
ChiSam=AGSamChi(100000,d,p);
xp=(d-1)^(1/p);
dx=(2*width)/100;
edges=xp-width:dx:xp+width;
TF=histfit(ChiSam);
CurveX = TF(2).XData;
CurveY = TF(2).YData;
In=trapz(CurveX,CurveY);

DX=(max(CurveX)-min(CurveX))/2;

histogram(qJ08,edges,'Normalization','pdf','FaceColor',colorJ0);
hold on
histogram(qLF8,edges,'Normalization','pdf','FaceColor',colorLF)
plot(CurveX,CurveY/In,'linewidth',3,'color','r')
hold off
set(gca,'XTick',[])
set(gca,'YTick',[])
xlim([xp-0.25 xp+0.25])
ylim([0 12])


nexttile(8)
hold on
for Ch=1:Chains
    samples = qJ0Ch8(Ch,(isfinite(qJ0Ch8(Ch,:))));
    Violin({samples'},Ch,'QuartileStyle','shadow','DataStyle', 'none','ShowMedian', true,'ShowMean', true, 'ViolinColor',{colorJ0});
    samples = qLFCh8(Ch,(isfinite(qLFCh8(Ch,:))));
    if max(samples)-min(samples)>10*eps
        Violin({samples'},10+Ch,'QuartileStyle','shadow','DataStyle', 'none','ShowMedian', true,'ShowMean', true, 'ViolinColor',{colorLF});
    end  
end
set(gca,'XTick',[])
hold off
ylim([min(CurveX)-DX/2 max(CurveX)+DX/2])
xlim([0 21])
ylabel({'$\xi$'},'Interpreter','latex')
view([90,90])


d=1200;
colorLF = [0, 0.4470, 0.7410];
colorJ0 = [0.4660, 0.6740, 0.1880];

nexttile(3)
xVal = linspace(0,4,200);
args = (d-1)*log(xVal)-abs(xVal).^p/p - (d/p-1)*log(p)-(log(sqrt(2*pi))+(d/p-1/2)*log(d/p)-d/p);%log(gamma(d/p));
yVal = exp(args);
hold on
plot(xVal,yVal,'linewidth',3,'color','r')
hold off
xlim([0 4])
ylim([0 12])
title("$d=1200$",'interpreter','latex');

nexttile(6)

ChiSam=AGSamChi(100000,d,p);
xp=(d-1)^(1/p);
s=linspace(xp-width,xp+width,1000);
dx=(2*width)/100;
edges=xp-width:dx:xp+width;
TF=histfit(ChiSam);
CurveX = TF(2).XData;
CurveY = TF(2).YData;
In=trapz(CurveX,CurveY);

DX=(max(CurveX)-min(CurveX))/2;
plot(CurveX,CurveY/In,'linewidth',3,'color','r')
hold on
histogram(qJ012,edges,'Normalization','pdf','FaceColor',colorJ0)
histogram(qLF12,edges,'Normalization','pdf','FaceColor',colorLF)
plot(CurveX,CurveY/In,'linewidth',3,'color','r')
hold off
set(gca,'XTick',[])
set(gca,'YTick',[])
legend({'Exact PDF','CHMC','HMC-LF'},'Interpreter','latex','Location','northeast')
xlim([xp-0.25 xp+0.25])
ylim([0 12])
leg = legend('Orientation', 'Horizontal');
leg.Layout.Tile = 'north';

nexttile(9)
hold on
for Ch=1:Chains
    samples = qJ0Ch12(Ch,(isfinite(qJ0Ch12(Ch,:))));
    Violin({samples'},Ch,'QuartileStyle','shadow','DataStyle', 'none','ShowMedian', true,'ShowMean', true, 'ViolinColor',{colorJ0});
    samples = qLFCh12(Ch,(isfinite(qLFCh12(Ch,:))));
    if max(samples)-min(samples)>10*eps
        Violin({samples'},10+Ch,'QuartileStyle','shadow','DataStyle', 'none','ShowMedian', true,'ShowMean', true, 'ViolinColor',{colorLF});
    end
end
set(gca,'XTick',[])
hold off
ylim([min(CurveX)-DX/2 max(CurveX)+DX/2])
xlim([0 21])
ylabel({'$\xi$'},'Interpreter','latex')
view([90,90])

title(tcl,'Comparison of sample distributions', 'interpreter','latex');

compactHistogramViolinPlotStr = strcat('compactPDFHistogramViolinPlot-d',num2str(d),'p',num2str(p),datestr(now,'_dd-mm-yy_HH-MM-SS'));
print('-dpng','-r400',compactHistogramViolinPlotStr)
