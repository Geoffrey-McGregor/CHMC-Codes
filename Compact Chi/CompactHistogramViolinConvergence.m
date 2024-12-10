format short
format compact

%Fixed parameters for Figure 1
Chains=10;
N=10000;
T=5;
dt=0.05;
p=6;

%Plot every 250th samples
plotStepSize=250;
totalStepNum=length(plotStepSize:plotStepSize:N);

%% Sample p-Chi Distribution using HMC and CHMC for d=400
d=400;
% a=d/p+1-1/p;
% b=p^(1/p-1)*gamma(d/p)/gamma(d/p+1-1/p);
[qJ0Ch4,qLFCh4]=ChiSampler(d,p,N,Chains,T,dt);
qJ04=reshape(qJ0Ch4,[],1);
qLF4=reshape(qLFCh4,[],1);

% Generating exact samples of p-Chi
ChiSam4=AGSamChi(100000,d,p);

%Setting up error variables
J0KSErrCh4=zeros(totalStepNum,Chains);
LFKSErrCh4=zeros(totalStepNum,Chains);
J0KSMeanErr4=zeros(totalStepNum,1);
LFKSMeanErr4=zeros(totalStepNum,1);
J0WErrCh4=zeros(totalStepNum,Chains);
LFWErrCh4=zeros(totalStepNum,Chains);
J0WMeanErr4=zeros(totalStepNum,1);
LFWMeanErr4=zeros(totalStepNum,1);
counter=0;
for j=plotStepSize:plotStepSize:N
    counter=counter+1;
    for i=1:Chains
        %Computing CDF from discrete samples for CHMC and HMC-Leapfrog
        [yJ0,xJ0]=ecdf(qJ0Ch4(i,1:j));
        [yLF,xLF]=ecdf(qLFCh4(i,1:j));

        %Computing Wasserstein distances
        J0WErrCh4(counter,i)=ws_distance(qJ0Ch4(i,1:j)',ChiSam4,1);
        LFWErrCh4(counter,i)=ws_distance(qLFCh4(i,1:j)',ChiSam4,1);

        %Computing Kolmogorov-Smirnov distances
        J0KSErrCh4(counter,i)=max(abs(yJ0-gamcdf(xJ0.^p/p,d/p,1)));
        LFKSErrCh4(counter,i)=max(abs(yLF-gamcdf(xLF.^p/p,d/p,1)));
    end
    %Computing means across chains for each norm
    J0KSMeanErr4(counter)=mean(J0KSErrCh4(counter,:));
    LFKSMeanErr4(counter)=mean(LFKSErrCh4(counter,:));
    J0WMeanErr4(counter)=mean(J0WErrCh4(counter,:));
    LFWMeanErr4(counter)=mean(LFWErrCh4(counter,:));
end

%% Sample p-Chi Distribution using HMC and CHMC for d=800
d=800;
[qJ0Ch8,qLFCh8]=ChiSampler(d,p,N,Chains,T,dt);
qJ08=reshape(qJ0Ch8,[],1);
qLF8=reshape(qLFCh8,[],1);

% Generating exact samples of p-Chi
ChiSam8=AGSamChi(100000,d,p);

%Setting up error variables
J0KSErrCh8=zeros(totalStepNum,Chains);
LFKSErrCh8=zeros(totalStepNum,Chains);
J0KSMeanErr8=zeros(totalStepNum,1);
LFKSMeanErr8=zeros(totalStepNum,1);
J0WErrCh8=zeros(totalStepNum,Chains);
LFWErrCh8=zeros(totalStepNum,Chains);
J0WMeanErr8=zeros(totalStepNum,1);
LFWMeanErr8=zeros(totalStepNum,1);
counter=0;
for j=plotStepSize:plotStepSize:N
    counter=counter+1;
    for i=1:Chains
        %Computing CDF from discrete samples for CHMC and HMC-Leapfrog
        [yJ0,xJ0]=ecdf(qJ0Ch8(i,1:j));
        [yLF,xLF]=ecdf(qLFCh8(i,1:j));

        %Computing Wasserstein distances
        J0WErrCh8(counter,i)=ws_distance(qJ0Ch8(i,1:j)',ChiSam8,1);
        LFWErrCh8(counter,i)=ws_distance(qLFCh8(i,1:j)',ChiSam8,1);

        %Computing Kolmogorov-Smirnov distances
        J0KSErrCh8(counter,i)=max(abs(yJ0-gamcdf(xJ0.^p/p,d/p,1)));
        LFKSErrCh8(counter,i)=max(abs(yLF-gamcdf(xLF.^p/p,d/p,1)));
    end
    %Computing means across chains for each norm
    J0KSMeanErr8(counter)=mean(J0KSErrCh8(counter,:));
    LFKSMeanErr8(counter)=mean(LFKSErrCh8(counter,:));
    J0WMeanErr8(counter)=mean(J0WErrCh8(counter,:));
    LFWMeanErr8(counter)=mean(LFWErrCh8(counter,:));
end

%% Sample p-Chi Distribution using HMC and CHMC for d=1200
d=1200;
[qJ0Ch12,qLFCh12]=ChiSampler(d,p,N,Chains,T,dt);
qJ012=reshape(qJ0Ch12,[],1);
qLF12=reshape(qLFCh12,[],1);

% Generating exact samples of p-Chi
ChiSam12=AGSamChi(100000,d,p);

%Setting up error variables
J0KSErrCh12=zeros(totalStepNum,Chains);
LFKSErrCh12=zeros(totalStepNum,Chains);
J0KSMeanErr12=zeros(totalStepNum,1);
LFKSMeanErr12=zeros(totalStepNum,1);
J0WErrCh12=zeros(totalStepNum,Chains);
LFWErrCh12=zeros(totalStepNum,Chains);
J0WMeanErr12=zeros(totalStepNum,1);
LFWMeanErr12=zeros(totalStepNum,1);
counter=0;
for j=plotStepSize:plotStepSize:N
    counter=counter+1;
    for i=1:Chains
        %Computing CDF from discrete samples for CHMC and HMC-Leapfrog
        [yJ0,xJ0]=ecdf(qJ0Ch12(i,1:j));
        [yLF,xLF]=ecdf(qLFCh12(i,1:j));

        %Computing Wasserstein distances
        J0WErrCh12(counter,i)=ws_distance(qJ0Ch12(i,1:j)',ChiSam12,1);
        LFWErrCh12(counter,i)=ws_distance(qLFCh12(i,1:j)',ChiSam12,1);

        %Computing Kolmogorov-Smirnov distances
        J0KSErrCh12(counter,i)=max(abs(yJ0-gamcdf(xJ0.^p/p,d/p,1)));
        LFKSErrCh12(counter,i)=max(abs(yLF-gamcdf(xLF.^p/p,d/p,1)));
    end
    %Computing means across chains for each norm
    J0KSMeanErr12(counter)=mean(J0KSErrCh12(counter,:));
    LFKSMeanErr12(counter)=mean(LFKSErrCh12(counter,:));
    J0WMeanErr12(counter)=mean(J0WErrCh12(counter,:));
    LFWMeanErr12(counter)=mean(LFWErrCh12(counter,:));
end

%% Plotting
fig = figure(1);
clf
tcl = tiledlayout(4,3,TileSpacing="tight");
alphaLevel = 0.1;
width = 0.3;
colorLF = [0, 0.4470, 0.7410];
colorCHMC = [0.4660, 0.6740, 0.1880];

%% First Column of Plots
nexttile(1)
d=400;
xVal = linspace(0,5,400);
args = (d-1)*log(xVal)-abs(xVal).^p/p - (d/p-1)*log(p)-log(gamma(d/p));
yVal = exp(args);
hold on
plot(xVal,yVal,'linewidth',3,'color','r')
hold off
xiStar = (d-1)^(1/p);
xlim([xiStar-1 xiStar+1])
xiStarMax = (1200-1)^(1/p);
args = (1200-1)*log(xiStarMax)-abs(xiStarMax).^p/p - (1200/p-1)*log(p)-gammaln(1200/p);
ylim([0 exp(args)])
yticks([])
title("$d=400$",'interpreter','latex');
xlabel('$\xi$','Interpreter','latex')
ylabel({'$6$-generalized $\chi$'; 'density function'},'interpreter','latex');

nexttile(4)
xp=(d-1)^(1/p);
dx=(2*width)/100;
edges=xp-width:dx:xp+width;
TF=histfit(ChiSam4);
CurveX = TF(2).XData;
CurveY = TF(2).YData;
In=trapz(CurveX,CurveY);
DX=(max(CurveX)-min(CurveX))/2;
histogram(qJ04,edges,'Normalization','pdf','FaceColor',colorCHMC)
hold on
histogram(qLF4,edges,'Normalization','pdf','FaceColor',colorLF)
plot(CurveX,CurveY/In,'linewidth',3,'color','r')
hold off
set(gca,'XTick',[])
xlim([xp-0.275 xp+0.275])
xticks([round(xp-0.25,2) round(xp,2) round(xp+0.25,2)])
xiStarMax = (1200-1)^(1/p);
args = (1200-1)*log(xiStarMax)-abs(xiStarMax).^p/p - (1200/p-1)*log(p)-gammaln(1200/p);
ylim([0 exp(args)+5])
yticks([]);
xlabel('$\xi$','Interpreter','latex')
ylabel({'Histogram of';'chain samples'},'interpreter','latex');

nexttile(7)
hold on
for Ch=1:Chains
    samples = qJ0Ch4(Ch,(isfinite(qJ0Ch4(Ch,:))));
    Violin({samples'},Ch,'QuartileStyle','shadow','DataStyle', 'none','ShowMedian', true,'ShowMean', true, 'ViolinColor',{colorCHMC});
    samples = qLFCh4(Ch,(isfinite(qLFCh4(Ch,:))));
    if max(samples)-min(samples)>10*eps
        Violin({samples'},10+Ch,'QuartileStyle','shadow','DataStyle', 'none','ShowMedian', true,'ShowMean', true, 'ViolinColor',{colorLF});
    end
end
set(gca,'XTick',[])
hold off
ylim([min(CurveX)-DX/2 max(CurveX)+DX/2])
xlim([0 21])
yticks([]);
xlabel({'Violin plot of';'chain samples'},'interpreter','latex')
xStr = strcat('$\xi^*$=', num2str(xp));
ylabel(xStr,'Interpreter','latex')
view([90,90])

nexttile(10)
hold on
semilogy(1,LFKSMeanErr4(1),'color',colorLF,'linewidth',2);
semilogy(1,LFWMeanErr4(1),'color',colorLF,'linewidth',2,'LineStyle','--');
semilogy(1,J0KSMeanErr4(1),'color',colorCHMC,'linewidth',2);
semilogy(1,J0WMeanErr4(1),'color',colorCHMC,'linewidth',2,'LineStyle','--');
for j=1:Chains
    ph = semilogy(plotStepSize:plotStepSize:N,LFKSErrCh4(:,j),'color',colorLF,'linewidth',2);
    ph.Color = [colorLF, alphaLevel];
    ph = semilogy(plotStepSize:plotStepSize:N,J0KSErrCh4(:,j),'color',colorCHMC,'linewidth',2);
    ph.Color = [colorCHMC, alphaLevel];
    ph = semilogy(plotStepSize:plotStepSize:N,LFWErrCh4(:,j),'color',colorLF,'linewidth',2,'LineStyle','--');
    ph.Color = [colorLF, alphaLevel];
    ph = semilogy(plotStepSize:plotStepSize:N,J0WErrCh4(:,j),'color',colorCHMC,'linewidth',2,'LineStyle','--');
    ph.Color = [colorCHMC, alphaLevel];
end
semilogy(plotStepSize:plotStepSize:N,LFKSMeanErr4,'color',colorLF,'linewidth',3);
semilogy(plotStepSize:plotStepSize:N,LFWMeanErr4,'color',colorLF,'linewidth',3,'LineStyle','--');
semilogy(plotStepSize:plotStepSize:N,J0KSMeanErr4,'color',colorCHMC,'linewidth',3);
semilogy(plotStepSize:plotStepSize:N,J0WMeanErr4,'color',colorCHMC,'linewidth',3,'LineStyle','--');
set(gca, 'YScale', 'log')
hold off
xlim([250 10000])
ylim([10^(-3.2) 10^(0.37)])
grid on

%% Second Column of Plots
nexttile(2)
d=800;
xVal = linspace(0,5,400);
args = (d-1)*log(xVal)-abs(xVal).^p/p - (d/p-1)*log(p)-log(gamma(d/p));
yVal = exp(args);
hold on
plot(xVal,yVal,'linewidth',3,'color','r')
hold off
xiStar = (d-1)^(1/p);
xlim([xiStar-1 xiStar+1])
xiStarMax = (1200-1)^(1/p);
args = (1200-1)*log(xiStarMax)-abs(xiStarMax).^p/p - (1200/p-1)*log(p)-gammaln(1200/p);
ylim([0 exp(args)])
yticks([])
title("$d=800$",'interpreter','latex');
xlabel('$\xi$','Interpreter','latex')

nexttile(5)
xp=(d-1)^(1/p);
dx=(2*width)/100;
edges=xp-width:dx:xp+width;
TF=histfit(ChiSam8);
CurveX = TF(2).XData;
CurveY = TF(2).YData;
In=trapz(CurveX,CurveY);
DX=(max(CurveX)-min(CurveX))/2;
histogram(qJ08,edges,'Normalization','pdf','FaceColor',colorCHMC);
hold on
histogram(qLF8,edges,'Normalization','pdf','FaceColor',colorLF)
plot(CurveX,CurveY/In,'linewidth',3,'color','r')
hold off
set(gca,'XTick',[])
set(gca,'YTick',[])
xlim([xp-0.275 xp+0.275])
xticks([round(xp-0.25,2) round(xp,2) round(xp+0.25,2)])
xiStarMax = (1200-1)^(1/p);
args = (1200-1)*log(xiStarMax)-abs(xiStarMax).^p/p - (1200/p-1)*log(p)-gammaln(1200/p);
ylim([0 exp(args)+5])
yticks([]);
xlabel('$\xi$','Interpreter','latex')

nexttile(8)
hold on
for Ch=1:Chains
    samples = qJ0Ch8(Ch,(isfinite(qJ0Ch8(Ch,:))));
    Violin({samples'},Ch,'QuartileStyle','shadow','DataStyle', 'none','ShowMedian', true,'ShowMean', true, 'ViolinColor',{colorCHMC});
    samples = qLFCh8(Ch,(isfinite(qLFCh8(Ch,:))));
    if max(samples)-min(samples)>10*eps
        Violin({samples'},10+Ch,'QuartileStyle','shadow','DataStyle', 'none','ShowMedian', true,'ShowMean', true, 'ViolinColor',{colorLF});
    end  
end
set(gca,'XTick',[])
hold off
ylim([min(CurveX)-DX/2 max(CurveX)+DX/2])
xlim([0 21])
yticks([]);
xStr = strcat('$\xi^*$=', num2str(xp));
ylabel(xStr,'Interpreter','latex')
view([90,90])

nexttile(11)
hold on
semilogy(1,LFKSMeanErr8(1),'color',colorLF,'linewidth',2);
semilogy(1,LFWMeanErr8(1),'color',colorLF,'linewidth',2,'LineStyle','--');
semilogy(1,J0KSMeanErr8(1),'color',colorCHMC,'linewidth',2);
semilogy(1,J0WMeanErr8(1),'color',colorCHMC,'linewidth',2,'LineStyle','--');
for j=1:Chains
    ph = semilogy(plotStepSize:plotStepSize:N,LFKSErrCh8(:,j),'color',colorLF,'linewidth',2);
    ph.Color = [colorLF, alphaLevel];
    ph = semilogy(plotStepSize:plotStepSize:N,J0KSErrCh8(:,j),'color',colorCHMC,'linewidth',2);
    ph.Color = [colorCHMC, alphaLevel];
    ph = semilogy(plotStepSize:plotStepSize:N,LFWErrCh8(:,j),'color',colorLF,'linewidth',2,'LineStyle','--');
    ph.Color = [colorLF, alphaLevel];
    ph = semilogy(plotStepSize:plotStepSize:N,J0WErrCh8(:,j),'color',colorCHMC,'linewidth',2,'LineStyle','--');
    ph.Color = [colorCHMC, alphaLevel];
end
semilogy(plotStepSize:plotStepSize:N,LFKSMeanErr8,'color',colorLF,'linewidth',3);
semilogy(plotStepSize:plotStepSize:N,LFWMeanErr8,'color',colorLF,'linewidth',3,'LineStyle','--');
semilogy(plotStepSize:plotStepSize:N,J0KSMeanErr8,'color',colorCHMC,'linewidth',3);
semilogy(plotStepSize:plotStepSize:N,J0WMeanErr8,'color',colorCHMC,'linewidth',3,'LineStyle','--');
set(gca, 'YScale', 'log')
hold off
xlim([250 10000])
ylim([10^(-3.2) 10^(0.37)])
yticklabels({})
grid on

%% Third Column of Plots
nexttile(3)
d=1200;
xVal = linspace(0,5,400);
args = (d-1)*log(xVal)-abs(xVal).^p/p - (d/p-1)*log(p)-(log(sqrt(2*pi))+(d/p-1/2)*log(d/p)-d/p);%log(gamma(d/p));
yVal = exp(args);
hold on
plot(xVal,yVal,'linewidth',3,'color','r')
hold off
xiStar = (d-1)^(1/p);
xlim([xiStar-1 xiStar+1])
xiStarMax = (1200-1)^(1/p);
args = (1200-1)*log(xiStarMax)-abs(xiStarMax).^p/p - (1200/p-1)*log(p)-gammaln(1200/p);
ylim([0 exp(args)])
yticks([])
title("$d=1200$",'interpreter','latex');
xlabel('$\xi$','Interpreter','latex')

nexttile(6)
xp=(d-1)^(1/p);
s=linspace(xp-width,xp+width,1000);
dx=(2*width)/100;
edges=xp-width:dx:xp+width;
TF=histfit(ChiSam12);
CurveX = TF(2).XData;
CurveY = TF(2).YData;
In=trapz(CurveX,CurveY);
DX=(max(CurveX)-min(CurveX))/2;
plot(CurveX,CurveY/In,'linewidth',3,'color','r')
hold on
histogram(qLF12,edges,'Normalization','pdf','FaceColor',colorLF,'FaceAlpha',0.5)
histogram(qJ012,edges,'Normalization','pdf','FaceColor',colorCHMC)
histogram(qLF12,edges,'Normalization','pdf','FaceColor',colorLF,'FaceAlpha',0.5)
plot(CurveX,CurveY/In,'linewidth',3,'color','r')
hold off
set(gca,'XTick',[])
set(gca,'YTick',[])
legend({'Exact PDF','HMC-LF','CHMC'},'Interpreter','latex')
xlim([xp-0.275 xp+0.275])
xticks([round(xp-0.25,2) round(xp,2) round(xp+0.25,2)])
xiStarMax = (1200-1)^(1/p);
args = (1200-1)*log(xiStarMax)-abs(xiStarMax).^p/p - (1200/p-1)*log(p)-gammaln(1200/p);
ylim([0 exp(args)+5])
yticks([]);
xlabel('$\xi$','Interpreter','latex')
leg = legend('Orientation', 'Horizontal');
leg.Layout.Tile = 'north';

nexttile(9)
hold on
for Ch=1:Chains
    samples = qJ0Ch12(Ch,(isfinite(qJ0Ch12(Ch,:))));
    Violin({samples'},Ch,'QuartileStyle','shadow','DataStyle', 'none','ShowMedian', true,'ShowMean', true, 'ViolinColor',{colorCHMC});
    samples = qLFCh12(Ch,(isfinite(qLFCh12(Ch,:))));
    if max(samples)-min(samples)>10*eps
        Violin({samples'},10+Ch,'QuartileStyle','shadow','DataStyle', 'none','ShowMedian', true,'ShowMean', true, 'ViolinColor',{colorLF});
    end
end
set(gca,'XTick',[])
hold off
ylim([min(CurveX)-DX/2 max(CurveX)+DX/2])
xlim([0 21])
yticks([]);
xStr = strcat('$\xi^*$=', num2str(xp));
ylabel(xStr,'Interpreter','latex')
view([90,90])

nexttile(12)
hold on
semilogy(1,LFKSMeanErr12(1),'color',colorLF,'linewidth',2);
semilogy(1,LFWMeanErr12(1),'color',colorLF,'linewidth',2,'LineStyle','--');
semilogy(1,J0KSMeanErr12(1),'color',colorCHMC,'linewidth',2);
semilogy(1,J0WMeanErr12(1),'color',colorCHMC,'linewidth',2,'LineStyle','--');
for j=1:Chains
    ph = semilogy(plotStepSize:plotStepSize:N,LFKSErrCh12(:,j),'color',colorLF,'linewidth',2);
    ph.Color = [colorLF, alphaLevel];
    ph = semilogy(plotStepSize:plotStepSize:N,J0KSErrCh12(:,j),'color',colorCHMC,'linewidth',2);
    ph.Color = [colorCHMC, alphaLevel];
    ph = semilogy(plotStepSize:plotStepSize:N,LFWErrCh12(:,j),'color',colorLF,'linewidth',2,'LineStyle','--');
    ph.Color = [colorLF, alphaLevel];
    ph = semilogy(plotStepSize:plotStepSize:N,J0WErrCh12(:,j),'color',colorCHMC,'linewidth',2,'LineStyle','--');
    ph.Color = [colorCHMC, alphaLevel];
end
semilogy(plotStepSize:plotStepSize:N,LFKSMeanErr12,'color',colorLF,'linewidth',3);
semilogy(plotStepSize:plotStepSize:N,LFWMeanErr12,'color',colorLF,'linewidth',3,'LineStyle','--');
semilogy(plotStepSize:plotStepSize:N,J0KSMeanErr12,'color',colorCHMC,'linewidth',3);
semilogy(plotStepSize:plotStepSize:N,J0WMeanErr12,'color',colorCHMC,'linewidth',3,'LineStyle','--');
set(gca, 'YScale', 'log')
hold off
xlim([250 10000])
ylim([10^(-3.2) 10^(0.37)])
yticklabels({})
grid on

ax = axes(tcl,'Visible','off');
hold(ax,'on');
labelLFKS = plot(ax,NaN,'DisplayName','HMC-LF (KS Distance)','color',colorLF,'linewidth',2.5);
labelLFW1 = plot(ax,NaN,'DisplayName','HMC-LF ($W_1$ Distance)','color',colorLF,'linewidth',2.5,'LineStyle','--');
labelCHMCKS = plot(ax,NaN,'DisplayName','CHMC (KS Distance)','color',colorCHMC,'linewidth',2.5);
labelCHMCW1 = plot(ax,NaN,'DisplayName','CHMC ($W_1$ Distance)','color',colorCHMC,'linewidth',2.5,'LineStyle','--');
hold(ax,'off');
leg = legend([labelLFKS, labelLFW1, labelCHMCKS, labelCHMCW1],'Interpreter','latex','Orientation','horizontal','Location','south');
leg.Layout.Tile = 'south'; 

set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
set(findall(gcf,'-property','FontSize'),'FontSize',14)

compactHistogramViolinPlotStr = strcat('compactHistogramViolinConvergencePlot-p',num2str(p),datestr(now,'_dd-mm-yy_HH-MM-SS'));
fig.Position = [0,0,1175,825];
print('-dpng','-r400',compactHistogramViolinPlotStr)
