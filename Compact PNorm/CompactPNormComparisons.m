%Compact Acceptance
clear all
format short
format compact

%Fixed parameters for Figure 1
Chains=1;
N=10000;
T=4;
p=4;
dimArray=[2560,10240,40960];

%Obtain HMC and CHMC samples for dt=0.1
dt=0.1;
[xiJ0dt1,xiFullJdt1,xiLFdt1,JacJ0dt1,JacFullJdt1,minJ0dt1,minFullJdt1,minLFdt1,deltaHJ0dt1,deltaHFullJdt1,deltaHLFdt1,RejectionTracker1]=PGaussSampler(dimArray,p,N,T,dt);

%Obtain HMC and CHMC samples for dt=0.05
dt=0.05;
[xiJ0dt2,xiFullJdt2,xiLFdt2,JacJ0dt2,JacFullJdt2,minJ0dt2,minFullJdt2,minLFdt2,deltaHJ0dt2,deltaHFullJdt2,deltaHLFdt2,RejectionTracker2]=PGaussSampler(dimArray,p,N,T,dt);

%Obtain HMC and CHMC samples for dt=0.025
dt=0.025;
[xiJ0dt3,xiFullJdt3,xiLFdt3,JacJ0dt3,JacFullJdt3,minJ0dt3,minFullJdt3,minLFdt3,deltaHJ0dt3,deltaHFullJdt3,deltaHLFdt3,RejectionTracker3]=PGaussSampler(dimArray,p,N,T,dt);


%% Generating Plots
fig = figure(1);
clf
tcl = tiledlayout(4,3,TileSpacing="tight");
histPlots = gobjects(9,1); 
warmupN = 10;

% Color scheme for different integrators
colorLF = [0, 0.4470, 0.7410];
colorJ0 = [0.4660, 0.6740, 0.1880];
colorFullJ = [1.0000, 0.7098, 0.1176];

% Color scheme for varying dimensions
colord{1} = [0.3010 0.7250 1.0000];
colord{2} = [0.8500, 0.3250, 0.0980];
colord{3} = [1.0000 0.9137 0.1176];

%% Histogram of exp(delta H)
binLimits = [0,3];
binNumber = 50;

nexttile(tcl)
hold on
for i=1:3
    histogram(exp(deltaHLFdt1(i,warmupN:N)),'Normalization','pdf','NumBins',binNumber,'BinLimits',binLimits,'FaceColor',colord{i})
    if(i==2)
        title("$\tau$ = 0.1",'Interpreter','latex')
    end
end
xlabel({'$\exp(-\Delta H)$'},'Interpreter','latex');
ylabel({'Histogram of $\exp(-\Delta H)$'},'Interpreter','latex')
hold off
yticks([]);

nexttile(tcl)
hold on
for i=1:3
    histogram(exp(deltaHLFdt2(i,warmupN:N)),'Normalization','pdf','NumBins',binNumber,'BinLimits',binLimits,'FaceColor',colord{i})
    if(i==2)
        title("$\tau$ = 0.05",'Interpreter','latex')
    end
end
hold off
xlabel({'$\exp(-\Delta H)$'},'Interpreter','latex');
yticks([]);

nexttile(tcl)
hold on
j=1;
for i=1:3
    histogram(exp(deltaHLFdt3(i,warmupN:N)),'Normalization','pdf','NumBins',binNumber,'BinLimits',binLimits,'FaceColor',colord{i})
    Legends{j}=strcat(' HMC--LF, $d$=', num2str(dimArray(i)));
    j = j+1;
    if(i==2)
        title("$\tau$ = 0.025",'Interpreter','latex')
    end
end
hold off
xlabel({'$\exp(-\Delta H)$'},'Interpreter','latex');
legend(Legends,'Interpreter','latex')
yticks([]);

%% Histogram of det(J)
binLimits=[0,3];
binNumber = 50;

nexttile(tcl)
hold on
for i=1:3
    histogram(JacFullJdt1(i,warmupN:N),'Normalization','pdf','NumBins',binNumber,'BinLimits',binLimits,'FaceColor',colord{i})
end
xlabel({'$\det J_{\Psi_{EP}}$'},'Interpreter','latex');
ylabel({'Histogram of $\det J_{\Psi_{EP}}$'},'Interpreter','latex')
hold off
yticks([]);

nexttile(tcl)
hold on
for i=1:3
    histogram(JacFullJdt2(i,warmupN:N),'Normalization','pdf','NumBins',binNumber,'BinLimits',binLimits,'FaceColor',colord{i})
end
hold off
xlabel({'$\det J_{\Psi_{EP}}$'},'Interpreter','latex');
yticks([]);

nexttile(tcl)
hold on
j=1;
for i=1:3
    histogram(JacFullJdt3(i,warmupN:N),'Normalization','pdf','NumBins',binNumber,'BinLimits',binLimits,'FaceColor',colord{i})
    Legends{j}=strcat(' CHMC--FullJ, $d$=', num2str(dimArray(i)));
    j = j+1;
end
hold off
xlabel({'$\det J_{\Psi_{EP}}$'},'Interpreter','latex');
legend(Legends,'Interpreter','latex')
yticks([]);

%% ViolinPlot for acceptance probability with det(J)
acceptFullJ = cell(3,1);
acceptLF = cell(3,1);
for i=1:3
    acceptFullJ{i} = minFullJdt1(i,warmupN:N);
    acceptLF{i} = minLFdt1(i,warmupN:N);
end

nexttile(tcl)
hold on
for k=1:3
    Violin({acceptLF{k}'},0.75*k,'ShowMedian', false,'ShowMean', true, 'HalfViolin','left','QuartileStyle','shadow','DataStyle', 'none','ViolinColor',{colorLF},'BandWidth', .03);
    Violin({acceptFullJ{k}'},0.75*k,'ShowMedian', false,'ShowMean', true, 'HalfViolin','right','QuartileStyle','shadow','DataStyle', 'none','ViolinColor',{colorFullJ},'BandWidth', .03);
end
set(gca,'XTick',[],'YAxisLocation', 'right')
hold off
ylim([0 1])
xlabel({'$d$'},'Interpreter','latex');
dimStr = dimArray';
xlim([0.25,2.75])
xticks(0.75:0.75:2.25)
xticklabels(dimStr())
ylabel({'Violin plot of $\alpha$'},'Interpreter','latex','Position', [0.1,0.5])

acceptFullJ = cell(3,1);
acceptLF = cell(3,1);
for i=1:3
    acceptFullJ{i} = minFullJdt2(i,warmupN:N);
    acceptLF{i} = minLFdt2(i,warmupN:N);
end

nexttile(tcl)
hold on
for k=1:3
    Violin({acceptLF{k}'},0.75*k,'ShowMedian', false,'ShowMean', true, 'HalfViolin','left','QuartileStyle','shadow','DataStyle', 'none','ViolinColor',{colorLF},'BandWidth', .02);
    Violin({acceptFullJ{k}'},0.75*k,'ShowMedian', false,'ShowMean', true, 'HalfViolin','right','QuartileStyle','shadow','DataStyle', 'none','ViolinColor',{colorFullJ},'BandWidth', .02);
end
set(gca,'XTick',[],'YAxisLocation', 'right')
hold off
ylim([0.35 1])
xlabel({'$d$'},'Interpreter','latex');
dimStr = dimArray';
xlim([0.25,2.75])
xticks(0.75:0.75:2.25)
xticklabels(dimStr)

acceptFullJ = cell(3,1);
acceptLF = cell(3,1);
for i=1:3
    acceptFullJ{i} = minFullJdt3(i,warmupN:N);
    acceptLF{i} = minLFdt3(i,warmupN:N);
end

nexttile(tcl)
hold on
for k=1:3
    Violin({acceptLF{k}'},0.75*k,'ShowMedian', false,'ShowMean', true, 'HalfViolin','left','QuartileStyle','shadow','DataStyle', 'none','ViolinColor',{colorLF},'BandWidth', .01);
    Violin({acceptFullJ{k}'},0.75*k,'ShowMedian', false,'ShowMean', true, 'HalfViolin','right','QuartileStyle','shadow','DataStyle', 'none','ViolinColor',{colorFullJ},'BandWidth', .01);
end
set(gca,'XTick',[],'YAxisLocation', 'right')
hold off
ylim([0.65 1])
xlabel({'$d$'},'Interpreter','latex','FontWeight','bold');
dimStr = dimArray';
xlim([0.25,2.75])
xticks(0.75:0.75:2.25)
xticklabels(dimStr)

%% Plotting pNorm Histogram
alphaLevel = 0.1;
width=.75;

dXi = [0.25, 0.225, 0.2, 0.175, 0.15];
binNumber = 25;

innerTile = tiledlayout(tcl,1,3,TileSpacing="tight");
for i=1:3
    d = dimArray(i);
    xiStar = (d-1)^(1/p);
    binLimits = [xiStar-dXi(i),xiStar+dXi(i)];
    innerTile.Layout.Tile = 10;
    nexttile(innerTile)
    hold on
    histogram(xiFullJdt1(i,warmupN:N),'Normalization','pdf','FaceColor',colorFullJ, 'NumBins',binNumber,'FaceAlpha',0.5,'BinLimits',binLimits)
    histogram(xiJ0dt1(i,warmupN:N),'Normalization','pdf','FaceColor',colorJ0, 'NumBins',binNumber,'FaceAlpha',0.5,'BinLimits',binLimits)
    histogram(xiLFdt1(i,warmupN:N),'Normalization','pdf','FaceColor',colorLF, 'NumBins',binNumber,'FaceAlpha',0.5,'BinLimits',binLimits)
    % Plot pdf of p-Chi distribution
    xVal = linspace(xiStar-0.25,xiStar+0.25,200);
    args = (d-1)*log(xVal)-abs(xVal).^p/p - (d/p-1)*log(p)-gammaln(d/p);
    yVal = exp(args);
    plot(xVal,yVal,'linewidth',2,'color','r')%,'LineStyle','-.')
    hold off
    xiStar = (dimArray(i)-1)^(1/p);
    xStr = strcat('$\xi^*$=', num2str(xiStar));
    xlabel(xStr,'Interpreter','latex');
    xticks([]);
    xlim(binLimits)
    xiStarMax = (dimArray(3)-1)^(1/p);
    args = (dimArray(3)-1)*log(xiStarMax)-abs(xiStarMax).^p/p - (dimArray(3)/p-1)*log(p)-gammaln(dimArray(3)/p);
    ylim([0 exp(args)])
    yticks([]);
    if(i==1)
        ylabel({'Histogram of $\xi$'},'Interpreter','latex')
    end
end

innerTile = tiledlayout(tcl,1,3,TileSpacing="tight");
for i=1:3
    d = dimArray(i);
    xiStar = (d-1)^(1/p);
    binLimits = [xiStar-dXi(i),xiStar+dXi(i)];
    innerTile.Layout.Tile = 11;
    nexttile(innerTile)
    hold on
    histogram(xiFullJdt2(i,warmupN:N),'Normalization','pdf','FaceColor',colorFullJ, 'NumBins',binNumber,'FaceAlpha',0.5,'BinLimits',binLimits)
    histogram(xiJ0dt2(i,warmupN:N),'Normalization','pdf','FaceColor',colorJ0, 'NumBins',binNumber,'FaceAlpha',0.5,'BinLimits',binLimits)
    histogram(xiLFdt2(i,warmupN:N),'Normalization','pdf','FaceColor',colorLF, 'NumBins',binNumber,'FaceAlpha',0.5,'BinLimits',binLimits)
    % Plot pdf of p-Chi distribution
    xVal = linspace(xiStar-0.25,xiStar+0.25,200);
    args = (d-1)*log(xVal)-abs(xVal).^p/p - (d/p-1)*log(p)-gammaln(d/p);
    yVal = exp(args);
    plot(xVal,yVal,'linewidth',2,'color','r')%,'LineStyle','-.')
    hold off
    xiStar = (dimArray(i)-1)^(1/p);
    xStr = strcat('$\xi^*$=', num2str(xiStar));
    xlabel(xStr,'Interpreter','latex');
    xticks([]);
    xlim([xiStar-dXi(i),xiStar+dXi(i)])
    xiStarMax = (dimArray(3)-1)^(1/p);
    args = (dimArray(3)-1)*log(xiStarMax)-abs(xiStarMax).^p/p - (dimArray(3)/p-1)*log(p)-gammaln(dimArray(3)/p);
    ylim([0 exp(args)])
    yticks([]);
end

innerTile = tiledlayout(tcl,1,3,TileSpacing="tight");
for i=1:3
    d = dimArray(i);
    xiStar = (d-1)^(1/p);
    binLimits = [xiStar-dXi(i),xiStar+dXi(i)];
    innerTile.Layout.Tile = 12;
    nexttile(innerTile)
    hold on
    histogram(xiFullJdt3(i,warmupN:N),'Normalization','pdf','FaceColor',colorFullJ, 'NumBins',binNumber,'FaceAlpha',0.5,'BinLimits',binLimits)
    histogram(xiJ0dt3(i,warmupN:N),'Normalization','pdf','FaceColor',colorJ0, 'NumBins',binNumber,'FaceAlpha',0.5,'BinLimits',binLimits)
    histogram(xiLFdt3(i,warmupN:N),'Normalization','pdf','FaceColor',colorLF, 'NumBins',binNumber,'FaceAlpha',0.5,'BinLimits',binLimits)
    % Plot pdf of p-Chi distribution
    xVal = linspace(xiStar-0.25,xiStar+0.25,200);
    args = (d-1)*log(xVal)-abs(xVal).^p/p - (d/p-1)*log(p)-gammaln(d/p);
    yVal = exp(args);
    plot(xVal,yVal,'linewidth',2,'color','r')%,'LineStyle','-.')
    hold off
    xStr = strcat('$\xi^*$=', num2str(xiStar));
    xlabel(xStr,'Interpreter','latex');
    xticks([]);
    xlim([xiStar-dXi(i),xiStar+dXi(i)])
    xiStarMax = (dimArray(3)-1)^(1/p);
    args = (dimArray(3)-1)*log(xiStarMax)-abs(xiStarMax).^p/p - (dimArray(3)/p-1)*log(p)-gammaln(dimArray(3)/p);
    ylim([0 exp(args)])
    yticks([]);
end

ax = axes(tcl,'Visible','off');
hold(ax,'on');
labelExact = plot(ax,NaN,'DisplayName','Exact PDF','color','r','linewidth',2);
labelLF = histogram(ax,NaN,'DisplayName','HMC--LF','FaceColor',colorLF, 'EdgeAlpha',0);
labelCHMC = histogram(ax,NaN,'DisplayName','CHMC','FaceColor',colorJ0, 'EdgeAlpha',0);
labelCHMCFullJ = histogram(ax,NaN,'DisplayName','CHMC--FullJ','FaceColor',colorFullJ, 'EdgeAlpha',0);
hold(ax,'off');
leg = legend([labelExact,labelLF, labelCHMCFullJ, labelCHMC],'Interpreter','latex','Orientation','horizontal','Location','south');
leg.Layout.Tile = 'south'; 

set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
set(findall(gcf,'-property','FontSize'),'FontSize',14)

% title("Comparison of $\alpha$, $\det(J)$ and $\xi$ versus $\tau$", 'interpreter','latex');
compactPNormComparisonPlotStr = strcat('CompactPNormComparison-d',num2str(d),'p',num2str(p),datestr(now,'_dd-mm-yy_HH-MM-SS'));
fig.Position = [0,0,1175,800];
print('-dpng','-r400',compactPNormComparisonPlotStr)
