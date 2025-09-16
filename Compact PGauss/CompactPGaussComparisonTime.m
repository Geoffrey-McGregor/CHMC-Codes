%Compact Acceptance
clear all
format short
format compact

%Fixed parameters for Figure 1
Chains=10;
N=5000;
numPlotPoints = 20;
plotStepSize = floor(N/numPlotPoints); % Number of points on convergence plot
T=4;
p=4;
dt=0.15;
dimArray=[1280, 2560, 5120, 10240, 20480, 40960];
%dimArray=[10, 20, 40, 80, 160, 640]; % Smaller dimensions for testing
sizeOfdimArray = length(dimArray);

% Sampling PGauss using HMC/CHMC and compute error metrics
[covMaxErrJ0,covMaxErrLF,KSMaxErrJ0,KSMaxErrLF,W1MaxErrJ0,W1MaxErrLF,covMaxErrFJ,KSMaxErrFJ,W1MaxErrFJ,timeLF,timeCHMCJ0,timeFJ]=PGaussSampler(dimArray,p,N,T,dt,Chains,numPlotPoints);

%% Generating Plots

%Computing error means across chains
meanCEJ0=zeros(sizeOfdimArray,numPlotPoints);
meanCELF=zeros(sizeOfdimArray,numPlotPoints);
meanCEFJ=zeros(sizeOfdimArray,numPlotPoints);

meanKSJ0=zeros(sizeOfdimArray,numPlotPoints);
meanKSLF=zeros(sizeOfdimArray,numPlotPoints);
meanKSFJ=zeros(sizeOfdimArray,numPlotPoints);

meanW1J0=zeros(sizeOfdimArray,numPlotPoints);
meanW1LF=zeros(sizeOfdimArray,numPlotPoints);
meanW1FJ=zeros(sizeOfdimArray,numPlotPoints);


for dimIdx=1:sizeOfdimArray
    for i=1:numPlotPoints
        meanCEJ0(dimIdx,i)=mean(covMaxErrJ0(dimIdx,i,:));
        meanCELF(dimIdx,i)=mean(covMaxErrLF(dimIdx,i,:));
        meanCEFJ(dimIdx,i)=mean(covMaxErrFJ(dimIdx,i,:));

        meanKSJ0(dimIdx,i)=mean(KSMaxErrJ0(dimIdx,i,:));
        meanKSLF(dimIdx,i)=mean(KSMaxErrLF(dimIdx,i,:));
        meanKSFJ(dimIdx,i)=mean(KSMaxErrFJ(dimIdx,i,:));

        meanW1J0(dimIdx,i)=mean(W1MaxErrJ0(dimIdx,i,:));
        meanW1LF(dimIdx,i)=mean(W1MaxErrLF(dimIdx,i,:));
        meanW1FJ(dimIdx,i)=mean(W1MaxErrFJ(dimIdx,i,:));

    end
end

% Setting plot colours for Leapfrog and CHMC
colorLF = [0, 0.4470, 0.7410];
colorCHMC = [0.4660, 0.6740, 0.1880];
colorFJ = [1.0000, 0.7098, 0.1176];
alphaLevel = 0.15;

%% Loglog plot
fig = figure(1);
clf
tcl = tiledlayout(2,3,TileSpacing="tight",Padding="compact");

for dimIdx = 1:sizeOfdimArray
    nexttile(tcl)
    hold on
    loglog(1,meanKSLF(dimIdx,:),'color',colorLF,'linewidth',3);
    loglog(1,meanW1LF(dimIdx,:),'color',colorLF,'linewidth',3,'LineStyle','--');
    loglog(1,meanCELF(dimIdx,:),'color',colorLF,'linewidth',3,'LineStyle',':');

    loglog(1,meanKSFJ(dimIdx,:),'color',colorFJ,'linewidth',3);
    loglog(1,meanW1FJ(dimIdx,:),'color',colorFJ,'linewidth',3,'LineStyle','--');
    loglog(1,meanCEFJ(dimIdx,:),'color',colorFJ,'linewidth',3,'LineStyle',':');
    
    loglog(1,meanKSJ0(dimIdx,:),'color',colorCHMC,'linewidth',3);
    loglog(1,meanW1J0(dimIdx,:),'color',colorCHMC,'linewidth',3,'LineStyle','--');
    loglog(1,meanCEJ0(dimIdx,:),'color',colorCHMC,'linewidth',3,'LineStyle',':');
    
    for j=1:Chains
        ph = loglog(timeLF(j,plotStepSize:plotStepSize:N),KSMaxErrLF(dimIdx,:,j),'color',colorLF,'linewidth',2);
        ph.Color = [colorLF, alphaLevel];
        ph = loglog(timeFJ(j,plotStepSize:plotStepSize:N),KSMaxErrFJ(dimIdx,:,j),'color',colorFJ,'linewidth',2);
        ph.Color = [colorFJ, alphaLevel];
        ph = loglog(timeCHMCJ0(j,plotStepSize:plotStepSize:N),KSMaxErrJ0(dimIdx,:,j),'color',colorCHMC,'linewidth',2);
        ph.Color = [colorCHMC, alphaLevel];

        ph = loglog(timeLF(j,plotStepSize:plotStepSize:N),W1MaxErrLF(dimIdx,:,j),'color',colorLF,'linewidth',2,'LineStyle','--');
        ph.Color = [colorLF, alphaLevel];
        ph = loglog(timeFJ(j,plotStepSize:plotStepSize:N),W1MaxErrFJ(dimIdx,:,j),'color',colorFJ,'linewidth',2,'LineStyle','--');
        ph.Color = [colorFJ, alphaLevel];
        ph = loglog(timeCHMCJ0(j,plotStepSize:plotStepSize:N),W1MaxErrJ0(dimIdx,:,j),'color',colorCHMC,'linewidth',2,'LineStyle','--');
        ph.Color = [colorCHMC, alphaLevel];

        ph = loglog(timeLF(j,plotStepSize:plotStepSize:N),covMaxErrLF(dimIdx,:,j),'color',colorLF,'linewidth',2,'LineStyle',':');
        ph.Color = [colorLF, alphaLevel];
        ph = loglog(timeFJ(j,plotStepSize:plotStepSize:N),covMaxErrFJ(dimIdx,:,j),'color',colorFJ,'linewidth',2,'LineStyle',':');
        ph.Color = [colorFJ, alphaLevel];
        ph = loglog(timeCHMCJ0(j,plotStepSize:plotStepSize:N),covMaxErrJ0(dimIdx,:,j),'color',colorCHMC,'linewidth',2,'LineStyle',':');
        ph.Color = [colorCHMC, alphaLevel];
    end
    loglog(timeLF(Chains,plotStepSize:plotStepSize:N),meanKSLF(dimIdx,:),'color',colorLF,'linewidth',3);
    loglog(timeLF(Chains,plotStepSize:plotStepSize:N),meanW1LF(dimIdx,:),'color',colorLF,'linewidth',3,'LineStyle','--');
    loglog(timeLF(Chains,plotStepSize:plotStepSize:N),meanCELF(dimIdx,:),'color',colorLF,'linewidth',3,'LineStyle',':');

    loglog(timeFJ(Chains,plotStepSize:plotStepSize:N),meanKSFJ(dimIdx,:),'color',colorFJ,'linewidth',3);
    loglog(timeFJ(Chains,plotStepSize:plotStepSize:N),meanW1FJ(dimIdx,:),'color',colorFJ,'linewidth',3,'LineStyle','--');
    loglog(timeFJ(Chains,plotStepSize:plotStepSize:N),meanCEFJ(dimIdx,:),'color',colorFJ,'linewidth',3,'LineStyle',':');

    loglog(timeCHMCJ0(Chains,plotStepSize:plotStepSize:N),meanKSJ0(dimIdx,:),'color',colorCHMC,'linewidth',3);
    loglog(timeCHMCJ0(Chains,plotStepSize:plotStepSize:N),meanW1J0(dimIdx,:),'color',colorCHMC,'linewidth',3,'LineStyle','--');
    loglog(timeCHMCJ0(Chains,plotStepSize:plotStepSize:N),meanCEJ0(dimIdx,:),'color',colorCHMC,'linewidth',3,'LineStyle',':');

    title(strcat('$d = $ ',' ',num2str(dimArray(dimIdx))),'Interpreter','latex')
    if( dimIdx <= 3)
        xticklabels({''})
    else
        xlabel({'Running Time (s)'},'Interpreter','latex')
    end
    if dimIdx == 1 || dimIdx == 4
        yticklabels({''})
    elseif dimIdx == 2 || dimIdx == 5 
        yticklabels({''})
    else
        set(gca,'YAxisLocation', 'right')
    end
    %xlim([250 5000])
    xlim([1 16])
    ylim([4e-2 1e-0])
    grid on
    set(gca, 'XScale','log','YScale', 'log')
    hold off
end

%Plot legend
ax = axes(tcl,'Visible','off');
hold(ax,'on');
labelLFKS = plot(ax,NaN,'DisplayName','HMC--LF (KS Distance)','color',colorLF,'linewidth',2.5);
labelFJKS = plot(ax,NaN,'DisplayName','CHMC Full-J (KS Distance)','color',colorFJ,'linewidth',2.5);
labelCHMCKS = plot(ax,NaN,'DisplayName','CHMC (KS Distance)','color',colorCHMC,'linewidth',2.5);

labelLFW1 = plot(ax,NaN,'DisplayName','HMC--LF ($W_1$ Distance)','color',colorLF,'linewidth',2.5,'LineStyle','--');
labelFJW1 = plot(ax,NaN,'DisplayName','CHMC Full-J ($W_1$ Distance)','color',colorFJ,'linewidth',2.5,'LineStyle','--');
labelCHMCW1 = plot(ax,NaN,'DisplayName','CHMC ($W_1$ Distance)','color',colorCHMC,'linewidth',2.5,'LineStyle','--');

labelLFCov = plot(ax,NaN,'DisplayName','HMC--LF (Covariance)','color',colorLF,'linewidth',2.5,'LineStyle',':');
labelFJCov = plot(ax,NaN,'DisplayName','CHMC Full-J (Covariance)','color',colorFJ,'linewidth',2.5,'LineStyle',':');
labelCHMCCov = plot(ax,NaN,'DisplayName','CHMC (Covariance)','color',colorCHMC,'linewidth',2.5,'LineStyle',':');

hold(ax,'off');
leg = legend([labelLFKS, labelFJKS, labelCHMCKS, labelLFW1, labelFJW1, labelCHMCW1, labelLFCov, labelFJCov, labelCHMCCov],'Interpreter','latex','Orientation','horizontal','Location','south');
leg.Layout.Tile = 'south';
leg.NumColumns = 3;
ylabel(tcl,{'Error in various metrics (KS Distance, $W_1$ Distance, Covariance)'},'Interpreter','latex')

set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
set(findall(gcf,'-property','FontSize'),'FontSize',14)

compactPGaussPlotStr = strcat('CompactPGaussComparison-p',num2str(p),datestr(now,'_dd-mm-yy_HH-MM-SS'));
fig.Position = [0,0,1125,500];
print('-dpng','-r400',compactPGaussPlotStr)

%% Semilog plot

fig = figure(2);
clf
tcl = tiledlayout(2,3,TileSpacing="tight",Padding="compact");

for dimIdx = 1:sizeOfdimArray
    currTile = nexttile(tcl);
    hold on
    semilogy(1,meanKSLF(dimIdx,:),'color',colorLF,'linewidth',3);
    semilogy(1,meanW1LF(dimIdx,:),'color',colorLF,'linewidth',3,'LineStyle','--');
    semilogy(1,meanCELF(dimIdx,:),'color',colorLF,'linewidth',3,'LineStyle',':');

    semilogy(1,meanKSFJ(dimIdx,:),'color',colorFJ,'linewidth',3);
    semilogy(1,meanW1FJ(dimIdx,:),'color',colorFJ,'linewidth',3,'LineStyle','--');
    semilogy(1,meanCEFJ(dimIdx,:),'color',colorFJ,'linewidth',3,'LineStyle',':');
    
    semilogy(1,meanKSJ0(dimIdx,:),'color',colorCHMC,'linewidth',3);
    semilogy(1,meanW1J0(dimIdx,:),'color',colorCHMC,'linewidth',3,'LineStyle','--');
    semilogy(1,meanCEJ0(dimIdx,:),'color',colorCHMC,'linewidth',3,'LineStyle',':');
    
    for j=1:Chains
        ph = semilogy(timeLF(j,plotStepSize:plotStepSize:N),KSMaxErrLF(dimIdx,:,j),'color',colorLF,'linewidth',2);
        ph.Color = [colorLF, alphaLevel];
        ph = semilogy(timeFJ(j,plotStepSize:plotStepSize:N),KSMaxErrFJ(dimIdx,:,j),'color',colorFJ,'linewidth',2);
        ph.Color = [colorFJ, alphaLevel];
        ph = semilogy(timeCHMCJ0(j,plotStepSize:plotStepSize:N),KSMaxErrJ0(dimIdx,:,j),'color',colorCHMC,'linewidth',2);
        ph.Color = [colorCHMC, alphaLevel];

        ph = semilogy(timeLF(j,plotStepSize:plotStepSize:N),W1MaxErrLF(dimIdx,:,j),'color',colorLF,'linewidth',2,'LineStyle','--');
        ph.Color = [colorLF, alphaLevel];
        ph = semilogy(timeFJ(j,plotStepSize:plotStepSize:N),W1MaxErrFJ(dimIdx,:,j),'color',colorFJ,'linewidth',2,'LineStyle','--');
        ph.Color = [colorFJ, alphaLevel];
        ph = semilogy(timeCHMCJ0(j,plotStepSize:plotStepSize:N),W1MaxErrJ0(dimIdx,:,j),'color',colorCHMC,'linewidth',2,'LineStyle','--');
        ph.Color = [colorCHMC, alphaLevel];

        ph = semilogy(timeLF(j,plotStepSize:plotStepSize:N),covMaxErrLF(dimIdx,:,j),'color',colorLF,'linewidth',2,'LineStyle',':');
        ph.Color = [colorLF, alphaLevel];
        ph = semilogy(timeFJ(j,plotStepSize:plotStepSize:N),covMaxErrFJ(dimIdx,:,j),'color',colorFJ,'linewidth',2,'LineStyle',':');
        ph.Color = [colorFJ, alphaLevel];
        ph = semilogy(timeCHMCJ0(j,plotStepSize:plotStepSize:N),covMaxErrJ0(dimIdx,:,j),'color',colorCHMC,'linewidth',2,'LineStyle',':');
        ph.Color = [colorCHMC, alphaLevel];
    end
    semilogy(timeLF(Chains,plotStepSize:plotStepSize:N),meanKSLF(dimIdx,:),'color',colorLF,'linewidth',3);
    semilogy(timeLF(Chains,plotStepSize:plotStepSize:N),meanW1LF(dimIdx,:),'color',colorLF,'linewidth',3,'LineStyle','--');
    semilogy(timeLF(Chains,plotStepSize:plotStepSize:N),meanCELF(dimIdx,:),'color',colorLF,'linewidth',3,'LineStyle',':');

    semilogy(timeFJ(Chains,plotStepSize:plotStepSize:N),meanKSFJ(dimIdx,:),'color',colorFJ,'linewidth',3);
    semilogy(timeFJ(Chains,plotStepSize:plotStepSize:N),meanW1FJ(dimIdx,:),'color',colorFJ,'linewidth',3,'LineStyle','--');
    semilogy(timeFJ(Chains,plotStepSize:plotStepSize:N),meanCEFJ(dimIdx,:),'color',colorFJ,'linewidth',3,'LineStyle',':');

    semilogy(timeCHMCJ0(Chains,plotStepSize:plotStepSize:N),meanKSJ0(dimIdx,:),'color',colorCHMC,'linewidth',3);
    semilogy(timeCHMCJ0(Chains,plotStepSize:plotStepSize:N),meanW1J0(dimIdx,:),'color',colorCHMC,'linewidth',3,'LineStyle','--');
    semilogy(timeCHMCJ0(Chains,plotStepSize:plotStepSize:N),meanCEJ0(dimIdx,:),'color',colorCHMC,'linewidth',3,'LineStyle',':');

    title(strcat('$d = $',{' '},num2str(dimArray(dimIdx))),'Interpreter','latex')
    if( dimIdx <= 3)
        xticklabels({''})
        grid(currTile,'minor')
    else
        xlabel({'Running Time (s)'},'Interpreter','latex')
        grid(currTile,'on')
    end
    if dimIdx == 1 || dimIdx == 4
        yticklabels({''})
    elseif dimIdx == 2 || dimIdx == 5 
        yticklabels({''})
    else
        set(gca,'YAxisLocation', 'right')
    end
    xlim([0.5 15.5])
    if( dimIdx <= 3)
        ylim([3.5e-2 1e-0])
    else
        ylim([4e-2 1e-0])
    end
    set(gca,'YScale', 'log')
    hold off
end

%Plot legend
ax = axes(tcl,'Visible','off');
hold(ax,'on');
labelLFKS = plot(ax,NaN,'DisplayName','HMC--LF (KS Distance)','color',colorLF,'linewidth',2.5);
labelFJKS = plot(ax,NaN,'DisplayName','CHMC--FullJ (KS Distance)','color',colorFJ,'linewidth',2.5);
labelCHMCKS = plot(ax,NaN,'DisplayName','CHMC (KS Distance)','color',colorCHMC,'linewidth',2.5);

labelLFW1 = plot(ax,NaN,'DisplayName','HMC--LF ($W_1$ Distance)','color',colorLF,'linewidth',2.5,'LineStyle','--');
labelFJW1 = plot(ax,NaN,'DisplayName','CHMC--FullJ ($W_1$ Distance)','color',colorFJ,'linewidth',2.5,'LineStyle','--');
labelCHMCW1 = plot(ax,NaN,'DisplayName','CHMC ($W_1$ Distance)','color',colorCHMC,'linewidth',2.5,'LineStyle','--');

labelLFCov = plot(ax,NaN,'DisplayName','HMC--LF (Covariance)','color',colorLF,'linewidth',2.5,'LineStyle',':');
labelFJCov = plot(ax,NaN,'DisplayName','CHMC--FullJ (Covariance)','color',colorFJ,'linewidth',2.5,'LineStyle',':');
labelCHMCCov = plot(ax,NaN,'DisplayName','CHMC (Covariance)','color',colorCHMC,'linewidth',2.5,'LineStyle',':');

hold(ax,'off');
leg = legend([labelLFKS, labelFJKS, labelCHMCKS, labelLFW1, labelFJW1, labelCHMCW1, labelLFCov, labelFJCov, labelCHMCCov],'Interpreter','latex','Orientation','horizontal','Location','south');
leg.Layout.Tile = 'south';
leg.NumColumns = 3;
ylabel(tcl,{'Error in various metrics (KS Distance, $W_1$ Distance, Covariance)'},'Interpreter','latex','FontSize',14)

set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
set(findall(gcf,'-property','FontSize'),'FontSize',14)

compactPGaussPlotStr = strcat('CompactPGaussComparisonSemiLog-p',num2str(p),datestr(now,'_dd-mm-yy_HH-MM-SS'));
%fig.Position = [0,0,1125,500];
fig.Position = [0,0,1025,500];
print('-dpng','-r400',compactPGaussPlotStr)