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
dt=0.1;
dimArray=[1280, 2560, 5120, 10240, 20480, 40960];
%dimArray=[10, 20, 40, 80, 160, 320]; % Smaller dimensions for testing
sizeOfdimArray = length(dimArray);

% Sampling PGauss using HMC/CHMC and compute error metrics
[covMaxErrJ0,covMaxErrLF,KSMaxErrJ0,KSMaxErrLF,W1MaxErrJ0,W1MaxErrLF]=PGaussSampler(dimArray,p,N,T,dt,Chains,numPlotPoints);

%% Generating Plots

%Computing error means across chains
meanCEJ0=zeros(sizeOfdimArray,numPlotPoints);
meanCELF=zeros(sizeOfdimArray,numPlotPoints);
meanKSJ0=zeros(sizeOfdimArray,numPlotPoints);
meanKSLF=zeros(sizeOfdimArray,numPlotPoints);
meanW1J0=zeros(sizeOfdimArray,numPlotPoints);
meanW1LF=zeros(sizeOfdimArray,numPlotPoints);

for dimIdx=1:sizeOfdimArray
    for i=1:numPlotPoints
        meanCEJ0(dimIdx,i)=mean(covMaxErrJ0(dimIdx,i,:));
        meanCELF(dimIdx,i)=mean(covMaxErrLF(dimIdx,i,:));
        meanKSJ0(dimIdx,i)=mean(KSMaxErrJ0(dimIdx,i,:));
        meanKSLF(dimIdx,i)=mean(KSMaxErrLF(dimIdx,i,:));
        meanW1J0(dimIdx,i)=mean(W1MaxErrJ0(dimIdx,i,:));
        meanW1LF(dimIdx,i)=mean(W1MaxErrLF(dimIdx,i,:));
    end
end

% Setting plot colours for Leapfrog and CHMC
colorLF = [0, 0.4470, 0.7410];
colorCHMC = [0.4660, 0.6740, 0.1880];
alphaLevel = 0.1;

fig = figure(1);
clf
tcl = tiledlayout(2,3,TileSpacing="tight",Padding="compact");

for dimIdx = 1:sizeOfdimArray
    nexttile(tcl)
    hold on
    semilogy(1,meanKSLF(dimIdx,:),'color',colorLF,'linewidth',3);
    semilogy(1,meanW1LF(dimIdx,:),'color',colorLF,'linewidth',3,'LineStyle','--');
    semilogy(1,meanCELF(dimIdx,:),'color',colorLF,'linewidth',3,'LineStyle',':');
    
    semilogy(1,meanKSJ0(dimIdx,:),'color',colorCHMC,'linewidth',3);
    semilogy(1,meanW1J0(dimIdx,:),'color',colorCHMC,'linewidth',3,'LineStyle','--');
    semilogy(1,meanCEJ0(dimIdx,:),'color',colorCHMC,'linewidth',3,'LineStyle',':');
    
    for j=1:Chains
        ph = semilogy(plotStepSize:plotStepSize:N,KSMaxErrLF(dimIdx,:,j),'color',colorLF,'linewidth',2);
        ph.Color = [colorLF, alphaLevel];
        ph = semilogy(plotStepSize:plotStepSize:N,KSMaxErrJ0(dimIdx,:,j),'color',colorCHMC,'linewidth',2);
        ph.Color = [colorCHMC, alphaLevel];
        ph = semilogy(plotStepSize:plotStepSize:N,W1MaxErrLF(dimIdx,:,j),'color',colorLF,'linewidth',2,'LineStyle','--');
        ph.Color = [colorLF, alphaLevel];
        ph = semilogy(plotStepSize:plotStepSize:N,W1MaxErrJ0(dimIdx,:,j),'color',colorCHMC,'linewidth',2,'LineStyle','--');
        ph.Color = [colorCHMC, alphaLevel];
        ph = semilogy(plotStepSize:plotStepSize:N,covMaxErrLF(dimIdx,:,j),'color',colorLF,'linewidth',2,'LineStyle',':');
        ph.Color = [colorLF, alphaLevel];
        ph = semilogy(plotStepSize:plotStepSize:N,covMaxErrJ0(dimIdx,:,j),'color',colorCHMC,'linewidth',2,'LineStyle',':');
        ph.Color = [colorCHMC, alphaLevel];
    end
    semilogy(plotStepSize:plotStepSize:N,meanKSLF(dimIdx,:),'color',colorLF,'linewidth',3);
    semilogy(plotStepSize:plotStepSize:N,meanW1LF(dimIdx,:),'color',colorLF,'linewidth',3,'LineStyle','--');
    semilogy(plotStepSize:plotStepSize:N,meanCELF(dimIdx,:),'color',colorLF,'linewidth',3,'LineStyle',':');
    semilogy(plotStepSize:plotStepSize:N,meanKSJ0(dimIdx,:),'color',colorCHMC,'linewidth',3);
    semilogy(plotStepSize:plotStepSize:N,meanW1J0(dimIdx,:),'color',colorCHMC,'linewidth',3,'LineStyle','--');
    semilogy(plotStepSize:plotStepSize:N,meanCEJ0(dimIdx,:),'color',colorCHMC,'linewidth',3,'LineStyle',':');

    title(strcat('$d = $ ',' ',num2str(dimArray(dimIdx))),'Interpreter','latex')
    if( dimIdx <= 3)
        xticklabels({''})
    else
        xlabel({'MCMC Iterations'},'Interpreter','latex')
    end
    if dimIdx == 1 || dimIdx == 4
        yticklabels({''})
    elseif dimIdx == 2 || dimIdx == 5 
        yticklabels({''})
    else
        set(gca,'YAxisLocation', 'right')
    end
    xlim([250 5000])
    ylim([3e-2 1e-0])
    grid on
    set(gca, 'YScale', 'log')
    hold off
end

%Plot legend
ax = axes(tcl,'Visible','off');
hold(ax,'on');
labelLFKS = plot(ax,NaN,'DisplayName','HMC--LF (KS Distance)','color',colorLF,'linewidth',2.5);
labelLFW1 = plot(ax,NaN,'DisplayName','HMC--LF ($W_1$ Distance)','color',colorLF,'linewidth',2.5,'LineStyle','--');
labelLFCov = plot(ax,NaN,'DisplayName','HMC--LF (Covariance)','color',colorLF,'linewidth',2.5,'LineStyle',':');
labelCHMCKS = plot(ax,NaN,'DisplayName','CHMC (KS Distance)','color',colorCHMC,'linewidth',2.5);
labelCHMCW1 = plot(ax,NaN,'DisplayName','CHMC ($W_1$ Distance)','color',colorCHMC,'linewidth',2.5,'LineStyle','--');
labelCHMCCov = plot(ax,NaN,'DisplayName','CHMC (Covariance)','color',colorCHMC,'linewidth',2.5,'LineStyle',':');
hold(ax,'off');
leg = legend([labelLFKS, labelLFW1, labelLFCov, labelCHMCKS, labelCHMCW1, labelCHMCCov],'Interpreter','latex','Orientation','horizontal','Location','south');
leg.Layout.Tile = 'south'; 
ylabel(tcl,{'Error in various metrics (KS Distance, $W_1$ Distance, Covariance)'},'Interpreter','latex')

set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
set(findall(gcf,'-property','FontSize'),'FontSize',14)

compactPGaussPlotStr = strcat('CompactPGaussComparison-p',num2str(p),datestr(now,'_dd-mm-yy_HH-MM-SS'));
fig.Position = [0,0,1125,500];
print('-dpng','-r400',compactPGaussPlotStr)