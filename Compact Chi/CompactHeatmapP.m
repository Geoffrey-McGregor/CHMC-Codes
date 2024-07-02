format long
format compact

%Set Fixed Parameters
d=400;
N=10000;
Chains=10;

%Computes HMC and CHMC Errors for p=2
[J0MeanWErrP2,LFMeanWErrP2,J0MeanKSErrP2,LFMeanKSErrP2]=HeatmapGen(d,2,N,Chains);
%Computes HMC and CHMC Errors for p=4
[J0MeanWErrP4,LFMeanWErrP4,J0MeanKSErrP4,LFMeanKSErrP4]=HeatmapGen(d,4,N,Chains);
%Computes HMC and CHMC Errors for p=6
[J0MeanWErrP6,LFMeanWErrP6,J0MeanKSErrP6,LFMeanKSErrP6]=HeatmapGen(d,6,N,Chains);


%Plotting
xvalues = {'$1/2$','$1/4$','$1/8$','$1/16$'};
yvalues = {'$1$','$2$','$3$','$4$','$5$'};

figure(1)
tcl = tiledlayout(2,3,TileSpacing="tight");
KShm = gobjects(6,1); 

nexttile
KShm(1) = heatmap(xvalues,yvalues,LFMeanKSErrP2,'Colormap',summer,'CellLabelColor','none','ColorbarVisible','off');
KShm(1).Interpreter = 'latex';
KShm(1).XDisplayLabels = repmat(' ',4,5);
KShm(1).YLabel = '$T$';
KShm(1).Title = '$p=2$';
clim([0 1])

nexttile
KShm(2) = heatmap(LFMeanKSErrP4,'Colormap',summer,'CellLabelColor','none','ColorbarVisible','off');
KShm(2).Interpreter = 'latex';
KShm(2).XDisplayLabels = repmat(' ',4,5);
KShm(2).YDisplayLabels = repmat(' ',5,4);
KShm(2).Title = '$p=4$';
clim([0 1])

nexttile
KShm(3) = heatmap(LFMeanKSErrP6,'Colormap',summer,'CellLabelColor','none','ColorbarVisible','off');
KShm(3).Interpreter = 'latex';
KShm(3).XDisplayLabels = repmat(' ',4,5);
KShm(3).YDisplayLabels = repmat(' ',5,4);
KShm(3).Title = '$p=6$';
clim([0 1])

nexttile
KShm(4) = heatmap(xvalues,yvalues,J0MeanKSErrP2,'Colormap',summer,'CellLabelColor','none','ColorbarVisible','off');
KShm(4).Interpreter = 'latex';
KShm(4).XLabel = '$\tau$';
KShm(4).YLabel = '$T$';
clim([0 1])

nexttile
KShm(5) = heatmap(xvalues,yvalues,J0MeanKSErrP4,'Colormap',summer,'CellLabelColor','none','ColorbarVisible','off');
KShm(5).Interpreter = 'latex';
KShm(5).YDisplayLabels = repmat(' ',5,4);
KShm(5).XLabel = '$\tau$';
clim([0 1])

nexttile
KShm(6) = heatmap(xvalues,yvalues,J0MeanKSErrP6,'Colormap',summer,'CellLabelColor','none','ColorbarVisible','off');
KShm(6).Interpreter = 'latex';
KShm(6).YDisplayLabels = repmat(' ',5,4);
KShm(6).XLabel = '$\tau$';
clim([0 1])

colorLims = vertcat(KShm.ColorLimits);
globalColorLim = [min(colorLims(:,1)), max(colorLims(:,2))];
set(KShm, 'ColorLimits', globalColorLim)
ax = axes(tcl,'visible','off','Colormap',KShm(1).Colormap,'CLim',globalColorLim);
cb = colorbar(ax);
cb.Layout.Tile = 'East';

tcl.XLabel.Visible='on';
title(tcl,'Error in Kolmogorov-Smirnov distance', 'interpreter','latex');
t1 = text(-0.3,0.35,'HMC-LF','interpreter','latex');
t1.Rotation = 90;
t2 = text(-0.3,-0.75,'CHMC','interpreter','latex');
t2.Rotation = 90;
fontsize(gcf,scale=1.2)
HeatMapPlotStr = strcat('compactHeatMapKS-d',datestr(now,'_dd-mm-yy_HH-MM-SS'));
print('-dpng','-r200',HeatMapPlotStr)

figure(2)
tcl = tiledlayout(2,3,TileSpacing="tight");
Whm = gobjects(6,1); 

nexttile
Whm(1) = heatmap(xvalues,yvalues,LFMeanWErrP2,'Colormap',summer,'CellLabelColor','none','ColorbarVisible','off');
Whm(1).Interpreter = 'latex';
Whm(1).XDisplayLabels = repmat(' ',4,5);
%Whm(1).YDisplayLabels = repmat(' ',5,4);
Whm(1).Title = '$p=2$';
Whm(1).YLabel = '$T$';
clim([0 1])

nexttile
Whm(2) = heatmap(LFMeanWErrP4,'Colormap',summer,'CellLabelColor','none','ColorbarVisible','off');
Whm(2).Interpreter = 'latex';
Whm(2).XDisplayLabels = repmat(' ',4,5);
Whm(2).YDisplayLabels = repmat(' ',5,4);
Whm(2).Title = '$p=4$';
clim([0 1])

nexttile
Whm(3) = heatmap(LFMeanWErrP6,'Colormap',summer,'CellLabelColor','none','ColorbarVisible','off');
Whm(3).Interpreter = 'latex';
Whm(3).XDisplayLabels = repmat(' ',4,5);
Whm(3).YDisplayLabels = repmat(' ',5,4);
Whm(3).Title = '$p=6$';
clim([0 1])

nexttile
Whm(4) = heatmap(xvalues,yvalues,J0MeanWErrP2,'Colormap',summer,'CellLabelColor','none','ColorbarVisible','off');
Whm(4).Interpreter = 'latex';
Whm(4).YLabel = '$T$';
Whm(4).XLabel = '$\tau$';
clim([0 1])

nexttile
Whm(5) = heatmap(xvalues,yvalues,J0MeanWErrP4,'Colormap',summer,'CellLabelColor','none','ColorbarVisible','off');
Whm(5).Interpreter = 'latex';
Whm(5).YDisplayLabels = repmat(' ',5,4);
Whm(5).XLabel = '$\tau$';
clim([0 1])

nexttile
Whm(6) = heatmap(xvalues,yvalues,J0MeanWErrP6,'Colormap',summer,'CellLabelColor','none','ColorbarVisible','off');
Whm(6).Interpreter = 'latex';
Whm(6).YDisplayLabels = repmat(' ',5,4);
Whm(6).XLabel = '$\tau$';
clim([0 1])

colorLims = vertcat(Whm.ColorLimits);
globalColorLim = [min(colorLims(:,1)), max(colorLims(:,2))];
set(Whm, 'ColorLimits', globalColorLim)
ax = axes(tcl,'visible','off','Colormap',Whm(1).Colormap,'CLim',globalColorLim);
cb = colorbar(ax);
cb.Layout.Tile = 'East';


title(tcl,'Error in Wasserstein $W_1$ distance', 'interpreter','latex');
t1 = text(-0.3,0.35,'HMC-LF','interpreter','latex');
t1.Rotation = 90;
t2 = text(-0.3,-0.75,'CHMC','interpreter','latex');
t2.Rotation = 90;
fontsize(gcf,scale=1.2)
HeatMapPlotStr = strcat('compactHeatMapW-d',datestr(now,'_dd-mm-yy_HH-MM-SS'));
print('-dpng','-r200',HeatMapPlotStr)

