format long
format compact

%Set Fixed Parameters
p=6;
N=10000;
Chains=10;

%Computes HMC and CHMC Errors for d=400
[J0MeanWErrd400,LFMeanWErrd400,J0MeanKSErrd400,LFMeanKSErrd400]=HeatmapGen(400,p,N,Chains);
%omputes HMC and CHMC Errors for d=800
[J0MeanWErrd800,LFMeanWErrd800,J0MeanKSErrd800,LFMeanKSErrd800]=HeatmapGen(800,p,N,Chains);
%omputes HMC and CHMC Errors for d=1200
[J0MeanWErrd1200,LFMeanWErrd1200,J0MeanKSErrd1200,LFMeanKSErrd1200]=HeatmapGen(1200,p,N,Chains);



%Plotting
xvalues = {'$1/2$','$1/4$','$1/8$','$1/16$'};
yvalues = {'$1$','$2$','$3$','$4$','$5$'};


figure(1)
tcl = tiledlayout(2,3,TileSpacing="tight");
KShm = gobjects(6,1); 

nexttile
KShm(1) = heatmap(xvalues,yvalues,LFMeanKSErrd400,'Colormap',summer,'CellLabelColor','none','ColorbarVisible','off');
KShm(1).Interpreter = 'latex';
KShm(1).XDisplayLabels = repmat(' ',4,5);
KShm(1).YLabel = '$T$';
KShm(1).Title = '$d=400$';
clim([0 1])

nexttile
KShm(2) = heatmap(LFMeanKSErrd800,'Colormap',summer,'CellLabelColor','none','ColorbarVisible','off');
KShm(2).Interpreter = 'latex';
KShm(2).XDisplayLabels = repmat(' ',4,5);
KShm(2).YDisplayLabels = repmat(' ',5,4);
KShm(2).Title = '$d=800$';
clim([0 1])

nexttile
KShm(3) = heatmap(LFMeanKSErrd1200,'Colormap',summer,'CellLabelColor','none','ColorbarVisible','off');
KShm(3).Interpreter = 'latex';
KShm(3).XDisplayLabels = repmat(' ',4,5);
KShm(3).YDisplayLabels = repmat(' ',5,4);
KShm(3).Title = '$d=1200$';
clim([0 1])

nexttile
KShm(4) = heatmap(xvalues,yvalues,J0MeanKSErrd400,'Colormap',summer,'CellLabelColor','none','ColorbarVisible','off');
KShm(4).Interpreter = 'latex';
KShm(4).XLabel = '$\tau$';
KShm(4).YLabel = '$T$';
clim([0 1])

nexttile
KShm(5) = heatmap(xvalues,yvalues,J0MeanKSErrd800,'Colormap',summer,'CellLabelColor','none','ColorbarVisible','off');
KShm(5).Interpreter = 'latex';
KShm(5).YDisplayLabels = repmat(' ',5,4);
KShm(5).XLabel = '$\tau$';
clim([0 1])

nexttile
KShm(6) = heatmap(xvalues,yvalues,J0MeanKSErrd1200,'Colormap',summer,'CellLabelColor','none','ColorbarVisible','off');
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
Whm(1) = heatmap(xvalues,yvalues,LFMeanWErrd400,'Colormap',summer,'CellLabelColor','none','ColorbarVisible','off');
Whm(1).Interpreter = 'latex';
Whm(1).XDisplayLabels = repmat(' ',4,5);
Whm(1).Title = '$d=400$';
Whm(1).YLabel = '$T$';
clim([0 1])

nexttile
Whm(2) = heatmap(LFMeanWErrd800,'Colormap',summer,'CellLabelColor','none','ColorbarVisible','off');
Whm(2).Interpreter = 'latex';
Whm(2).XDisplayLabels = repmat(' ',4,5);
Whm(2).YDisplayLabels = repmat(' ',5,4);
Whm(2).Title = '$d=800$';
clim([0 1])

nexttile
Whm(3) = heatmap(LFMeanWErrd1200,'Colormap',summer,'CellLabelColor','none','ColorbarVisible','off');
Whm(3).Interpreter = 'latex';
Whm(3).XDisplayLabels = repmat(' ',4,5);
Whm(3).YDisplayLabels = repmat(' ',5,4);
Whm(3).Title = '$d=1200$';
clim([0 1])

nexttile
Whm(4) = heatmap(xvalues,yvalues,J0MeanWErrd400,'Colormap',summer,'CellLabelColor','none','ColorbarVisible','off');
Whm(4).Interpreter = 'latex';
Whm(4).YLabel = '$T$';
Whm(4).XLabel = '$\tau$';
clim([0 1])

nexttile
Whm(5) = heatmap(xvalues,yvalues,J0MeanWErrd800,'Colormap',summer,'CellLabelColor','none','ColorbarVisible','off');
Whm(5).Interpreter = 'latex';
Whm(5).YDisplayLabels = repmat(' ',5,4);
Whm(5).XLabel = '$\tau$';
clim([0 1])

nexttile
Whm(6) = heatmap(xvalues,yvalues,J0MeanWErrd1200,'Colormap',summer,'CellLabelColor','none','ColorbarVisible','off');
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

