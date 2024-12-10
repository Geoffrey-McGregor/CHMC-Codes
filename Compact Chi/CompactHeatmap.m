format long
format compact

%% Generate Heatmaps with varying d
p=6;
N=10000;
Chains=10;

%Calls HeatMapGen to computes HMC and CHMC Errors for d=400
[J0MeanWErrd400,LFMeanWErrd400,J0MeanKSErrd400,LFMeanKSErrd400]=HeatmapGen(400,p,N,Chains);
%Calls HeatMapGen to computes HMC and CHMC Errors for d=800
[J0MeanWErrd800,LFMeanWErrd800,J0MeanKSErrd800,LFMeanKSErrd800]=HeatmapGen(800,p,N,Chains);
%Calls HeatMapGen to computes HMC and CHMC Errors for d=1200
[J0MeanWErrd1200,LFMeanWErrd1200,J0MeanKSErrd1200,LFMeanKSErrd1200]=HeatmapGen(1200,p,N,Chains);

%% Generate Heatmaps with varying p
d=400;
%Calls HeatMapGen to computes HMC and CHMC Errors for p=2
[J0MeanWErrP2,LFMeanWErrP2,J0MeanKSErrP2,LFMeanKSErrP2]=HeatmapGen(d,2,N,Chains);
%Calls HeatMapGen to computes HMC and CHMC Errors for p=4
[J0MeanWErrP4,LFMeanWErrP4,J0MeanKSErrP4,LFMeanKSErrP4]=HeatmapGen(d,4,N,Chains);
%Calls HeatMapGen to computes HMC and CHMC Errors for p=6
[J0MeanWErrP6,LFMeanWErrP6,J0MeanKSErrP6,LFMeanKSErrP6]=HeatmapGen(d,6,N,Chains);

%% Plots Heatmap error in KS distance
xvalues = {'$2^{-1}$','$2^{-2}$','$2^{-3}$','$2^{-4}$'};
yvalues = {'$1$','$2$','$3$','$4$','$5$'};

fig = figure(1);
clf
tcl = tiledlayout(2,6,TileSpacing="tight");
KShm = gobjects(6,1); 

%% First row of heatmap (HMC)
nexttile
KShm(1) = heatmap(xvalues,yvalues,LFMeanKSErrd400,'Colormap',sky,'CellLabelColor','none','ColorbarVisible','off');
KShm(1).Interpreter = 'latex';
KShm(1).XDisplayLabels = repmat(' ',4,5);
KShm(1).YLabel = '$T$';
KShm(1).Title = '$d=400$';
clim([0 1])

nexttile
KShm(2) = heatmap(LFMeanKSErrd800,'Colormap',sky,'CellLabelColor','none','ColorbarVisible','off');
KShm(2).Interpreter = 'latex';
KShm(2).XDisplayLabels = repmat(' ',4,5);
KShm(2).YDisplayLabels = repmat(' ',5,4);
KShm(2).Title = '$d=800$';
clim([0 1])

nexttile
KShm(3) = heatmap(LFMeanKSErrd1200,'Colormap',sky,'CellLabelColor','none','ColorbarVisible','off');
KShm(3).Interpreter = 'latex';
KShm(3).XDisplayLabels = repmat(' ',4,5);
KShm(3).YDisplayLabels = repmat(' ',5,4);
KShm(3).Title = '$d=1200$';
clim([0 1])

nexttile
KShm(4) = heatmap(xvalues,yvalues,LFMeanKSErrP2,'Colormap',sky,'CellLabelColor','none','ColorbarVisible','off');
KShm(4).Interpreter = 'latex';
KShm(4).XDisplayLabels = repmat(' ',4,5);
KShm(4).YDisplayLabels = repmat(' ',5,4);
KShm(4).Title = '$p=2$';
clim([0 1])

nexttile
KShm(5) = heatmap(LFMeanKSErrP4,'Colormap',sky,'CellLabelColor','none','ColorbarVisible','off');
KShm(5).Interpreter = 'latex';
KShm(5).XDisplayLabels = repmat(' ',4,5);
KShm(5).YDisplayLabels = repmat(' ',5,4);
KShm(5).Title = '$p=4$';
clim([0 1])

nexttile
KShm(6) = heatmap(LFMeanKSErrP6,'Colormap',sky,'CellLabelColor','none','ColorbarVisible','off');
KShm(6).Interpreter = 'latex';
KShm(6).XDisplayLabels = repmat(' ',4,5);
KShm(6).YDisplayLabels = repmat(' ',5,4);
KShm(6).Title = '$p=6$';
clim([0 1])

%% Second row of heatmap (CHMC)
nexttile
KShm(7) = heatmap(xvalues,yvalues,J0MeanKSErrd400,'Colormap',sky,'CellLabelColor','none','ColorbarVisible','off');
KShm(7).Interpreter = 'latex';
KShm(7).XLabel = '$\tau$';
KShm(7).YLabel = '$T$';
clim([0 1])

nexttile
KShm(8) = heatmap(xvalues,yvalues,J0MeanKSErrd800,'Colormap',sky,'CellLabelColor','none','ColorbarVisible','off');
KShm(8).Interpreter = 'latex';
KShm(8).YDisplayLabels = repmat(' ',5,4);
KShm(8).XLabel = '$\tau$';
clim([0 1])

nexttile
KShm(9) = heatmap(xvalues,yvalues,J0MeanKSErrd1200,'Colormap',sky,'CellLabelColor','none','ColorbarVisible','off');
KShm(9).Interpreter = 'latex';
KShm(9).YDisplayLabels = repmat(' ',5,4);
KShm(9).XLabel = '$\tau$';
clim([0 1])

nexttile
KShm(10) = heatmap(xvalues,yvalues,J0MeanKSErrP2,'Colormap',sky,'CellLabelColor','none','ColorbarVisible','off');
KShm(10).Interpreter = 'latex';
KShm(10).YDisplayLabels = repmat(' ',5,4);
KShm(10).XLabel = '$\tau$';
clim([0 1])

nexttile
KShm(11) = heatmap(xvalues,yvalues,J0MeanKSErrP4,'Colormap',sky,'CellLabelColor','none','ColorbarVisible','off');
KShm(11).Interpreter = 'latex';
KShm(11).YDisplayLabels = repmat(' ',5,4);
KShm(11).XLabel = '$\tau$';
clim([0 1])

nexttile
KShm(12) = heatmap(xvalues,yvalues,J0MeanKSErrP6,'Colormap',sky,'CellLabelColor','none','ColorbarVisible','off');
KShm(12).Interpreter = 'latex';
KShm(12).YDisplayLabels = repmat(' ',5,4);
KShm(12).XLabel = '$\tau$';
clim([0 1])

colorLims = vertcat(KShm.ColorLimits);
globalColorLim = [min(colorLims(:,1)), max(colorLims(:,2))];
set(KShm, 'ColorLimits', globalColorLim)
ax = axes(tcl,'visible','off','Colormap',KShm(1).Colormap,'CLim',globalColorLim);
cb = colorbar(ax);
cb.Layout.Tile = 'East';
set(cb,'TickLabelInterpreter', 'latex');

tcl.XLabel.Visible='on';
t1 = text(-0.35,0.15,'HMC--LF','interpreter','latex');
t1.Rotation = 90;
t2 = text(-0.35,-1.05,'CHMC','interpreter','latex');
t2.Rotation = 90;
set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
set(findall(gcf,'-property','FontSize'),'FontSize',14)

HeatMapPlotStr = strcat('compactHeatMapKS',datestr(now,'_dd-mm-yy_HH-MM-SS'));
fig.Position = [0,0,1000,260];
print('-dpng','-r400',HeatMapPlotStr)

%% Plots Heatmap error in W1 distance
fig = figure(2);
clf
tcl = tiledlayout(2,6,TileSpacing="tight");
Whm = gobjects(6,1); 

%% First row of heatmap (HMC)
nexttile
Whm(1) = heatmap(xvalues,yvalues,LFMeanWErrd400,'Colormap',sky,'CellLabelColor','none','ColorbarVisible','off');
Whm(1).Interpreter = 'latex';
Whm(1).XDisplayLabels = repmat(' ',4,5);
Whm(1).Title = '$d=400$';
Whm(1).YLabel = '$T$';
clim([0 1])

nexttile
Whm(2) = heatmap(LFMeanWErrd800,'Colormap',sky,'CellLabelColor','none','ColorbarVisible','off');
Whm(2).Interpreter = 'latex';
Whm(2).XDisplayLabels = repmat(' ',4,5);
Whm(2).YDisplayLabels = repmat(' ',5,4);
Whm(2).Title = '$d=800$';
clim([0 1])

nexttile
Whm(3) = heatmap(LFMeanWErrd1200,'Colormap',sky,'CellLabelColor','none','ColorbarVisible','off');
Whm(3).Interpreter = 'latex';
Whm(3).XDisplayLabels = repmat(' ',4,5);
Whm(3).YDisplayLabels = repmat(' ',5,4);
Whm(3).Title = '$d=1200$';
clim([0 1])

nexttile
Whm(4) = heatmap(xvalues,yvalues,LFMeanWErrP2,'Colormap',sky,'CellLabelColor','none','ColorbarVisible','off');
Whm(4).Interpreter = 'latex';
Whm(4).XDisplayLabels = repmat(' ',4,5);
Whm(4).YDisplayLabels = repmat(' ',5,4);
Whm(4).Title = '$p=2$';
clim([0 1])

nexttile
Whm(5) = heatmap(LFMeanWErrP4,'Colormap',sky,'CellLabelColor','none','ColorbarVisible','off');
Whm(5).Interpreter = 'latex';
Whm(5).XDisplayLabels = repmat(' ',4,5);
Whm(5).YDisplayLabels = repmat(' ',5,4);
Whm(5).Title = '$p=4$';
clim([0 1])

nexttile
Whm(6) = heatmap(LFMeanWErrP6,'Colormap',sky,'CellLabelColor','none','ColorbarVisible','off');
Whm(6).Interpreter = 'latex';
Whm(6).XDisplayLabels = repmat(' ',4,5);
Whm(6).YDisplayLabels = repmat(' ',5,4);
Whm(6).Title = '$p=6$';
clim([0 1])

%% Second row of heatmap (CHMC)
nexttile
Whm(7) = heatmap(xvalues,yvalues,J0MeanWErrd400,'Colormap',sky,'CellLabelColor','none','ColorbarVisible','off');
Whm(7).Interpreter = 'latex';
Whm(7).YLabel = '$T$';
Whm(7).XLabel = '$\tau$';
clim([0 1])

nexttile
Whm(8) = heatmap(xvalues,yvalues,J0MeanWErrd800,'Colormap',sky,'CellLabelColor','none','ColorbarVisible','off');
Whm(8).Interpreter = 'latex';
Whm(8).YDisplayLabels = repmat(' ',5,4);
Whm(8).XLabel = '$\tau$';
clim([0 1])

nexttile
Whm(9) = heatmap(xvalues,yvalues,J0MeanWErrd1200,'Colormap',sky,'CellLabelColor','none','ColorbarVisible','off');
Whm(9).Interpreter = 'latex';
Whm(9).YDisplayLabels = repmat(' ',5,4);
Whm(9).XLabel = '$\tau$';
clim([0 1])

nexttile
Whm(10) = heatmap(xvalues,yvalues,J0MeanWErrP2,'Colormap',sky,'CellLabelColor','none','ColorbarVisible','off');
Whm(10).Interpreter = 'latex';
Whm(10).YDisplayLabels = repmat(' ',5,4);
Whm(10).XLabel = '$\tau$';
clim([0 1])

nexttile
Whm(11) = heatmap(xvalues,yvalues,J0MeanWErrP4,'Colormap',sky,'CellLabelColor','none','ColorbarVisible','off');
Whm(11).Interpreter = 'latex';
Whm(11).YDisplayLabels = repmat(' ',5,4);
Whm(11).XLabel = '$\tau$';
clim([0 1])

nexttile
Whm(12) = heatmap(xvalues,yvalues,J0MeanWErrP6,'Colormap',sky,'CellLabelColor','none','ColorbarVisible','off');
Whm(12).Interpreter = 'latex';
Whm(12).YDisplayLabels = repmat(' ',5,4);
Whm(12).XLabel = '$\tau$';
clim([0 1])

colorLims = vertcat(Whm.ColorLimits);
globalColorLim = [min(colorLims(:,1)), max(colorLims(:,2))];
set(Whm, 'ColorLimits', globalColorLim)
ax = axes(tcl,'visible','off','Colormap',Whm(1).Colormap,'CLim',globalColorLim);
cb = colorbar(ax);
cb.Layout.Tile = 'East';
set(cb,'TickLabelInterpreter', 'latex');
anno = annotation('rectangle',[921/1000 16/260 12/1000 12/260], 'FaceColor','black');
anno.Parent.HandleVisibility = 'on';
anno = annotation('textbox',[935/1000 22/260 10/1000 10/260],'String','NaN','LineStyle','none','Interpreter','latex');
anno.Parent.HandleVisibility = 'on';

t1 = text(-0.35,0.15,'HMC--LF','interpreter','latex');
t1.Rotation = 90;
t2 = text(-0.35,-1.05,'CHMC','interpreter','latex');
t2.Rotation = 90;
set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
set(findall(gcf,'-property','FontSize'),'FontSize',14)

HeatMapPlotStr = strcat('compactHeatMapW1',datestr(now,'_dd-mm-yy_HH-MM-SS'));
fig.Position = [0,0,1000,260];
print('-dpng','-r400',HeatMapPlotStr)