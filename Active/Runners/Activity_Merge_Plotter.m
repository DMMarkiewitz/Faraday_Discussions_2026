%Merging Activities
close all
clear all

%Loading data for merging
%21m
TFSI_21_p_2 = load("S_Sticky_WiSE_21m_01_a_03_p_02.mat");
TFSI_21_p_2 = TFSI_21_p_2.data;
TFSI_21_m_2 = load("S_Sticky_WiSE_21m_01_a_03_m_02.mat");
TFSI_21_m_2 = TFSI_21_m_2.data;

%15m
TFSI_15_p_2 = load("S_Sticky_WiSE_15m_01_a_03_p_02.mat");
TFSI_15_p_2 = TFSI_15_p_2.data;
TFSI_15_m_2 = load("S_Sticky_WiSE_15m_01_a_03_m_02.mat");
TFSI_15_m_2 = TFSI_15_m_2.data;

%12m
TFSI_12_p_2 = load("S_Sticky_WiSE_12m_01_a_03_p_02.mat");
TFSI_12_p_2 = TFSI_12_p_2.data;
TFSI_12_m_2 = load("S_Sticky_WiSE_12m_01_a_03_m_02.mat");
TFSI_12_m_2 = TFSI_12_m_2.data;

%0.5m bulk reference
TFSI_05_0 = load("S_Sticky_WiSE_05m_01_a_03_0_02.mat");
TFSI_05_0 = TFSI_05_0.data;

%For ploting
xrange = 2;

figure('Renderer', 'painters', 'units','inches','Position',[.01 .01 4.5 6])
ft = tiledlayout(2,1,'TileSpacing','none','padding','compact') %tight
letters = 18;
titlesize = 20;

%Plotting Lithium activity using bulk 0.5m as reference state
nexttile
hold on
plot(TFSI_21_p_2.spatial_nm,TFSI_21_p_2.mu_p-TFSI_05_0.mu_p(end),'-','Color','#940000','linewidth',1.5)
hold on
plot(TFSI_15_p_2.spatial_nm,TFSI_15_p_2.mu_p-TFSI_05_0.mu_p(end),'--','Color','#FF3333','linewidth',1.5)
plot(TFSI_12_p_2.spatial_nm,TFSI_12_p_2.mu_p-TFSI_05_0.mu_p(end),':','Color','#FF8670','linewidth',1.5)
xlim([0,xrange])
ylim([-0.5,16.45])
box on
leg=legend({'21m','15m','12m'},'location','north','interpreter','latex','Orientation','horizontal','box','off');   
leg.ItemTokenSize = [18,18];
ylabel(['Ln($\bar{a}_i$/$a_i^\theta$)'],'interpreter','latex', 'fontsize', 16)
text(.1,14.75,'\bf{a)}','interpreter','latex','FontSize',letters)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'DefaultTextInterpreter', 'latex')
ax = gca;
ax.FontSize = 16;
title('Li$^+$ Activity Shift', 'interpreter','latex', 'FontSize', titlesize)
ax.YTick = 0:5:15;
set(gca,'XTickLabel',[]);
g=gcf;
g.Renderer='painters';
set(gca, 'Xcolor', 'k');
set(gca, 'Ycolor', 'k');
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1.5);
set(gca, 'Layer', 'Top');
set(gca, 'color', 'none');
set(g,'InvertHardcopy','on');

nexttile
hold on
plot(TFSI_21_m_2.spatial_nm,TFSI_21_m_2.mu_p-TFSI_05_0.mu_p(end),'-','Color','#940000','linewidth',1.5)
hold on
plot(TFSI_15_m_2.spatial_nm,TFSI_15_m_2.mu_p-TFSI_05_0.mu_p(end),'--','Color','#FF3333','linewidth',1.5)
plot(TFSI_12_m_2.spatial_nm,TFSI_12_m_2.mu_p-TFSI_05_0.mu_p(end),':','Color','#FF8670','linewidth',1.5)
xlim([0,xrange])
ylim([6.5,14.75])
box on
xlabel(['Distance from electrode, nm'],'interpreter','latex', 'fontsize', 16)
ylabel(['Ln($\bar{a}_i$/$a_i^\theta$)'],'interpreter','latex', 'fontsize', 16)
text(.1,13.75,'\bf{b)}','interpreter','latex','FontSize',letters)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'DefaultTextInterpreter', 'latex')
ax = gca;
ax.FontSize = 16;
ax.YTick = 1:2:15;
g=gcf;
g.Renderer='painters';
set(gca, 'Xcolor', 'k');
set(gca, 'Ycolor', 'k');
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1.5);
set(gca, 'Layer', 'Top');
set(gca, 'color', 'none');
set(g,'InvertHardcopy','on');

%Saving Figure
exportgraphics(ft,['Activity_05m_Li_V2.eps'],'Resolution',600)
exportgraphics(ft,['Activity_05m_Li_V2.jpeg'],'Resolution',600)
saveas(ft,['Activity_05m_Li_V2.fig'])

figure('Renderer', 'painters', 'units','inches','Position',[.01 .01 4.5 6])
ft = tiledlayout(2,1,'TileSpacing','none','padding','compact')
letters = 18;
titlesize = 20;

%Plotting TFSI activity using bulk 0.5m as reference state
nexttile
hold on
plot(TFSI_21_p_2.spatial_nm,TFSI_21_p_2.mu_m-TFSI_05_0.mu_m(end),'-','Color','#000094','linewidth',1.5)
hold on
plot(TFSI_15_p_2.spatial_nm,TFSI_15_p_2.mu_m-TFSI_05_0.mu_m(end),'--','Color','#1957FF','linewidth',1.5)
plot(TFSI_12_p_2.spatial_nm,TFSI_12_p_2.mu_m-TFSI_05_0.mu_m(end),':','Color','#85B1FF','linewidth',1.5)
xlim([0,xrange])
ylim([-0.75,4.75])
box on
leg=legend({'21m','15m','12m'},'location','north','interpreter','latex','Orientation','horizontal','box','off');   
leg.ItemTokenSize = [18,18];
title('TFSI$^-$ Activity Shift', 'interpreter','latex', 'FontSize', titlesize)
text(.1,4.35,'\bf{e)}','interpreter','latex','FontSize',letters)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'DefaultTextInterpreter', 'latex')
ax = gca;
ax.FontSize = 16;
ax.YTick = -1:1:4;
set(gca,'XTickLabel',[]);
g=gcf;
g.Renderer='painters';
set(gca, 'Xcolor', 'k');
set(gca, 'Ycolor', 'k');
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1.5);
set(gca, 'Layer', 'Top');
set(gca, 'color', 'none');
set(g,'InvertHardcopy','on');

nexttile
hold on
plot(TFSI_21_m_2.spatial_nm,TFSI_21_m_2.mu_m-TFSI_05_0.mu_m(end),'-','Color','#000094','linewidth',1.5)
hold on
plot(TFSI_15_m_2.spatial_nm,TFSI_15_m_2.mu_m-TFSI_05_0.mu_m(end),'--','Color','#1957FF','linewidth',1.5)
plot(TFSI_12_m_2.spatial_nm,TFSI_12_m_2.mu_m-TFSI_05_0.mu_m(end),':','Color','#85B1FF','linewidth',1.5)
xlim([0,xrange])
ylim([-4.25,1.75])
box on
xlabel(['Distance from electrode, nm'],'interpreter','latex', 'fontsize', 16)
text(.1,1.25,'\bf{f)}','interpreter','latex','FontSize',letters)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'DefaultTextInterpreter', 'latex')
ax = gca;
ax.FontSize = 16;
ax.YTick = -4:1:1;
g=gcf;
g.Renderer='painters';
set(gca, 'Xcolor', 'k');
set(gca, 'Ycolor', 'k');
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1.5);
set(gca, 'Layer', 'Top');
set(gca, 'color', 'none');
set(g,'InvertHardcopy','on');

%Saving Figure
exportgraphics(ft,['Activity_05m_TFSI_V2.eps'],'Resolution',600)
exportgraphics(ft,['Activity_05m_TFSI_V2.jpeg'],'Resolution',600)
saveas(ft,['Activity_05m_TFSI_V2.fig'])

figure('Renderer', 'painters', 'units','inches','Position',[.01 .01 4.5 6])
ft = tiledlayout(2,1,'TileSpacing','none','padding','compact')
letters = 18;
titlesize = 20;

%Plotting Water activity using bulk 0.5m as reference state
nexttile
hold on
plot(TFSI_21_p_2.spatial_nm,TFSI_21_p_2.mu_0-TFSI_05_0.mu_0(end),'-','Color','#000000','linewidth',1.5)
hold on
plot(TFSI_15_p_2.spatial_nm,TFSI_15_p_2.mu_0-TFSI_05_0.mu_0(end),'--','Color','#636363','linewidth',1.5)
plot(TFSI_12_p_2.spatial_nm,TFSI_12_p_2.mu_0-TFSI_05_0.mu_0(end),':','Color','#858585','linewidth',1.5)
xlim([0,xrange])
ylim([-3.75,-0.5])
box on
leg=legend({'21m','15m','12m'},'location','north','interpreter','latex','Orientation','horizontal','box','off'); 
leg.ItemTokenSize = [18,18];
title('H$_2$O Activity Shift', 'interpreter','latex', 'FontSize', titlesize)
text(.1,-0.8,'\bf{c)}','interpreter','latex','FontSize',letters)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'DefaultTextInterpreter', 'latex')
ax = gca;
ax.FontSize = 16;
ax.YTick = -3:1:0;
set(gca,'XTickLabel',[]);
g=gcf;
g.Renderer='painters';
set(gca, 'Xcolor', 'k');
set(gca, 'Ycolor', 'k');
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1.5);
set(gca, 'Layer', 'Top');
set(gca, 'color', 'none');
set(g,'InvertHardcopy','on');

nexttile
hold on
plot(TFSI_21_m_2.spatial_nm,TFSI_21_m_2.mu_0-TFSI_05_0.mu_0(end),'-','Color','#000000','linewidth',1.5)
hold on
plot(TFSI_15_m_2.spatial_nm,TFSI_15_m_2.mu_0-TFSI_05_0.mu_0(end),'--','Color','#636363','linewidth',1.5)
plot(TFSI_12_m_2.spatial_nm,TFSI_12_m_2.mu_0-TFSI_05_0.mu_0(end),':','Color','#858585','linewidth',1.5)
xlim([0,xrange])
ylim([-3.75,-0.5])
box on
xlabel(['Distance from electrode, nm'],'interpreter','latex', 'fontsize', 16)
text(.1,-0.8,'\bf{d)}','interpreter','latex','FontSize',letters)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'DefaultTextInterpreter', 'latex')
ax = gca;
ax.FontSize = 16;
ax.YTick = -3:1:0;
g=gcf;
g.Renderer='painters';
set(gca, 'Xcolor', 'k');
set(gca, 'Ycolor', 'k');
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1.5);
set(gca, 'Layer', 'Top');
set(gca, 'color', 'none');
set(g,'InvertHardcopy','on');

%Saving Figure
exportgraphics(ft,['Activity_05m_Water_V2.eps'],'Resolution',600)
exportgraphics(ft,['Activity_05m_Water_V2.jpeg'],'Resolution',600)
saveas(ft,['Activity_05m_Water_V2.fig'])