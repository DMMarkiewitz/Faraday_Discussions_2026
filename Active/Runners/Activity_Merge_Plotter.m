close all
clearvars

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultTextInterpreter', 'latex');

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end

xRange = 2;
labelFontSize = 18;
titleFontSize = 20;
lineWidth = 1.5;
outputDpi = 600;

% Load data once.
data.p21 = loadActivityData(fullfile(scriptDir, 'S_Sticky_WiSE_21m_01_a_03_p_02.mat'));
data.m21 = loadActivityData(fullfile(scriptDir, 'S_Sticky_WiSE_21m_01_a_03_m_02.mat'));
data.p15 = loadActivityData(fullfile(scriptDir, 'S_Sticky_WiSE_15m_01_a_03_p_02.mat'));
data.m15 = loadActivityData(fullfile(scriptDir, 'S_Sticky_WiSE_15m_01_a_03_m_02.mat'));
data.p12 = loadActivityData(fullfile(scriptDir, 'S_Sticky_WiSE_12m_01_a_03_p_02.mat'));
data.m12 = loadActivityData(fullfile(scriptDir, 'S_Sticky_WiSE_12m_01_a_03_m_02.mat'));
ref = loadActivityData(fullfile(scriptDir, 'S_Sticky_WiSE_05m_01_a_03_0_02.mat'));
refMuP = ref.mu_p(end);
refMuM = ref.mu_m(end);
refMu0 = ref.mu_0(end);

series = [
    struct('name', '21m', 'style', '-',  'color', '#940000', 'short', 'p21', 'long', 'm21'), ...
    struct('name', '15m', 'style', '--', 'color', '#FF3333', 'short', 'p15', 'long', 'm15'), ...
    struct('name', '12m', 'style', ':',  'color', '#FF8670', 'short', 'p12', 'long', 'm12')
];

plotActivityFigure( ...
    'Activity_05m_Li_V2', ...
    'Li$^+$ Activity Shift', ...
    'Ln($\bar{a}_i$/$a_i^\theta$)', ...
    'mu_p', refMuP, ...
    [ -0.5, 16.45 ], [6.5, 14.75], ...
    [0 5 15], [1 2 15], ...
    'a)', 'b)', ...
    '#000000', '#000000', '#000000', ...
    series, data, xRange, labelFontSize, titleFontSize, lineWidth, outputDpi);

plotActivityFigure( ...
    'Activity_05m_TFSI_V2', ...
    'TFSI$^-$ Activity Shift', ...
    'Ln($\bar{a}_i$/$a_i^\theta$)', ...
    'mu_m', refMuM, ...
    [ -0.75, 4.75 ], [ -4.25, 1.75 ], ...
    [-1 1 4], [-4 -3 -2 -1 0 1], ...
    'e)', 'f)', ...
    '#000094', '#1957FF', '#85B1FF', ...
    series, data, xRange, labelFontSize, titleFontSize, lineWidth, outputDpi);

plotActivityFigure( ...
    'Activity_05m_Water_V2', ...
    'H$_2$O Activity Shift', ...
    'Ln($\bar{a}_i$/$a_i^\theta$)', ...
    'mu_0', refMu0, ...
    [ -3.75, -0.5 ], [ -3.75, -0.5 ], ...
    [-3 -2 -1 0], [-3 -2 -1 0], ...
    'c)', 'd)', ...
    '#000000', '#636363', '#858585', ...
    series, data, xRange, labelFontSize, titleFontSize, lineWidth, outputDpi);

function plotActivityFigure(fileStem, figureTitle, yLabelText, valueField, refValue, yLimTop, yLimBottom, yTicksTop, yTicksBottom, topLabel, bottomLabel, topColor, midColor, bottomColor, series, data, xRange, labelFontSize, titleFontSize, lineWidth, outputDpi)
    fig = figure('Renderer', 'painters', 'Units', 'inches', 'Position', [0.01 0.01 4.5 6]);
    tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'none', 'Padding', 'compact');

    plotPanel(tl, data, series, valueField, refValue, xRange, yLimTop, yTicksTop, ...
        topLabel, false, true, figureTitle, yLabelText, topColor, midColor, bottomColor, ...
        labelFontSize, titleFontSize, lineWidth);

    plotPanel(tl, data, series, valueField, refValue, xRange, yLimBottom, yTicksBottom, ...
        bottomLabel, true, false, '', '', topColor, midColor, bottomColor, ...
        labelFontSize, titleFontSize, lineWidth);

    exportgraphics(fig, [fileStem '.eps'], 'Resolution', outputDpi);
    exportgraphics(fig, [fileStem '.jpeg'], 'Resolution', outputDpi);
    saveas(fig, [fileStem '.fig']);
end

function plotPanel(tl, data, series, valueField, refValue, xRange, yLim, yTicks, panelLabel, showXLabel, addLegend, figureTitle, yLabelText, topColor, midColor, bottomColor, labelFontSize, titleFontSize, lineWidth)
    ax = nexttile(tl);
    hold(ax, 'on');

    colors = {topColor, midColor, bottomColor};
    for k = 1:numel(series)
        s = series(k);
        x = data.(s.short).spatial_nm;
        y = data.(s.short).(valueField) - refValue;
        plot(ax, x, y, s.style, 'Color', colors{k}, 'LineWidth', lineWidth);
    end

    xlim(ax, [0, xRange]);
    ylim(ax, yLim);
    box(ax, 'on');
    set(ax, 'XColor', 'k', 'YColor', 'k', 'FontSize', 14, 'LineWidth', 1.5, 'Layer', 'Top', 'Color', 'none');
    ax.YTick = yTicks;

    if showXLabel
        xlabel(ax, 'Distance from electrode, nm', 'Interpreter', 'latex', 'FontSize', 16);
    else
        ax.XTickLabel = [];
    end

    if ~isempty(yLabelText)
        ylabel(ax, yLabelText, 'Interpreter', 'latex', 'FontSize', 16);
    end

    if ~isempty(figureTitle)
        title(ax, figureTitle, 'Interpreter', 'latex', 'FontSize', titleFontSize);
    end

    text(ax, 0.1, yLim(2) - 0.12 * range(yLim), ['\bf{' panelLabel '}'], 'Interpreter', 'latex', 'FontSize', labelFontSize);

    if addLegend
        leg = legend(ax, {series.name}, 'Location', 'north', 'Interpreter', 'latex', 'Orientation', 'horizontal', 'Box', 'off');
        leg.ItemTokenSize = [18, 18];
    end
end

function data = loadActivityData(matFile)
    if ~isfile(matFile)
        error('Activity_Merge_Plotter_clean:MissingFile', 'Could not find required file: %s', matFile);
    end

    loaded = load(matFile);
    if ~isfield(loaded, 'data')
        error('Activity_Merge_Plotter_clean:MissingVariable', 'File %s does not contain a variable named "data".', matFile);
    end
    data = loaded.data;
end
