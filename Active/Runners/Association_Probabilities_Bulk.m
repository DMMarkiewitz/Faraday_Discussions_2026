clearvars
clc

%% System properties and bulk quantities
params.fa  = 3;        % Anion functionality
params.fc  = 4;        % Cation functionality
params.xic = 0.4;      % Cation volume ratio
params.xia = 10.8;     % Anion volume ratio
params.MMw = 18.02e-3; % kg/mol, molar mass of water
params.L   = 0.2527;   % Association constant at 21 m

% Mole fraction grid used in the original script
mol_f_ref = 1 / ((1 / params.MMw / 0.5) + 1);
mol_f = [mol_f_ref:0.001:0.02-0.001, 0.02:0.0025:0.35];

% Convert to molality
molality = (1 / params.MMw) .* (1 ./ (1 ./ mol_f - 1));

% Preallocate outputs
n = numel(mol_f);
phi_p_b = zeros(n, 1);
phi_m_b = zeros(n, 1);
phi_0_b = zeros(n, 1);
pc   = zeros(n, 1);
pa   = zeros(n, 1);
pc0  = zeros(n, 1);
p0   = zeros(n, 1);

%% Compute bulk association probabilities
for i = 1:n
    [phi_p_b(i), phi_m_b(i), phi_0_b(i), pc(i), pa(i), pc0(i), p0(i)] = ...
        bulkProbabilities(mol_f(i), params);
end

%% Plot bulk association probabilities
fig = figure('Color', 'w');
ax = axes(fig);
hold(ax, 'on');

plot(ax, molality, pc,  'r-',  'LineWidth', 1.5, 'DisplayName', '$p_{+-}$');
plot(ax, molality, pc0, 'r--', 'LineWidth', 1.5, 'DisplayName', '$p_{+0}$');
plot(ax, molality, pa,  'b:',  'LineWidth', 1.5, 'DisplayName', '$p_{-+}$');
plot(ax, molality, p0,  'k-.', 'LineWidth', 1.5, 'DisplayName', '$p_{0+}$');

xlim(ax, [0, 21]);
ylim(ax, [-0.025, 1.175]);
xticks(ax, get(ax, 'XTick'));
yticks(ax, 0:0.2:1);
box(ax, 'on');

xlabel(ax, 'LiTFSI molality, m', 'Interpreter', 'latex', 'FontSize', 16);
ylabel(ax, 'Association probability, $p_{ij}$', 'Interpreter', 'latex', 'FontSize', 16);
pbaspect(ax, [(1 + sqrt(5))/2, 1, 1]);

set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultTextInterpreter', 'latex');

leg = legend(ax, 'Location', 'north', 'Orientation', 'horizontal', ...
    'Interpreter', 'latex', 'Color', 'none', 'Box', 'off');
leg.ItemTokenSize = [18, 18];

set(ax, 'XColor', 'k', 'YColor', 'k', 'FontSize', 14, 'LineWidth', 1.5, ...
    'Layer', 'top', 'Color', 'none');
fig.Renderer = 'painters';
set(fig, 'InvertHardcopy', 'on');

% Export
baseName = 'Sticky_WiSE_Bulk_21m_L_assoc_v1';
exportgraphics(fig, [baseName, '.eps'],  'Resolution', 600);
exportgraphics(fig, [baseName, '.jpeg'], 'Resolution', 600);
saveas(fig, [baseName, '.fig']);

%% Local functions
function [phi_p, phi_m, phi_0, pc, pa, pc0, p0] = bulkProbabilities(mf, params)
    % Volume fractions in bulk
    phi_p = params.xic / (1 / mf - 1 + params.xic + params.xia);
    phi_m = params.xia / params.xic * phi_p;
    phi_0 = 1 - phi_p - phi_m;

    % Reduced number densities used by the original expressions
    psi_p = params.fc * phi_p / params.xic;
    psi_m = params.fa * phi_m / params.xia;

    % Association probabilities
    pc  = p_c(phi_0, psi_p, psi_m, params.L);
    pa  = p_a(phi_0, psi_p, psi_m, params.L);
    pc0 = p_c0(phi_0, psi_p, psi_m, params.L);
    p0  = p_0(phi_0, psi_p, psi_m, params.L);
end

function value = p_c(phi_0, psi_p, psi_m, L)
    radicand = 4 * phi_0 * psi_p * (L - 1) + (L * (psi_m - psi_p) + psi_p + phi_0)^2;
    value = (phi_0 - psi_p + L * (psi_p + psi_m) - sqrt(radicand)) / (2 * (L - 1) * psi_p);
end

function value = p_a(phi_0, psi_p, psi_m, L)
    radicand = 4 * phi_0 * psi_p * (L - 1) + (L * (psi_m - psi_p) + psi_p + phi_0)^2;
    value = (phi_0 - psi_p + L * (psi_p + psi_m) - sqrt(radicand)) / (2 * (L - 1) * psi_m);
end

function value = p_c0(phi_0, psi_p, psi_m, L)
    radicand = 4 * phi_0 * psi_p * (L - 1) + (L * (psi_m - psi_p) + psi_p + phi_0)^2;
    value = (phi_0 + psi_p + L * (-psi_p + psi_m) - sqrt(radicand)) / (2 * (1 - L) * psi_p);
end

function value = p_0(phi_0, psi_p, psi_m, L)
    radicand = 4 * phi_0 * psi_p * (L - 1) + (L * (psi_m - psi_p) + psi_p + phi_0)^2;
    value = (phi_0 + psi_p + L * (-psi_p + psi_m) - sqrt(radicand)) / (2 * (1 - L) * phi_0);
end
