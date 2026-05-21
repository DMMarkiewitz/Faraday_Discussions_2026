clearvars
clc

%% System properties and bulk quantities
fa  = 3;          % Anion functionality
fc  = 4;          % Cation functionality
xic = 0.4;        % Cation volume ratio
xia = 10.8;       % Anion volume ratio
MMw = 18.02/1000;  % Water molar mass [kg/mol]
L   = 0.2527;     % Association constant at 21 m

% Composition grid
mol_f_ref = 1 / ((1 / MMw / 0.5) + 1);
mol_f = [mol_f_ref:0.001:0.02-0.001, 0.02:0.0025:0.35];
molality = (1 / MMw) .* (1 ./ (1 ./ mol_f - 1));

nPts = numel(mol_f);
phi_p_b = zeros(nPts, 1);
phi_m_b = zeros(nPts, 1);
phi_0_b = zeros(nPts, 1);
pc   = zeros(nPts, 1);
pa   = zeros(nPts, 1);
pc0  = zeros(nPts, 1);
p0   = zeros(nPts, 1);

%% Bulk composition and association probabilities
for k = 1:nPts
    [phi_p_b(k), phi_m_b(k), phi_0_b(k)] = bulkFractions(mol_f(k), xic, xia);

    psi_p = fc * phi_p_b(k) / xic;
    psi_m = fa * phi_m_b(k) / xia;

    pc(k)  = f_pc(phi_0_b(k), psi_p, psi_m, L);
    pa(k)  = f_pa(phi_0_b(k), psi_p, psi_m, L);
    pc0(k) = f_pc0(phi_0_b(k), psi_p, psi_m, L);
    p0(k)  = f_p0(phi_0_b(k), psi_p, psi_m, L);
end

%% Activity predictions
mu_p = 1 + fc + log((1 + fc / xic) .* phi_p_b .* (pc0 .^ fc));
mu_m = 1 + log(phi_m_b .* (1 - pa) .^ fa);
mu_0 = 1 + log(phi_0_b .* (1 - p0));

% Reference-state shift used in the original script
mu_p_hydrated = mu_p - mu_p(1) - (4 * (mu_0 - mu_0(1)));
mu_p_unhydrated = mu_p - mu_p(1);
mu_m_rel = mu_m - mu_m(1);
mu_0_rel = mu_0 - mu_0(1);

%% Plot
fig = figure('Color', 'w');
ax = axes(fig);
hold(ax, 'on');

plot(ax, molality, mu_p_hydrated, 'r-',  'LineWidth', 1.5);
plot(ax, molality, mu_p_unhydrated, 'r--', 'LineWidth', 1.5);
plot(ax, molality, mu_m_rel, 'b:', 'LineWidth', 1.5);
plot(ax, molality, mu_0_rel, 'k-.', 'LineWidth', 1.5);

xlim(ax, [0, 21]);
ylim(ax, [-5.5, 16]);
box(ax, 'on');

leg = legend(ax, {'Li$^+$', 'Hydrated Li$^+$', 'TFSI$^-$', 'H$_2$O'}, ...
    'Location', 'north', 'Interpreter', 'latex', 'Orientation', 'horizontal', 'Box', 'off');
leg.ItemTokenSize = [18, 18];

pbaspect(ax, [(1 + sqrt(5)) / 2, 1, 1]);

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultTextInterpreter', 'latex');

xlabel(ax, 'LiTFSI molality, m', 'Interpreter', 'latex', 'FontSize', 16);
ylabel(ax, 'Ln($a_i$/$a_i^\theta$)', 'Interpreter', 'latex', 'FontSize', 16);

set(ax, ...
    'XColor', 'k', ...
    'YColor', 'k', ...
    'FontSize', 14, ...
    'LineWidth', 1.5, ...
    'Layer', 'top', ...
    'Color', 'none');

set(fig, 'Renderer', 'painters', 'InvertHardcopy', 'on');

% Export outputs
baseName = 'Dual_Sticky_WiSE_Bulk_21m_L_05m_ref';
exportgraphics(fig, [baseName '.eps'], 'Resolution', 600);
exportgraphics(fig, [baseName '.jpeg'], 'Resolution', 600);
savefig(fig, [baseName '.fig']);

%% Local functions
function [phi_p_b, phi_m_b, phi_0_b] = bulkFractions(mol_f, xic, xia)
    phi_p_b = xic / (1 / mol_f - 1 + xic + xia);
    phi_m_b = (xia / xic) * phi_p_b;
    phi_0_b = 1 - phi_p_b - phi_m_b;
end

function pc = f_pc(phi_0, psi_p, psi_m, L)
    discriminant = 4 * phi_0 * psi_p * (L - 1) + ...
        (L * (psi_m - psi_p) + psi_p + phi_0)^2;
    pc = (phi_0 - psi_p + L * (psi_p + psi_m) - sqrt(discriminant)) / ...
        (2 * (L - 1) * psi_p);
end

function pa = f_pa(phi_0, psi_p, psi_m, L)
    discriminant = 4 * phi_0 * psi_p * (L - 1) + ...
        (L * (psi_m - psi_p) + psi_p + phi_0)^2;
    pa = (phi_0 - psi_p + L * (psi_p + psi_m) - sqrt(discriminant)) / ...
        (2 * (L - 1) * psi_m);
end

function pc0 = f_pc0(phi_0, psi_p, psi_m, L)
    discriminant = 4 * phi_0 * psi_p * (L - 1) + ...
        (L * (psi_m - psi_p) + psi_p + phi_0)^2;
    pc0 = (phi_0 + psi_p + L * (-psi_p + psi_m) - sqrt(discriminant)) / ...
        (2 * (1 - L) * psi_p);
end

function p0 = f_p0(phi_0, psi_p, psi_m, L)
    discriminant = 4 * phi_0 * psi_p * (L - 1) + ...
        (L * (psi_m - psi_p) + psi_p + phi_0)^2;
    p0 = (phi_0 + psi_p + L * (-psi_p + psi_m) - sqrt(discriminant)) / ...
        (2 * (1 - L) * phi_0);
end
