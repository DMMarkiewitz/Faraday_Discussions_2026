clearvars
close all
clc

%%%%%%%%%%%%%%%%%%%%%
% Defining parameters
%%%%%%%%%%%%%%%%%%%%%
T = 300;                       % Temperature [K]
kB = 1.38e-23;                 % Boltzmann constant [J/K]
eCharge = 1.6e-19;             % Elementary charge [C]
beta = 1 / (kB * T);           % Inverse thermal energy [1/J]
fa = 3;                        % Anion functionality
fc = 4;                        % Cation functionality
xic = 0.4;                     % Cation volume ratio
xia = 10.8;                    % Anion volume ratio
repulsionCoeff = 0.4;          % Short-range repulsion coefficient
eps_r = 10.1;                  % Relative dielectric constant
eps0 = 8.85e-12;               % Vacuum permittivity [F/m]
epsilon = eps_r * eps0;
MMw = 18.02 / 1000;            % Water molar mass [kg/mol]

% Molality cases to run
xMolalValues = 12;

%%%%%%%%%%%%%%%%%%%%%%%%%
% Potential / field grids
%%%%%%%%%%%%%%%%%%%%%%%%%
phiMax = 0.1 * eCharge * beta;
dPhiStep = -0.01;
Phi = (phiMax:dPhiStep:-phiMax).';
nPhi = numel(Phi);
DPhi = (0:0.01:0.1).';
nDPhi = numel(DPhi);

% Start continuation from the index nearest Phi = 0 and move outward
[~, centerIdx] = min(abs(Phi));
phiTraversalOrder = [centerIdx:-1:1, centerIdx+1:nPhi];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solver / numerical parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solverTol = 1e-5;
optionsFS = optimoptions('fsolve', ...
    'Display', 'off', ...
    'StepTolerance', solverTol, ...
    'FunctionTolerance', solverTol, ...
    'MaxFunctionEvaluations', 1e4, ...
    'MaxIterations', 1e4);

% Diagnostic helper for the binding-law check
bindingLawResidual = @(x,y,z) z - x .* (1 - y) ./ (y .* (1 - x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main loop over molality
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m = xMolalValues

    [phi_p_b, v0, L] = get_composition_parameters(m, MMw, xia, xic);

    % Bulk quantities
    phi_m_b = xia * phi_p_b / xic;
    psi_p_b = fc * phi_p_b / xic;
    psi_m_b = fa * phi_m_b / xia;
    phi_0_b = 1 - phi_p_b - phi_m_b;

    pDipole = 1.67e-29 / sqrt(v0 * epsilon / (beta * (phi_p_b / xic + phi_m_b / xia)));

    % Bulk association probabilities
    pcBulk  = f_pc(phi_0_b, psi_p_b, psi_m_b, L);
    paBulk  = f_pa(phi_0_b, psi_p_b, psi_m_b, L);
    pc0Bulk = f_pc0(phi_0_b, psi_p_b, psi_m_b, L);
    p0Bulk  = f_p0(phi_0_b, psi_p_b, psi_m_b, L);

    % Initial guess in the bulk
    bulkGuess = [1, phi_p_b, phi_m_b, phi_0_b];

    % Allocate storage
    phiVals   = nan(nPhi, nDPhi, 3);
    tauVals   = nan(nPhi, nDPhi);
    etaVals   = nan(nPhi, nDPhi);   
    tauStore  = nan(nPhi, nDPhi);
    epsrVals  = nan(nPhi, nDPhi);
    paVals    = nan(nPhi, nDPhi);   % p_{-+}
    p0Vals    = nan(nPhi, nDPhi);   % p_{0+}
    rhoVals   = nan(nPhi, nDPhi);

    % Diagnostics
    exitflagMap            = nan(nPhi, nDPhi);
    residualNormMap        = nan(nPhi, nDPhi);
    isGelMask              = false(nPhi, nDPhi);
    isConstraintViolation  = false(nPhi, nDPhi);
    isConsistencyWarning   = false(nPhi, nDPhi);
    isSolveFailure         = false(nPhi, nDPhi);

    numSolveFailures       = 0;
    numConstraintViol      = 0;
    numGelPoints           = 0;
    numConsistencyWarnings = 0;

    lastSuccessfulGuess = bulkGuess;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Finding charge density map
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:nDPhi

        Lbar = L * sinch(pDipole * DPhi(j));

        for k = 1:numel(phiTraversalOrder)

            i = phiTraversalOrder(k);

            residualFun = @(x) [
                x(2) + x(3) + x(4) - 1;
                x(2) - phi_p_b * (1 - pcBulk)^fc * exp(-repulsionCoeff * Phi(i)) * x(1)^(xic + fc) / ...
                    f_pc0(x(4), fc * x(2) / xic, fa * x(3) / xia, Lbar)^fc;
                x(3) - phi_m_b * (1 - paBulk)^fa * exp(repulsionCoeff * Phi(i)) * x(1)^xia / ...
                    (1 - f_pa(x(4), fc * x(2) / xic, fa * x(3) / xia, Lbar))^fa;
                x(4) - sinch(pDipole * DPhi(j)) * phi_0_b * (1 - p0Bulk) * x(1) / ...
                    (1 - f_p0(x(4), fc * x(2) / xic, fa * x(3) / xia, Lbar))
            ];

            guessList = build_guess_list(i, j, k, phiTraversalOrder, tauVals, phiVals, lastSuccessfulGuess, bulkGuess);

            [keys, flag, resnorm] = solve_with_fallbacks(residualFun, guessList, bulkGuess, optionsFS);

            exitflagMap(i, j) = flag;
            residualNormMap(i, j) = resnorm;

            if flag <= 0 || any(~isfinite(keys)) || ~isreal(keys)
                isSolveFailure(i, j) = true;
                numSolveFailures = numSolveFailures + 1;
                continue
            end

            tauVal = keys(1);
            phi_c  = keys(2);
            phi_a  = keys(3);
            phi_0  = keys(4);

            [isPhysical, isConstraintBad] = validate_local_solution(tauVal, phi_c, phi_a, phi_0, fc, fa, xic, xia, solverTol);

            if ~isPhysical
                isSolveFailure(i, j) = true;
                numSolveFailures = numSolveFailures + 1;
                continue
            end

            if isConstraintBad
                isConstraintViolation(i, j) = true;
                numConstraintViol = numConstraintViol + 1;
                continue
            end

            % Local probabilities
            pc0Local = f_pc0(phi_0, fc * phi_c / xic, fa * phi_a / xia, Lbar); % p_{+0}
            paLocal  = f_pa(phi_0,  fc * phi_c / xic, fa * phi_a / xia, Lbar); % p_{-+}
            p0Local  = f_p0(phi_0,  fc * phi_c / xic, fa * phi_a / xia, Lbar); % p_{0+}
            pcPair   = 1 - pc0Local;                                           % p_{+-}

            if any(~isfinite([pc0Local, paLocal, p0Local])) || any(~isreal([pc0Local, paLocal, p0Local]))
                isSolveFailure(i, j) = true;
                numSolveFailures = numSolveFailures + 1;
                continue
            end

            % Store primary outputs
            phiVals(i, j, 1) = phi_c;
            phiVals(i, j, 2) = phi_a;
            phiVals(i, j, 3) = phi_0;

            tauVals(i, j) = tauVal;
            rhoVals(i, j) = phi_c / xic - phi_a / xia;

            etaVals(i, j)  = pc0Local;
            paVals(i, j)   = paLocal;
            p0Vals(i, j)   = p0Local;
            tauStore(i, j) = tauVal;
            epsrVals(i, j) = f_t_epsr(phi_p_b / xic, phi_m_b / xia, DPhi(j), pDipole, phi_0, p0Local);

            % Consistency checks
            check1 = abs(fa * phi_a * paLocal / xia - fc * phi_c * pcPair / xic);
            check2 = abs(fc * phi_c * pc0Local / xic - phi_0 * p0Local);
            check3 = abs(bindingLawResidual(paLocal, p0Local, Lbar));
            check4 = abs(phi_c + phi_a + phi_0 - 1);

            if any([check1, check2, check3, check4] > 10 * solverTol)
                isConsistencyWarning(i, j) = true;
                numConsistencyWarnings = numConsistencyWarnings + 1;
            end

            % Gelation check
            if (fa - 1) * (fc - 1) * paLocal * pcPair - 1 > 0
                isGelMask(i, j) = true;
                numGelPoints = numGelPoints + 1;
            end

            % Update continuation state
            lastSuccessfulGuess = [tauVal, phi_c, phi_a, phi_0];
        end
    end

    % Derived quantity
    pcVals = 1 - etaVals;

    % File tag
    caseTag = sprintf('WiSE_fa%d_fc%d_alpha_%0.2f_%dm', fa, fc, repulsionCoeff, m);

    % Save results and metadata
    save([caseTag '.mat'], ...
        'Phi', 'DPhi', 'rhoVals', 'phiVals', 'pcVals', 'etaVals', 'paVals', 'p0Vals', ...
        'tauStore', 'epsrVals', 'tauVals', 'exitflagMap', 'residualNormMap', ...
        'isGelMask', 'isConstraintViolation', 'isConsistencyWarning', 'isSolveFailure', ...
        'T', 'kB', 'eCharge', 'beta', 'fa', 'fc', 'xic', 'xia', 'repulsionCoeff', ...
        'eps_r', 'epsilon', 'm', 'v0', 'L', 'phi_p_b', 'phi_m_b', 'phi_0_b', ...
        'pcBulk', 'paBulk', 'pc0Bulk', 'p0Bulk', 'pDipole');

    % Summary
    fprintf('\nCase: %s\n', caseTag);
    fprintf('  Solve failures         : %d\n', numSolveFailures);
    fprintf('  Constraint violations  : %d\n', numConstraintViol);
    fprintf('  Gel points             : %d\n', numGelPoints);
    fprintf('  Consistency warnings   : %d\n', numConsistencyWarnings);

    % Plotting
    proximityToGel = 1 - (fc - 1) * (fa - 1) * (pcVals .* paVals);
    plot_results(Phi, DPhi, rhoVals, phiVals, etaVals, paVals, p0Vals, tauStore, epsrVals, proximityToGel, caseTag);
end

%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%

function [phi_p_b, v0, L] = get_composition_parameters(m, MMw, xia, xic)

    switch m
        case 12
            phi_p_b = 0.02526;
            v0 = 2.2863e-29;
            L = 0.2259;
        case 15
            phi_p_b = 0.02684;
            v0 = 2.2456e-29;
            L = 0.2312;
        case 21
            phi_p_b = 0.02888;
            v0 = 2.1744e-29;
            L = 0.2527;
        otherwise
            phi_p_b = 1 / (1 + (1 / (m * MMw) + xia) / xic);
            v0 = mean([2.1744e-29, 2.2863e-29, 2.2456e-29]);
            L = mean([0.2259, 0.2312, 0.2527]);
    end
end

function guessList = build_guess_list(i, j, k, phiTraversalOrder, tauVals, phiVals, lastSuccessfulGuess, bulkGuess)

    guessList = [];

    % Previous successful point along the continuation path at the same DPhi
    if k > 1
        prevI = phiTraversalOrder(k - 1);
        prevGuess = pack_guess(prevI, j, tauVals, phiVals);
        if ~isempty(prevGuess)
            guessList = [guessList; sanitize_initial_guess(prevGuess, bulkGuess)];
        end
    end

    % Same Phi from the previous DPhi slice
    if j > 1
        prevFieldGuess = pack_guess(i, j - 1, tauVals, phiVals);
        if ~isempty(prevFieldGuess)
            guessList = [guessList; sanitize_initial_guess(prevFieldGuess, bulkGuess)];
        end
    end

    % Previous path point from the previous DPhi slice
    if k > 1 && j > 1
        prevI = phiTraversalOrder(k - 1);
        crossGuess = pack_guess(prevI, j - 1, tauVals, phiVals);
        if ~isempty(crossGuess)
            guessList = [guessList; sanitize_initial_guess(crossGuess, bulkGuess)];
        end
    end

    % Last successful solution anywhere in the current run
    guessList = [guessList; sanitize_initial_guess(lastSuccessfulGuess, bulkGuess)];

    % Bulk fallback
    guessList = [guessList; sanitize_initial_guess(bulkGuess, bulkGuess)];

    % Remove duplicate rows
    guessList = unique(guessList, 'rows', 'stable');
end

function packed = pack_guess(i, j, tauVals, phiVals)

    tauVal = tauVals(i, j);
    phi_c  = phiVals(i, j, 1);
    phi_a  = phiVals(i, j, 2);
    phi_0  = phiVals(i, j, 3);

    if all(isfinite([tauVal, phi_c, phi_a, phi_0]))
        packed = [tauVal, phi_c, phi_a, phi_0];
    else
        packed = [];
    end
end

function guess = sanitize_initial_guess(guess, bulkGuess)

    if isempty(guess) || numel(guess) ~= 4 || any(~isfinite(guess))
        guess = bulkGuess;
        return
    end

    tauVal = real(guess(1));
    phis = real(guess(2:4));

    if ~isfinite(tauVal) || tauVal <= 0
        tauVal = max(bulkGuess(1), 1e-6);
    end

    badPhi = ~isfinite(phis);
    phis = bulkGuess(1, 2:4);
    phis = max(phis, 1e-10);

    phiSum = sum(phis);
    if phiSum <= 0
        phis = bulkGuess(1, 2:4);
    else
        phis = phis / phiSum;
    end

    guess = [tauVal, phis];
end

function [keys, bestFlag, bestResnorm] = solve_with_fallbacks(residualFun, guessList, bulkGuess, optionsFS)

    keys = bulkGuess;
    bestFlag = -1;
    bestResnorm = inf;

    for n = 1:size(guessList, 1)
        thisGuess = guessList(n, :);
        [trialKeys, ~, trialFlag] = fsolve(residualFun, thisGuess, optionsFS);

        if any(~isfinite(trialKeys)) || ~isreal(trialKeys)
            continue
        end

        trialResnorm = norm(residualFun(trialKeys), inf);

        if trialFlag > 0
            keys = trialKeys;
            bestFlag = trialFlag;
            bestResnorm = trialResnorm;
            return
        end

        if trialResnorm < bestResnorm
            keys = trialKeys;
            bestFlag = trialFlag;
            bestResnorm = trialResnorm;
        end
    end
end

function [isPhysical, isConstraintBad] = validate_local_solution(tauVal, phi_c, phi_a, phi_0, fc, fa, xic, xia, tolVal)

    values = [tauVal, phi_c, phi_a, phi_0];

    isPhysical = isreal(values) && all(isfinite(values)) && ...
        tauVal > 0 && ...
        phi_c >= -tolVal && phi_a >= -tolVal && phi_0 >= -tolVal && ...
        phi_c <= 1 + tolVal && phi_a <= 1 + tolVal && phi_0 <= 1 + tolVal;

    isConstraintBad = (fc * phi_c / xic - fa * phi_a / xia - phi_0) > tolVal;
end

function plot_results(Phi, DPhi, rhoVals, phiVals, etaVals, paVals, p0Vals, tauVals, epsrVals, proximityToGel, caseTag)

    [PhiGrid, DPhiGrid] = meshgrid(Phi, DPhi);

    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'DefaultTextInterpreter', 'latex');

    export_contour(PhiGrid, DPhiGrid, rhoVals',        '$\tilde{\rho}$',           [caseTag '_rho']);
    export_contour(PhiGrid, DPhiGrid, phiVals(:,:,1)', '$\bar{\phi}_+$',           [caseTag '_phi_c']);
    export_contour(PhiGrid, DPhiGrid, phiVals(:,:,2)', '$\bar{\phi}_-$',           [caseTag '_phi_a']);
    export_contour(PhiGrid, DPhiGrid, phiVals(:,:,3)', '$\bar{\phi}_0$',           [caseTag '_phi_0']);
    export_contour(PhiGrid, DPhiGrid, etaVals',        '$\bar{p}_{+0}$',           [caseTag '_p_p0']);
    export_contour(PhiGrid, DPhiGrid, (1 - etaVals)',  '$\bar{p}_{+-}$',           [caseTag '_p_pm']);
    export_contour(PhiGrid, DPhiGrid, paVals',         '$\bar{p}_{-+}$',           [caseTag '_p_mp']);
    export_contour(PhiGrid, DPhiGrid, p0Vals',         '$\bar{p}_{0+}$',           [caseTag '_p_0p']);
    export_contour(PhiGrid, DPhiGrid, tauVals',        '$\tau$',                   [caseTag '_tau']);
    export_contour(PhiGrid, DPhiGrid, epsrVals',       '$\tilde{\epsilon}$',       [caseTag '_epsr']);
    export_contour(PhiGrid, DPhiGrid, proximityToGel', 'Proximity to Gelation',    [caseTag '_prox_gel']);
end

function export_contour(PhiGrid, DPhiGrid, Z, colorbarLabel, baseName)

    fig = figure('Renderer', 'painters');
    contourf(PhiGrid, DPhiGrid, Z, 250, 'edgecolor', 'none');
    hold on

    xlabel('$u$', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('$\tilde{\nabla}u$', 'Interpreter', 'latex', 'FontSize', 16);

    cb = colorbar();
    ylabel(cb, colorbarLabel, 'Interpreter', 'latex', 'FontSize', 16);

    ax = gca;
    ax.FontSize = 14;
    ax.LineWidth = 1.5;
    ax.Layer = 'top';
    ax.XColor = 'k';
    ax.YColor = 'k';
    ax.Color = 'none';
    ax.Units = 'normalized';
    ax.Position = [0.1 0.15 0.675 0.8];

    set(fig, 'InvertHardcopy', 'on');

    exportgraphics(fig, [baseName '.eps'], 'Resolution', 600);
    exportgraphics(fig, [baseName '.jpeg'], 'Resolution', 600);
    savefig(fig, [baseName '.fig']);
end

function out = sinch(x)

    out = ones(size(x));

    mask = abs(x) > 1e-3;
    out(mask) = sinh(x(mask)) ./ x(mask);

    xs = x(~mask);
    out(~mask) = 1 + xs.^2 / 6 + xs.^4 / 120 + xs.^6 / 5040;
end

function pc = f_pc(phi_0, psi_p, psi_m, L)
    pc = (phi_0 - psi_p + L * (psi_p + psi_m) - ...
        sqrt(4 * phi_0 * psi_p * (L - 1) + (L * (psi_m - psi_p) + psi_p + phi_0).^2)) ...
        ./ (2 * (L - 1) * psi_p);
end

function pa = f_pa(phi_0, psi_p, psi_m, L)
    pa = (phi_0 - psi_p + L * (psi_p + psi_m) - ...
        sqrt(4 * phi_0 * psi_p * (L - 1) + (L * (psi_m - psi_p) + psi_p + phi_0).^2)) ...
        ./ (2 * (L - 1) * psi_m);
end

function pc0 = f_pc0(phi_0, psi_p, psi_m, L)
    pc0 = (phi_0 + psi_p + L * (-psi_p + psi_m) - ...
        sqrt(4 * phi_0 * psi_p * (L - 1) + (L * (psi_m - psi_p) + psi_p + phi_0).^2)) ...
        ./ (2 * (1 - L) * psi_p);
end

function p0 = f_p0(phi_0, psi_p, psi_m, L)
    p0 = (phi_0 + psi_p + L * (-psi_p + psi_m) - ...
        sqrt(4 * phi_0 * psi_p * (L - 1) + (L * (psi_m - psi_p) + psi_p + phi_0).^2)) ...
        ./ (2 * (1 - L) * phi_0);
end

function Lf = f_L_f(p, up)

    Lf = zeros(size(up));
    k = p .* up;

    mask = abs(up) > 1e-3;
    Lf(mask) = (coth(k(mask)) - 1 ./ k(mask)) ./ up(mask);

    ks = k(~mask);
    Lf(~mask) = p * (1/3 - ks.^2 / 45 + 2 * ks.^4 / 945 - ks.^6 / 4725);
end

function t_epsr = f_t_epsr(cpb, cmb, du, p, phi_0, p0p)
    t_epsr = cpb + cmb + phi_0 .* (1 - p0p) .* p .* f_L_f(p, du);
end