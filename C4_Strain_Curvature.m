clear; close all; clc; tic;

% Parameters
H = 50;
lambdaH = 8;

% α–ψ pairs (third mode for each α)
alphas = [0.5       1E-5   -0.5];
psis   = [sqrt(0.8474), sqrt(0.9637), sqrt(1.0480)];

z = 0:0.2:H;
z_H = z / H;
nZ = length(z);
nCases = length(alphas);

% Results
Displacement = NaN(nZ, nCases);
StrainRatio = NaN(nZ, nCases);
AvgStrain = NaN(nZ, nCases);
Curvature = NaN(nZ, nCases);
Stress = NaN(nZ, nCases);
for j = 1:nCases
    alpha = alphas(j);
    psi = psis(j);

    % Constants
    S1 = 0.5 * (1 + sqrt(1 - 4 * psi^2));
    S2 = 0.5 * (1 - sqrt(1 - 4 * psi^2));
    a1 = (1 - 2i*psi + sqrt(1 - 4*psi^2))/2;
    b1 = (1 + 2i*psi + sqrt(1 - 4*psi^2))/2;
    c1 = 1 + sqrt(1 - 4*psi^2);
    a2 = (1 - 2i*psi - sqrt(1 - 4*psi^2))/2;
    b2 = (1 + 2i*psi - sqrt(1 - 4*psi^2))/2;
    c2 = 1 - sqrt(1 - 4*psi^2);
    a3 = a1 + 1; b3 = b1 + 1; c3 = c1 + 1;
    a4 = a2 + 1; b4 = b2 + 1; c4 = c2 + 1;

    Z2 = -((exp(lambdaH) - alpha) / alpha);
    Z3 = -((1 - alpha) / alpha);

    % Surface
    F1SURF = hypergeom([a1, b1], c1, Z3);
    F2SURF = hypergeom([a2, b2], c2, Z3);
    PH1SURF = (-(lambdaH/H)) / (2 * alpha) * hypergeom([a3, b3], c3, Z3);
    PH2SURF = (-(lambdaH/H)) / (2 * alpha) * hypergeom([a4, b4], c4, Z3);
    E_surf = (lambdaH/H)*S1*(1 - alpha)^(-S2)*F1SURF + (1 - alpha)^S1*PH1SURF;
    M_surf = (lambdaH/H)*S2*(1 - alpha)^(-S1)*F2SURF + (1 - alpha)^S2*PH2SURF;
    uSurf = E_surf - ((1 - alpha)^sqrt(1 - 4*psi^2)) * (F1SURF/F2SURF) * M_surf;

    % Base
    F1Base = hypergeom([a1, b1], c1, Z2);
    F2Base = hypergeom([a2, b2], c2, Z2);
    PH1Base = (-(lambdaH/H) * exp(lambdaH)) / (2 * alpha) * hypergeom([a3, b3], c3, Z2);
    PH2Base = (-(lambdaH/H) * exp(lambdaH)) / (2 * alpha) * hypergeom([a4, b4], c4, Z2);
    E_base = (lambdaH/H) * exp(lambdaH) * S1 * (exp(lambdaH) - alpha)^(-S2) * F1Base + ...
             (exp(lambdaH) - alpha)^S1 * PH1Base;
    M_base = (lambdaH/H) * exp(lambdaH) * S2 * (exp(lambdaH) - alpha)^(-S1) * F2Base + ...
             (exp(lambdaH) - alpha)^S2 * PH2Base;
    uBase = E_base - ((1 - alpha)^sqrt(1 - 4*psi^2)) * (F1SURF/F2SURF) * M_base;

    % Loop over depth
    for i = 1:nZ
        zi = z(i);
        Z1 = -((exp(lambdaH * zi / H) - alpha) / alpha);
        F1 = hypergeom([a1, b1], c1, Z1);
        F2 = hypergeom([a2, b2], c2, Z1);
        PH1 = (-(lambdaH/H) * exp(lambdaH*zi/H)) / (2 * alpha) * hypergeom([a3, b3], c3, Z1);
        PH2 = (-(lambdaH/H) * exp(lambdaH*zi/H)) / (2 * alpha) * hypergeom([a4, b4], c4, Z1);

        E = (lambdaH/H) * exp(lambdaH*zi/H) * S1 * (exp(lambdaH*zi/H) - alpha)^(-S2) * F1 + ...
            (exp(lambdaH*zi/H) - alpha)^S1 * PH1;
        M = (lambdaH/H) * exp(lambdaH*zi/H) * S2 * (exp(lambdaH*zi/H) - alpha)^(-S1) * F2 + ...
            (exp(lambdaH*zi/H) - alpha)^S2 * PH2;
        u = E - ((1 - alpha)^sqrt(1 - 4*psi^2)) * (F1SURF/F2SURF) * M;

        % Displacement ratio
        Displacement(i, j) =u / uSurf;

        % Strain ratio
        numG = -psi^2 * (lambdaH/H)^2 * exp(2*lambdaH*zi/H) * ...
               ( (exp(lambdaH*zi/H) - alpha)^(S1-2) * F1 - ...
                 (exp(lambdaH*zi/H) - alpha)^(S2-2) * ...
                 (1 - alpha)^sqrt(1 - 4*psi^2) * (F1SURF * F2 / F2SURF) );
        denG = -psi^2 * (lambdaH/H)^2 * exp(2*lambdaH) * ...
               ( (exp(lambdaH) - alpha)^(S1-2) * F1Base - ...
                 (exp(lambdaH) - alpha)^(S2-2) * ...
                 (1 - alpha)^sqrt(1 - 4*psi^2) * (F1SURF * F2Base / F2SURF) );
        StrainRatio(i, j) = numG / denG;

        
         %Stress ratio
        SqrtTerm = (1-alpha)^sqrt(1 - 4*psi^2);
        ExpTerm_depth = exp(lambdaH * zi / H) - alpha;
        ExpTerm_base = exp(lambdaH) - alpha;

        Stress_depth = ExpTerm_depth^(S1) * F1 - SqrtTerm* (F1SURF/F2SURF) * ExpTerm_depth^(S2)* F2;

        Stress_Base = ExpTerm_base^(S1) * F1BASE - SqrtTerm* (F1SURF/F2SURF) *ExpTerm_base^(S2)*F2BASE;

        Stress(i, j) = Stress_depth / Stress_Base;

       % Curvature ratio
        Gamma_Surface = -psi^2 * (lambdaH/H)^2 * ...
               ( (1 - alpha)^(S1 - 2) * F1SURF ...
               - (1 - alpha)^(S2 - 2) * ...
                 (1 - alpha)^sqrt(1 - 4 * psi^2) * ...
                 (F1SURF * F2SURF / F2SURF) );  % simplifies

        numC = -(((psi^2 * (lambdaH/H)^2 * exp(2 * lambdaH * zi / H)) / ...
                 (exp(lambdaH * zi / H) - alpha)^2) * u + ...
                 ((2 * alpha * lambdaH / H) / (exp(lambdaH * zi / H) - alpha)) * numG);

        denC = -(((psi^2 * (lambdaH/H)^2) / ((1 - alpha)^2)) * uSurf + ...
                 ((2 * alpha * lambdaH / H) / (1 - alpha)) * Gamma_Surface);

        Curvature(i, j) = numC / denC;

    end
end

%  Plotting
figure;
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

titles  = {'u(z) / u(0)', '\gamma(z) / \gamma(H)', ...
           '\tau(z) / \tau(H)', '1/R(z) / 1/R(0)'};

letters = {'(a)','(b)','(c)','(d)'};
dataAll = {Displacement, StrainRatio, Stress, Curvature};
alphaLabels = {'\alpha = 0.5', '\alpha = 0.000001', '\alpha = -0.5'};

lineStyles = {'-','--',':'};

for k = 1:4
    ax = nexttile;
    hold(ax,'on');

    for j = 1:nCases
        plot(ax, dataAll{k}(:, j), z_H, ...
             'LineWidth',1.5, ...
             'LineStyle',lineStyles{j});
    end

    set(ax,'YDir','reverse', ...
           'XAxisLocation','top', ...
           'XLim',[-3 3], ...
           'YLim',[0 1], ...
           'YTick',0:0.2:1);


    xlabel(ax, titles{k});

    box(ax,'on');
end



toc;