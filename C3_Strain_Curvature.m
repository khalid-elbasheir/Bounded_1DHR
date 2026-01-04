clear all; close all; clc; tic;

% Common Parameters
zeta = 0.05;
H = 50;
Vinf = 200;
Vinf_complex = Vinf * sqrt(1 + (2 * 1i * zeta));
lambdaH = 8;
alpha = 0.5;
psi = [sqrt(0.0378) sqrt(0.3172) sqrt(0.8474)];

z = 0:0.2:H;
z_H = z / H;
nModes = length(psi);
nPoints = length(z);

% Initialize result holders
Displacement = NaN(nPoints, nModes);
StrainRatio = NaN(nPoints, nModes);
AvgStrain = NaN(nPoints, nModes);
Curvature = NaN(nPoints, nModes);

% Begin mode loop
for j = 1:nModes
    % Precompute constants per mode
    S1 = 0.5*(1 + sqrt(1 - 4 * psi(j)^2));
    S2 = 0.5*(1 - sqrt(1 - 4 * psi(j)^2));
    a1 = (1 - 2i*psi(j) + sqrt(1 - 4*psi(j)^2))/2;
    b1 = (1 + 2i*psi(j) + sqrt(1 - 4*psi(j)^2))/2;
    c1 = 1 + sqrt(1 - 4*psi(j)^2);
    a2 = (1 - 2i*psi(j) - sqrt(1 - 4*psi(j)^2))/2;
    b2 = (1 + 2i*psi(j) - sqrt(1 - 4*psi(j)^2))/2;
    c2 = 1 - sqrt(1 - 4*psi(j)^2);
    a3 = a1 + 1; b3 = b1 + 1; c3 = c1 + 1;
    a4 = a2 + 1; b4 = b2 + 1; c4 = c2 + 1;

    % Surface and base evaluations (constants)
    Z2 = -((exp(lambdaH) - alpha) / alpha);
    Z3 = -((1 - alpha) / alpha);
    F1SURF = hypergeom([a1, b1], c1, Z3);
    F2SURF = hypergeom([a2, b2], c2, Z3);
    PH1SURF = (-(lambdaH/H)) / (2 * alpha) * hypergeom([a3, b3], c3, Z3);
    PH2SURF = (-(lambdaH/H)) / (2 * alpha) * hypergeom([a4, b4], c4, Z3);
    E_surf = (lambdaH/H) * S1 * (1 - alpha)^(-S2) * F1SURF + (1 - alpha)^S1 * PH1SURF;
    M_surf = (lambdaH/H) * S2 * (1 - alpha)^(-S1) * F2SURF + (1 - alpha)^S2 * PH2SURF;
    uSurf = E_surf - ((1 - alpha)^sqrt(1 - 4*psi(j)^2)) * (F1SURF/F2SURF) * M_surf;
    Den_u = uSurf;

    for i = 1:nPoints
        zi = z(i);
        Z1 = -((exp(lambdaH * zi / H) - alpha) / alpha);

        % Hypergeom evaluations
        F1Depth = hypergeom([a1, b1], c1, Z1);
        PH1Depth = (-(lambdaH/H) * exp(lambdaH * zi / H)) / (2 * alpha) * hypergeom([a3, b3], c3, Z1);
        F2Depth = hypergeom([a2, b2], c2, Z1);
        PH2Depth = (-(lambdaH/H) * exp(lambdaH * zi / H)) / (2 * alpha) * hypergeom([a4, b4], c4, Z1);

        % E and M
        E_depth = (lambdaH/H) * exp(lambdaH * zi / H) * S1 * (exp(lambdaH * zi / H) - alpha)^(-S2) * F1Depth ...
                  + (exp(lambdaH * zi / H) - alpha)^S1 * PH1Depth;
        M_depth = (lambdaH/H) * exp(lambdaH * zi / H) * S2 * (exp(lambdaH * zi / H) - alpha)^(-S1) * F2Depth ...
                  + (exp(lambdaH * zi / H) - alpha)^S2 * PH2Depth;

        % Displacement
        Numer_u = E_depth - ((1 - alpha)^sqrt(1 - 4*psi(j)^2)) * (F1SURF/F2SURF) * M_depth;
        Displacement(i, j) = Numer_u / Den_u;

        % Gamma (strain)
        NumStrain = -psi(j)^2 * (lambdaH/H)^2 * exp(2 * lambdaH * zi / H) * ...
                   ( (exp(lambdaH * zi / H) - alpha)^(S1 - 2) * F1Depth ...
                   - (exp(lambdaH * zi / H) - alpha)^(S2 - 2) * (1 - alpha)^sqrt(1 - 4*psi(j)^2) ...
                   * (F1SURF * F2Depth / F2SURF) );
        BaseStrain = -psi(j)^2 * (lambdaH/H)^2 * exp(2 * lambdaH) * ...
                   ( (exp(lambdaH) - alpha)^(S1 - 2) * hypergeom([a1, b1], c1, Z2) ...
                   - (exp(lambdaH) - alpha)^(S2 - 2) * (1 - alpha)^sqrt(1 - 4*psi(j)^2) ...
                   * (F1SURF * hypergeom([a2, b2], c2, Z2) / F2SURF) );
        StrainRatio(i, j) = NumStrain / BaseStrain;

        % Average Strain
        E_base = (lambdaH/H) * exp(lambdaH) * S1 * (exp(lambdaH) - alpha)^(-S2) * hypergeom([a1, b1], c1, Z2) ...
                 + (exp(lambdaH) - alpha)^S1 * ((-(lambdaH/H) * exp(lambdaH)) / (2 * alpha)) * hypergeom([a3, b3], c3, Z2);
        M_base = (lambdaH/H) * exp(lambdaH) * S2 * (exp(lambdaH) - alpha)^(-S1) * hypergeom([a2, b2], c2, Z2) ...
                 + (exp(lambdaH) - alpha)^S2 * ((-(lambdaH/H) * exp(lambdaH)) / (2 * alpha)) * hypergeom([a4, b4], c4, Z2);
        uBase = E_base - ((1 - alpha)^sqrt(1 - 4*psi(j)^2)) * (F1SURF/F2SURF) * M_base;
        Den_avg = (uSurf - uBase) / H;
        AvgStrain(i, j) = NumStrain / Den_avg;

        % Curvature
        Gamma_Depth = NumStrain;
        uDepth = Numer_u;
        NumCurv = -(((psi(j)^2 * (lambdaH/H)^2 * exp(2 * lambdaH * zi / H)) / ...
                  (exp(lambdaH * zi / H) - alpha)^2) * uDepth + ...
                  ((2 * alpha * lambdaH / H) / (exp(lambdaH * zi / H) - alpha)) * Gamma_Depth);
        Gamma_Surface = -psi(j)^2 * (lambdaH/H)^2 * 1 * ...
                   ( (1 - alpha)^(S1 - 2) * F1SURF ...
                   - (1 - alpha)^(S2 - 2) * (1 - alpha)^sqrt(1 - 4*psi(j)^2) ...
                   * (F1SURF * F2SURF / F2SURF) );
        DenCurv = -(  (((psi(j)^(2) *(lambdaH/H)^(2))/((1 - alpha)^(2)) )*uSurf ) + ( ((2*alpha*lambdaH/H)/((1 - alpha)))*Gamma_Surface)  );

        Curvature(i, j) = NumCurv / DenCurv;
    end
end

%% Plotting
figure;
titles = {'u(z) / u(0)', '\gamma(z) / \gamma(H)', ...
          '\gamma(z)H / [u(0)-u(H)]', '1/R(z) / 1/R(0)'};
dataAll = {Displacement, StrainRatio, AvgStrain, Curvature};
for k = 1:4
    subplot(2,2,k);
    hold on;
    for j = 1:nModes
        plot(dataAll{k}(:, j), z_H, 'LineWidth', 1.5);
    end
    set(gca, 'YDir', 'reverse', 'XAxisLocation', 'top');
    xlabel(titles{k});
    ylabel('z/H');
    title(titles{k});
    legend('\psi_{1}','\psi_{2}','\psi_{3}','Location','best');
    grid on;
end

toc;
