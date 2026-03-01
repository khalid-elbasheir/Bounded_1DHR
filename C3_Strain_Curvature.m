clear ; 
close all; 
clc; 
tic;

% Common Parameters
H = 50;
lambdaH = 8;
alpha = 0.5;
psi = [sqrt(0.0378) sqrt(0.3172) sqrt(0.8474)];

z = 0:0.2:H;
z_H = z / H;
nModes = length(psi);
nPoints = length(z);

% Initialize
Displacement = NaN(nPoints, nModes);
StrainRatio = NaN(nPoints, nModes);
AvgStrain = NaN(nPoints, nModes);
Curvature = NaN(nPoints, nModes);
Stress = NaN(nPoints, nModes);

% loop
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

    % Surface and base functions
    Z2 = -((exp(lambdaH) - alpha) / alpha);
    Z3 = -((1 - alpha) / alpha);
    F1SURF = hypergeom([a1, b1], c1, Z3);
    F2SURF = hypergeom([a2, b2], c2, Z3);
    F1BASE = hypergeom([a1, b1], c1, Z2);
    F2BASE = hypergeom([a2, b2], c2, Z2);
    PH1SURF = (-(lambdaH/H)) / (2 * alpha) * hypergeom([a3, b3], c3, Z3);
    PH2SURF = (-(lambdaH/H)) / (2 * alpha) * hypergeom([a4, b4], c4, Z3);
    E_surf = (lambdaH/H) * S1 * (1 - alpha)^(-S2) * F1SURF + (1 - alpha)^S1 * PH1SURF;
    M_surf = (lambdaH/H) * S2 * (1 - alpha)^(-S1) * F2SURF + (1 - alpha)^S2 * PH2SURF;
    uSurf = E_surf - ((1 - alpha)^sqrt(1 - 4*psi(j)^2)) * (F1SURF/F2SURF) * M_surf;
    Den_u = uSurf;

    for i = 1:nPoints
        zi = z(i);
        Z1 = -((exp(lambdaH * zi / H) - alpha) / alpha);

        % Hypergeom with depth
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

        % Gamma
        NumStrain = -psi(j)^2 * (lambdaH/H)^2 * exp(2 * lambdaH * zi / H) * ...
                   ( (exp(lambdaH * zi / H) - alpha)^(S1 - 2) * F1Depth ...
                   - (exp(lambdaH * zi / H) - alpha)^(S2 - 2) * (1 - alpha)^sqrt(1 - 4*psi(j)^2) ...
                   * (F1SURF * F2Depth / F2SURF) );
        BaseStrain = -psi(j)^2 * (lambdaH/H)^2 * exp(2 * lambdaH) * ...
                   ( (exp(lambdaH) - alpha)^(S1 - 2) * hypergeom([a1, b1], c1, Z2) ...
                   - (exp(lambdaH) - alpha)^(S2 - 2) * (1 - alpha)^sqrt(1 - 4*psi(j)^2) ...
                   * (F1SURF * hypergeom([a2, b2], c2, Z2) / F2SURF) );
        StrainRatio(i, j) = NumStrain / BaseStrain;
        
        %Stress
        SqrtTerm = (1-alpha)^sqrt(1 - 4*psi(j)^2);
        ExpTerm_depth = exp(lambdaH * zi / H) - alpha;
        ExpTerm_base = exp(lambdaH) - alpha;

        Stress_depth = ExpTerm_depth^(S1) * F1Depth - SqrtTerm* (F1SURF/F2SURF) * ExpTerm_depth^(S2)* F2Depth;

        Stress_Base = ExpTerm_base^(S1) * F1BASE - SqrtTerm* (F1SURF/F2SURF) *ExpTerm_base^(S2)*F2BASE;

        Stress(i, j) = Stress_depth / Stress_Base;

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

% Plotting
figure;
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

titles  = {'u(z) / u(0)', '\gamma(z) / \gamma(H)', ...
           '\tau(z) / \tau(H)', '1/R(z) / 1/R(0)'};
letters = {'(a)','(b)','(c)','(d)'};
dataAll = {Displacement, StrainRatio, Stress, Curvature};

lineStyles = {'-','--',':'};

for k = 1:4
    ax = nexttile;
    hold(ax,'on');


    for j = 1:nModes
        plot(ax, dataAll{k}(:,j), z_H, ...
             'LineWidth',1.5, 'LineStyle',lineStyles{j});
    end

    % axes
    set(ax,'YDir','reverse', ...          
           'XAxisLocation','top', ...
           'XLim',[-3 3], ...
           'YLim',[0 1], ...
           'YTick',0:0.2:1);
    xlabel(ax, titles{k});

    box(ax,'on');
end
toc;
