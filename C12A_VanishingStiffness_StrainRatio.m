clear
close all
clc
tic

% Alpha is Fixed at 1 - 4 types of LambdaH - 3 Frequencies
H     = 50;


lambdaH_list = [5, 8, 10, 15];
psi_values_list = {
    [sqrt(0.0829), sqrt(0.3367), sqrt(0.5646)];
    [sqrt(0.0366), sqrt(0.2402), sqrt(0.9333)];
    [sqrt(0.024),  sqrt(0.1806), sqrt(0.9333)];
    [sqrt(0.0109), sqrt(0.2271), sqrt(0.9151)];
};

z   = 1E-12:0.02:H;
z_H = z ./ H;
nZ  = length(z);

figure('Color','w');

tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

letters    = {'(a)','(b)','(c)','(d)'};      
lineStyles = {'-','--',':'};                
% ------------------------------------------------------------

for lambdaIdx = 1:length(lambdaH_list)
    lambdaH = lambdaH_list(lambdaIdx);
    psi     = psi_values_list{lambdaIdx};
    nModes  = length(psi);

    STRAIN_MATRIX = NaN(nZ, nModes);

    % ====== Compute strain ratio profiles ======
    for j = 1:nModes
        psi_j = psi(j);

        S1 = 0.5 * (1 + sqrt(1 - 4 * psi_j^2));
        a1 = (1 - 2*1i*psi_j + sqrt(1 - 4*psi_j^2)) / 2;
        b1 = (1 + 2*1i*psi_j + sqrt(1 - 4*psi_j^2)) / 2;
        c1 = 1 + sqrt(1 - 4*psi_j^2);

        Z2 = -((exp(lambdaH) - 1) / 1);

        % Base value
        F1Base = hypergeom([a1, b1], c1, Z2);
        Denominator = - psi_j^2 * (lambdaH/H)^2 * exp(2*lambdaH) * ...
                      (exp(lambdaH) - 1)^(S1 - 2) * F1Base;

        % Strain ratio profile
        GAMMA_ratio = NaN(1, nZ);
        for i = 1:nZ
            Zi = -((exp(lambdaH * z(i)/H) - 1) / 1);
            F1Depth = hypergeom([a1, b1], c1, Zi);

            Numerator = - psi_j^2 * (lambdaH/H)^2 * exp(2 * lambdaH * z(i)/H) * ...
                        (exp(lambdaH * z(i)/H) - 1)^(S1 - 2) * F1Depth;

            GAMMA_ratio(i) = Numerator / Denominator;
        end

        STRAIN_MATRIX(:, j) = GAMMA_ratio.';
    end

    % ====== PLOT ======
    ax = nexttile;
    hold(ax,'on');

    % curves
    for j = 1:nModes
        plot(ax, STRAIN_MATRIX(:, j), z_H, ...
             'LineWidth',1.5, 'LineStyle',lineStyles{j});
    end

    % axes style
    set(ax,'YDir','reverse', ...          
           'XAxisLocation','top', ...
           'XLim',[-8 8], ...
           'YLim',[0 1], ...
           'YTick',0:0.2:1, ...
           'FontSize',12);

    titles = { ...
        sprintf('\\lambdaH = %g', lambdaH_list(1)), ...
        sprintf('\\lambdaH = %g', lambdaH_list(2)), ...
        sprintf('\\lambdaH = %g', lambdaH_list(3)), ...
        sprintf('\\lambdaH = %g', lambdaH_list(4))};

    xlabel(ax, titles{lambdaIdx});




    box(ax,'on');
end

toc;
