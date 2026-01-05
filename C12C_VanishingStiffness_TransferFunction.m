clear
close all
clc

% Soil and wave parameters
H = 50;  % depth (m)
Vinf = 200; % m/s
alpha = 0.9995; % Fixed alpha

C_m = [sqrt(0.0829), sqrt(0.0366), sqrt(0.024), sqrt(0.0109)];

% Frequency range
f=linspace(0,80,2E3);
LENGTH = length(f);

% Parameter ranges
lambda_vals = [5/H, 8/H, 10/H, 15/H];  % Four lambda values
xi_vals = [0.03, 0.05, 0.1];       % Three damping values

% Create all (lambda, zeta) combinations
[LambdaGrid, xiGrid] = meshgrid(lambda_vals, xi_vals);
lambda_list = LambdaGrid(:);
xi_list = xiGrid(:);
num_cases = length(lambda_list);

% Preallocate results
TRANSFER_functions = cell(num_cases, 1);
f_normalized_list = cell(num_cases, 1);

for m=1:length(C_m)
for k = 1:num_cases
    lam = lambda_list(k);
    xi = xi_list(k);
    Vinf_complex = Vinf * sqrt(1 + (2 * 1i * xi));   
    f1n=(Vinf_complex / (4 * H)) * (2 / pi) * lam * H * C_m(m);
    
    f_normalized = f ./ f1n;

    
    z1 = -((1 - alpha) / alpha);
    z2 = -((exp(lam * H) - alpha) / alpha);

    TRANSFERF = zeros(1, LENGTH);

    parfor j = 1:LENGTH
        k_infty = (2 * pi * f_normalized(j)) / Vinf_complex;
        psi = (k_infty) / (lam);
        sq = sqrt(1 - 4 * psi^2);
        rt = 2 * 1i * psi;
        
        S1 = 0.5 * (1 + sq);
        S2 = 0.5 * (1 - sq);

        a1 = (1 - rt + sq) / 2;
        b1 = (1 + rt + sq) / 2;
        c1 = 1 + sq;
        a2 = (1 - rt - sq) / 2;
        b2 = (1 + rt - sq) / 2;
        c2 = 1 - sq;
        a3 = a1 + 1;
        b3 = b1 + 1;
        c3 = c1 + 1;
        a4 = a2 + 1;
        b4 = b2 + 1;
        c4 = c2 + 1;

        F1base = hypergeom([a1, b1], c1, z2);
        F1surf = hypergeom([a1, b1], c1, z1);
        F2base = hypergeom([a2, b2], c2, z2);
        F2surf = hypergeom([a2, b2], c2, z1);

        PH1base = ((-lam * exp(lam * H)) / (2 * alpha)) * hypergeom([a3, b3], c3, z2);
        PH1surf = ((-lam) / (2 * alpha)) * hypergeom([a3, b3], c3, z1);
        PH2base = ((-lam * exp(lam * H)) / (2 * alpha)) * hypergeom([a4, b4], c4, z2);
        PH2surf = (-lam / (2 * alpha)) * hypergeom([a4, b4], c4, z1);

        
        
        E_SURF = lam * 1 * S1 * ...
            (1 - alpha)^(-S2) * F1surf + ...
            (1 - alpha)^(S1) * PH1surf;

        E_BASE = lam * exp(lam * H) * S1 * ...
            (exp(lam * H) - alpha)^(-S2) * F1base + ...
            (exp(lam * H) - alpha)^(S1) * PH1base;
        
        M_SURF = lam * 1 * (S2) * ...
            (1 - alpha)^(-S1) * F2surf + ...
            (1 - alpha)^(S2) * PH2surf;

        M_BASE = lam * exp(lam * H) * (S2) * ...
            (exp(lam * H) - alpha)^(-S1) * F2base + ...
            (exp(lam * H) - alpha)^(S2) * PH2base;

        CONSTANT = (1 - alpha)^(sq) * (F1surf/F2surf);

        TRANSFERF(j) = (E_SURF - CONSTANT*M_SURF) / ...
        (E_BASE - CONSTANT*M_BASE);
    end

    TRANSFER_functions{k} = abs(TRANSFERF);
    f_normalized_list{k} = f_normalized;
end
end

% Plotting

for i = 1:length(lambda_vals)
    subplot(2,2,i);
    hold on;

    for j = 1:length(xi_vals)
        idx = (i - 1) * length(xi_vals) + j;
        plot(f_normalized_list{idx}, TRANSFER_functions{idx}, ...
             'DisplayName', sprintf('\\xi = %.2f', xi_vals(j)));
    end

    title(sprintf('\\lambdaH = %.2f', lambda_vals(i) * H), 'FontWeight', 'bold');
    xlabel('f/f_{1}');
    ylabel('|F(\omega)|');
    legend('Location', 'best');
    xlim([0 6])
end






