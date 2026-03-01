clear
close all
clc

% parameters
xi = 0.05;
H = 50;  
Vinf = 200; 
lambda_vals = [3/H 1/H 3/H 8/H 3/H 1/H];
alpha = [0.9 0.9 0.8 0.5 0.5 1E-5]; 
C_m = [sqrt(0.1864) sqrt(0.5909) sqrt(0.1997) sqrt(0.0378) sqrt(0.2308) sqrt(2.4662)];%Charestrstic Equation first Root

%%% For decreasing stiffness profile remove the comment from the code lines below
%lambda_vals = [8/H 5/H 3/H 10/H 1/H];
%alpha = [-0.5 -0.9 -1 -3 1E-5]; 
%C_m = [sqrt(0.0391) sqrt(0.1073) sqrt(0.3531) sqrt(0.0256) sqrt(2.4662)];

% Frequency vector
f=linspace(0, 80, 2E3);
LENGTH = length(f);

     
num_cases = length(alpha);

% Preallocate results
TRANSFER_functions = cell(num_cases, 1);
f_normalized_list = cell(num_cases, 1);


for k = 1:num_cases
    a = alpha(k);
    lam = lambda_vals(k);
    CharqRoot = C_m(k);
    Vinf_complex = Vinf * sqrt(1 + (2 * 1i * xi));   
    f1n = (Vinf_complex / (4 * H)) * (2 / pi) * lam * H * CharqRoot;
    
    f_normalized = f ./ f1n;

    
    z1 = -((1 - a) / a);
    z2 = -((exp(lam * H) - a) / a);

    TRANSFERF = zeros(1, LENGTH);

    parfor j = 1:LENGTH
        k_infty = (2 * pi * f(j)) / Vinf_complex;
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

        PH1base = ((-lam * exp(lam * H)) / (2 * a)) * hypergeom([a3, b3], c3, z2);
        PH1surf = ((-lam) / (2 * a)) * hypergeom([a3, b3], c3, z1);
        PH2base = ((-lam * exp(lam * H)) / (2 * a)) * hypergeom([a4, b4], c4, z2);
        PH2surf = (-lam / (2 * a)) * hypergeom([a4, b4], c4, z1);

      
        E_SURF = lam * 1 * S1 * ...
            (1 - a)^(-S2) * F1surf + ...
            (1 - a)^(S1) * PH1surf;

        E_BASE = lam * exp(lam * H) * S1 * ...
            (exp(lam * H) - a)^(-S2) * F1base + ...
            (exp(lam * H) - a)^(S1) * PH1base;
        
        M_SURF = lam * 1 * (S2) * ...
            (1 - a)^(-S1) * F2surf + ...
            (1 - a)^(S2) * PH2surf;

        M_BASE = lam * exp(lam * H) * (S2) * ...
            (exp(lam * H) - a)^(-S1) * F2base + ...
            (exp(lam * H) - a)^(S2) * PH2base;

        CONSTANT = (1 - a)^(sq) * (F1surf/F2surf);

        TRANSFERF(j) = (E_SURF - CONSTANT*M_SURF) / ...
        (E_BASE - CONSTANT*M_BASE);
    end

    TRANSFER_functions{k} = abs(TRANSFERF);
    f_normalized_list{k} = f_normalized;
end




% Plotting

    figure;
    hold on; 

    for j = 1:length(alpha)
       
        plot(f_normalized_list{j}, TRANSFER_functions{j}, ...
             'LineWidth', 2);
        legendEntries = arrayfun(@(j) sprintf('\\lambdaH = %.2f, \\alpha = %.2f', ...
    lambda_vals(j)*H, alpha(j) ), 1:length(alpha), 'UniformOutput', false);
    end

    xlabel('f/f_{1}');
    ylabel('|F(\omega)|');
    legend(legendEntries, 'Location', 'Best');
    xlim([0 6])
    set(gca, 'FontSize', 12);







