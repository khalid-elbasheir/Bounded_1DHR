% Clear variables and close figures
clear 
close all
clc

% Soil Deposit
xi = 0.05;
H = 50; 
Vinf = 200; % m/s

% Excitation Frequencies
omega = linspace(0, 900, 10E3);% Increase the Number of points for accurate values
LENGTH = length(omega);
f = omega / (2 * pi);

lambda = 3/H; 
alpha = [0.5 0.9 -0.5 -0.9];
CharqRoot = [0.2308 0.1864 0.3142 0.3453];


F1sin_Norm=ones(1, length(alpha)); 
for id=1:length(alpha)
 F1sin_Norm(id) = (H * sqrt((exp(-2*H*lambda) * (4 * exp(H*lambda) * alpha(id) * (pi^4 + 6 * H^2 * pi^2 * lambda^2 + 8 * H^4 * lambda^4) - alpha(id)^2 * (pi^4 + 9 * H^2 * pi^2 * lambda^2 + 8 * H^4 * lambda^4) + exp(2 * H * lambda) * (8 * H^5 * lambda^5 + pi^4 * ((-4 + alpha(id)) * alpha(id) + 2 * H * lambda) + H^2 * pi^2 * lambda^2 * ((-16 + alpha(id)) * alpha(id) + 10 * H * lambda)))) / (H^3 * lambda * (pi^4 + 5 * H^2 * pi^2 * lambda^2 + 4 * H^4 * lambda^4)))) / sqrt(2);
end   

% Create figure for subplots
figure;

% Loop over the alpha values to create subplots
for idx = 1:length(alpha)
    a = alpha(idx);
    % Calculate Fundamental Frequency and WaveNumber Function
    Vinf_complex = Vinf * sqrt(1 + (2 * 1i * xi));   
    f1n = (Vinf_complex / (4 * H)) * (2 / pi) * lambda * H * CharqRoot(idx);
    f_Normalized = f / f1n;

    k_infty = (2 * pi * f) ./ Vinf_complex;
    psi = (k_infty) ./ (lambda);


    % Initialize arrays
    ps = ones(1, LENGTH); 
    a1 = ones(1, LENGTH); 
    b1 = ones(1, LENGTH); 
    c1 = ones(1, LENGTH); 
    a2 = ones(1, LENGTH);
    b2 = ones(1, LENGTH); 
    c2 = ones(1, LENGTH);
    a3 = ones(1, LENGTH);
    b3 = ones(1, LENGTH); 
    c3 = ones(1, LENGTH);
    a4 = ones(1, LENGTH);
    b4 = ones(1, LENGTH); 
    c4 = ones(1, LENGTH);
    z1 = -((1 - a) / a);
    z2 = -((exp(lambda * H) - a) / a);
    F1baseAsymp = ones(1, LENGTH);
    F1surfAsymp = ones(1, LENGTH);
    F2baseAsymp = ones(1, LENGTH);
    F2surfAsymp = ones(1, LENGTH);
    PH1baseAsymp = ones(1, LENGTH);
    PH1surfAsymp = ones(1, LENGTH);
    PH2baseAsymp = ones(1, LENGTH);
    PH2surfAsymp = ones(1, LENGTH);
    D = ones(1, LENGTH);
    E = ones(1, LENGTH);
    M = ones(1, LENGTH);
    TRANSFERF_exact = ones(1, LENGTH);
    TRANSFER_ABSOLUTE_exact = ones(1, LENGTH);
    TRANSFERF_HIGHF = zeros(1, LENGTH);
    TRANSFER_ABSOLUTE_HIGHF = zeros(1, LENGTH);
    TRANSFERF_LOWF = zeros(1, LENGTH);
    TRANSFER_ABSOLUTE_LOWF = zeros(1, LENGTH);

    % Exact Solution Calculation
    parfor j = 1:LENGTH
        ps = psi(j);        
        sq = sqrt(1 - 4 * ps^2);
        rt = 2 * 1i * ps;
        
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

        PH1base = ((-lambda * exp(lambda * H)) / (2 * a)) * hypergeom([a3, b3], c3, z2);
        PH1surf = ((-lambda) / (2 * a)) * hypergeom([a3, b3], c3, z1);
        PH2base = ((-lambda * exp(lambda * H)) / (2 * a)) * hypergeom([a4, b4], c4, z2);
        PH2surf = (-lambda / (2 * a)) * hypergeom([a4, b4], c4, z1);

      
        E_SURF = lambda * 1 * S1 * ...
            (1 - a)^(-S2) * F1surf + ...
            (1 - a)^(S1) * PH1surf;

        E_BASE = lambda * exp(lambda * H) * S1 * ...
            (exp(lambda * H) - a)^(-S2) * F1base + ...
            (exp(lambda * H) - a)^(S1) * PH1base;
        
        M_SURF = lambda * 1 * (S2) * ...
            (1 - a)^(-S1) * F2surf + ...
            (1 - a)^(S2) * PH2surf;

        M_BASE = lambda * exp(lambda * H) * (S2) * ...
            (exp(lambda * H) - a)^(-S1) * F2base + ...
            (exp(lambda * H) - a)^(S2) * PH2base;

        CONSTANT = (1 - a)^(sq) * (F1surf/F2surf);

        TRANSFERF_exact(j) = (E_SURF - CONSTANT*M_SURF) / ...
        (E_BASE - CONSTANT*M_BASE);

        TRANSFER_ABSOLUTE_exact(j) = abs(TRANSFERF_exact(j)); 
    end

    % High Frequency Solution 
    for j = 1:LENGTH
        TRANSFERF_HIGHF(j) = sqrt((1 - a * exp(-H * lambda))/(1 - a)) * 1/cos(psi(j) * (-log(1 - a) + log(exp(H * lambda) - a)));
        
        
        TRANSFER_ABSOLUTE_HIGHF(j) = abs(TRANSFERF_HIGHF(j));
    end

    % Low Frequency Solution 
    for j = 1:LENGTH
        F_sinsoidal = F1sin_Norm(idx) * Vinf_complex / (4 * H);
        TRANSFERF_LOWF(j) = 1/cos(2*pi*f(j)/(4*F_sinsoidal));
        TRANSFER_ABSOLUTE_LOWF(j) = abs(TRANSFERF_LOWF(j));
    end

    % Exclude the first resonance part of the HIGHF solution
    [~, first_peak_index_HIGHF] = max(TRANSFERF_HIGHF(1:floor(LENGTH/2)));
    second_resonance_index_HIGHF = first_peak_index_HIGHF + find(TRANSFERF_HIGHF(first_peak_index_HIGHF+1:end) < TRANSFERF_HIGHF(first_peak_index_HIGHF), 1);
    % Find the second resonance peak after the first resonance
    [~, second_peak_index_HIGHF] = max(TRANSFERF_HIGHF(second_resonance_index_HIGHF:end));
    second_peak_index_HIGHF = second_resonance_index_HIGHF + second_peak_index_HIGHF - 1;
    % Add an offset to start just a bit before the second resonance
    offset_HIGHF = 30; % Adjust this value as needed
    adjusted_start_index_HIGHF = max(1, second_peak_index_HIGHF + offset_HIGHF);


    % Exclude the first resonance part of the LOWF solution
    [~, first_peak_index_LOWF] = max(TRANSFERF_LOWF(1:floor(LENGTH/2)));
    second_resonance_index_LOWF = first_peak_index_LOWF + find(TRANSFERF_LOWF(first_peak_index_LOWF+1:end) < TRANSFERF_LOWF(first_peak_index_LOWF), 1);
    offset_LOWF = 30; % Adjust this value as needed
    adjusted_end_index_LOWF = second_resonance_index_LOWF + offset_LOWF;

    % Create subplot for each alpha
    subplot(2, 2, idx);
    hold on;
    plot(f_Normalized, abs(TRANSFERF_exact), 'k', 'LineWidth', 3, 'DisplayName', 'Exact Solution');
    plot(f_Normalized(adjusted_start_index_HIGHF:end), abs(TRANSFERF_HIGHF(adjusted_start_index_HIGHF:end)), 'k--', 'LineWidth', 3, 'DisplayName', 'Asymptotics - High Frequency');
    plot(f_Normalized(1:adjusted_end_index_LOWF), abs(TRANSFERF_LOWF(1:adjusted_end_index_LOWF)), 'k:', 'LineWidth', 3, 'DisplayName', 'Asymptotics - Low Frequency');
    hold off;
    set(gca, 'FontSize', 22);
    set(findall(gcf, 'Type', 'Text'), 'FontWeight', 'bold');
    xlabel('f/f_{1}', 'Interpreter', 'tex');
    ylabel('F(\omega)');
    legend('Location', 'Best', 'FontSize', 20);
    xlim([0 6]); % Adjust this range as necessary
    ylim([0 18]); % Adjust this range as necessary
    title(['\alpha = ' num2str(a)]);
end

set(gcf, 'WindowState', 'maximized');
