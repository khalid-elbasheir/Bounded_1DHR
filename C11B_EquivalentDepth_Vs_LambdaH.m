clear;
close all;
clc;

% constants
H     = 50;
xi  = 0.05;
Vzero = 100;
Vzero_complex = Vzero * sqrt(1 + (2 * 1i * xi)); 

% parameters
alphas   = [-0.9, -0.5, 0.5, 0.9];
lambda   = 0.0005:0.01:0.3;          
lambda_H = lambda .* H;

% Fine grid for smoothing with spline
lambda_H_fine = linspace(lambda_H(1), lambda_H(end), 100);

figure('Color','w');
for a_idx = 1:length(alphas)
    alpha = alphas(a_idx);

    % Arrays to store equivalent depths
    zeq_H_high = zeros(length(lambda), 1);
    zeq_H_low  = zeros(length(lambda), 1);
    fsin       = zeros(1, length(lambda));

    % Precompute terms
    Vinf         = Vzero / (1 - alpha);
    Vinf_complex = Vinf * sqrt(1 + (2 * 1i * xi));
    fhomo        = Vinf_complex / (4 * H);

    
    for idx = 1:length(lambda)
        lam = lambda(idx);
        fsin(idx) = (H * sqrt((exp(-2*H*lam) * ...
            (4 * exp(H*lam) * alpha * (pi^4 + 6 * H^2 * pi^2 * lam^2 + 8 * H^4 * lam^4) ...
            - alpha^2 * (pi^4 + 9 * H^2 * pi^2 * lam^2 + 8 * H^4 * lam^4) ...
            + exp(2 * H * lam) * (8 * H^5 * lam^5 + pi^4 * ((-4 + alpha) * alpha + 2 * H * lam) ...
            + H^2 * pi^2 * lam^2 * ((-16 + alpha) * alpha + 10 * H * lam)))) ...
            / (H^3 * lam * (pi^4 + 5 * H^2 * pi^2 * lam^2 + 4 * H^4 * lam^4)))) / sqrt(2);
    end

   
    for i = 1:length(lambda)
        lh = lambda_H(i);

        % High frequencies
        zeq_H_high(i) = log( ...
            alpha - ((lh * alpha) / (lh + log(1 - alpha) - log(exp(lh) - alpha))) ) ...
            / lh;

        % Low frequencies
        zeq_H_low(i)  = log( (fhomo * alpha) / (fhomo - fsin(i)*fhomo) ) / lh;
    end

    % Spline
    zeq_H_high_s = spline(lambda_H, zeq_H_high, lambda_H_fine);
    zeq_H_low_s  = spline(lambda_H, zeq_H_low,  lambda_H_fine);

    % ---- Plotting ----
    subplot(2,2,a_idx);
    hold on;


    plot(lambda_H_fine, zeq_H_high_s, 'k-',  'LineWidth', 1.5); 
    plot(lambda_H_fine, zeq_H_low_s,  'k--', 'LineWidth', 1.5); 

    
    xlim([0 15]);
    ylim([0 1]);
    xlabel('\lambda H','Interpreter','tex');
    ylabel('z_{eq}/H','Interpreter','tex');

    set(gca, ...
             'FontSize', 12, ...
             'LineWidth',1, ...
             'Box','on', ...
             'TickDir','in');

    hold off;
end