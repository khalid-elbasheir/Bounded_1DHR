clear;
close all;
clc;

% Constants
H = 50;
xi = 0.05;
Vzero = 100; 
Vzero_complex = Vzero * sqrt(1 + (2 * 1i * xi));

% Parameters
alpha = -10:1e-5:1;        
lambdaH_vec = [3, 5, 10, 15];

% Precompute terms
A      = 1 ./ (1 - alpha);
ratio  = 1 - alpha;
fhomo  = (Vzero_complex .* A) ./ (4 * H);


% Figure:  
figure(1); 
clf;
for k = 1:numel(lambdaH_vec)
     lambda_H = lambdaH_vec(k);     
     lambda   = lambda_H / H;       

     fsin   = (H .* sqrt((exp(-2 * H * lambda) .* ...
        (4 .* exp(H * lambda) .* alpha .* (pi^4 + 6 * H^2 * pi^2 * lambda^2 + 8 * H^4 * lambda^4) ...
        - alpha.^2 .* (pi^4 + 9 * H^2 * pi^2 * lambda^2 + 8 * H^4 * lambda^4) ...
        + exp(2 * H * lambda) .* (8 * H^5 * lambda^5 + pi^4 .* ((-4 + alpha) .* alpha + 2 * H * lambda) ...
        + H^2 * pi^2 * lambda^2 .* ((-16 + alpha) .* alpha + 10 * H * lambda)))) ./ ...
        (H^3 * lambda .* (pi^4 + 5 * H^2 * pi^2 * lambda^2 + 4 * H^4 * lambda^4)))) ./ sqrt(2);

    den_common = fhomo - (fsin .* fhomo);


    % High-frequency
    zeq_H_high = log(alpha + ((lambda_H .* alpha) ./ (-lambda_H - log(-1 + alpha) + log(-exp(lambda_H) + alpha)))) ./ lambda_H;
    zeq_H_high(isnan(zeq_H_high)) = -1e6;

    % Low-frequency
    den = den_common;
    near_zero = abs(den) < eps;        
    zeq_H_low = zeros(size(alpha));
    zeq_H_low(~near_zero) = log((fhomo(~near_zero) .* alpha(~near_zero)) ./ den(~near_zero)) ./ lambda_H;
    zeq_H_low(near_zero) = -1e6;        



    % Plot
    subplot(2,2,k); hold on; box on;
    plot(alpha, real(zeq_H_high), 'k-',  'DisplayName','High Frequency');
    plot(alpha, real(zeq_H_low),  'k--', 'DisplayName','Low Frequency');
    xlabel('\alpha'); ylabel('z_{eq}/H');
    title(sprintf('\\lambdaH = %g', lambda_H));
    xlim([-1 1]); ylim([0 1]);
    legend('show','Location','southwest');
end