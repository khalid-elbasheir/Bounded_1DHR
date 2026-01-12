% Clear variables and close figures
clear ;
close all;
clc;

% Define constants
H = 50;
zeta = 0.05;
Vzero = 100; % m/s
Vzero_complex = Vzero * sqrt(1 + (2 * 1i * zeta));

% Define alpha values and lambda range
alphas = [-0.9, -0.5, 0.5, 0.9];
lambda = 0.0005:0.01:0.3;
lambda_H = lambda .* H;

% Loop over each alpha and create subplots
figure;
for a_idx = 1:length(alphas)
    alpha = alphas(a_idx);
    
    % Initialize arrays to store equivalent depths
    zeq_H_high = zeros(length(lambda), 1);
    zeq_H_low = zeros(length(lambda), 1);
    fhomo = zeros(1, length(lambda));
    fsin = zeros(1, length(lambda));
    
    for idx = 1:length(lambda)
        Vinf = Vzero / (1 - alpha);
        Vinf_complex = Vinf * sqrt(1 + (2 * 1i * zeta));
        fhomo = Vinf_complex / (4 * H);
        fsin(idx) = (H * sqrt((exp(-2*H*lambda(idx)) * (4 * exp(H*lambda(idx)) * alpha * (pi^4 + 6 * H^2 * pi^2 * lambda(idx)^2 + 8 * H^4 * lambda(idx)^4) - alpha^2 * (pi^4 + 9 * H^2 * pi^2 * lambda(idx)^2 + 8 * H^4 * lambda(idx)^4) + exp(2 * H * lambda(idx)) * (8 * H^5 * lambda(idx)^5 + pi^4 * ((-4 + alpha) * alpha + 2 * H * lambda(idx)) + H^2 * pi^2 * lambda(idx)^2 * ((-16 + alpha) * alpha + 10 * H * lambda(idx))))) / (H^3 * lambda(idx) * (pi^4 + 5 * H^2 * pi^2 * lambda(idx)^2 + 4 * H^4 * lambda(idx)^4)))) / sqrt(2);
    end

    for i = 1:length(lambda)
        % Define the function for NormalizedEqDepth_HighFrequency
        zeq_H_high(i) =   (log(alpha - ((lambda_H(i) * alpha) / (lambda_H(i) + log(1 - alpha) - log(exp(lambda_H(i)) - alpha)))) / (lambda_H(i)));
         zeq_H_low(i) =  log( (fhomo * alpha) / (fhomo - fsin(i)*fhomo) ) / lambda_H(i);
        %zeq_H_low(i) = ( log(fhomo * alpha) - log(fhomo - fsin(i)) ) / lambda_H(i);
    end

    % Create subplot
    subplot(2, 2, a_idx);
    hold on;
    plot(lambda_H, zeq_H_high, 'k-', 'DisplayName', 'High Frequency');
    plot(lambda_H, zeq_H_low, 'k--', 'DisplayName', 'Low Frequency (small \alpha)');
 
   
    hold off;
    xlabel('\lambda H');
    ylabel('z_{eq}/H');
    legend('show');
    title(['\alpha = ', num2str(alpha)]);
    ylim([0 1]); % Adjust the y-limits to show the trend
end
