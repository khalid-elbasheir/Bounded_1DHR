clearvars; close all; clc;

% ----- Given parameters -----
H = 50;  % m
alpha_vec    = [0.9 0.9 0.8 0.5 0.5 -0.5 -0.9 -1 -3];
lambdaH_vec  = [3   1   3   8   3    8    5    3  10];

% Depth vector (physical z for the function), and normalized depth for plotting
z  = linspace(0, H, 1000);     % m
zn = z / H;                    % normalized depth z/H

% Preallocate for min/max x-limits calculation
all_x = [];

figure(1); clf; hold on; box on;

% Plot each (alpha, lambdaH) pair
for k = 1:numel(alpha_vec)
    alpha    = alpha_vec(k);
    lambdaH  = lambdaH_vec(k);
    lambda   = lambdaH / H;           % inverse length

    % Vs/Vinf = 1 - alpha * exp(-lambda * z)  (use physical z here)
    Vs_over_Vinf = 1 - alpha .* exp(-lambda .* z);

    % Store for global x-limits
    all_x = [all_x, Vs_over_Vinf];

    % Plot Vs/Vinf vs normalized depth
    plot(Vs_over_Vinf, zn, 'LineWidth', 1, ...
        'DisplayName', sprintf('\\alpha=%.2g, \\lambdaH=%.2g', alpha, lambdaH));
end

% Axis formatting (match your old style: x on top, depth downward)
set(gca, 'XAxisLocation', 'top');
set(gca, 'YDir', 'reverse');

% Labels
xlabel('V_s(z)/V_{\infty}');
ylabel('z/H');

% % Tight x-limits from data (with a small padding)
% xmin = min(all_x); xmax = max(all_x);
% pad  = 0.02 * (xmax - xmin + eps);
% xlim([xmin - pad, xmax + pad]);

% y-limits 0..1 (top=surface, bottom=base)
ylim([0 1]);

legend('Location','bestoutside');
hold off;
