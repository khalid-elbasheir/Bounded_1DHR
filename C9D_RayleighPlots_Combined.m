clear 
close all
clc
tic

H = 50;
lambda = [2e-04 0.01 0.02 0.04 0.06 0.08 0.1 0.12 0.14 0.16 0.18 0.2];
Alpha = [-0.9 -0.5 -0.3 -0.1 0.1 0.4 0.6 0.9];
sheetNames = {'Alpha -0.9','Alpha -0.5','Alpha -0.3','Alpha -0.1', 'Alpha 0.1', 'Alpha 0.4', 'Alpha 0.6', 'Alpha 0.9'};

data = cell(8, 1);

for i = 1:8
    data{i} = readtable('C9BC_RayleighPlots_FrequencyTables.xlsx', 'Sheet', sheetNames{i});
end

figure;
set(gcf, 'Position', [100, 100, 1200, 900]); % Set figure size

Finhomo_Linear_All = zeros(length(lambda), length(Alpha));
Finhomo_Parabolic_All = zeros(length(lambda), length(Alpha));
Finhomo_Sinsusoidal_All = zeros(length(lambda), length(Alpha));
Finhomo_Exact_All = zeros(length(lambda), length(Alpha));

for i = 1:8
    Finhomo_Linear_All(:, i) = data{i}{:, 1};
    Finhomo_Parabolic_All(:, i) = data{i}{:, 2};
    Finhomo_Sinsusoidal_All(:, i) = data{i}{:, 3};
    Finhomo_Exact_All(:, i) = data{i}{:, 4};
end

% Generate a finer grid for Alpha for smooth curves
Alpha_fine = linspace(min(Alpha), max(Alpha), 500);

plot_indices = [3, 6, 10, 12];
subplot_positions = 1:4;

for k = 1:length(plot_indices)
    j = plot_indices(k);
    subplot(2, 2, subplot_positions(k));
    hold on;
    
    % Spline interpolation for each data vector
    spline_Exact = spline(Alpha, Finhomo_Exact_All(j, :), Alpha_fine);
    spline_Linear = spline(Alpha, Finhomo_Linear_All(j, :), Alpha_fine);
    spline_Parabolic = spline(Alpha, Finhomo_Parabolic_All(j, :), Alpha_fine);
    spline_Sinsusoidal = spline(Alpha, Finhomo_Sinsusoidal_All(j, :), Alpha_fine);
    
    
    % Plot the spline-interpolated curves
    plot(Alpha_fine, spline_Exact, 'r-', 'LineWidth', 2, 'DisplayName', 'Exact');
    plot(Alpha_fine, spline_Linear, 'k-', 'LineWidth', 1, 'DisplayName', 'Linear');
    plot(Alpha_fine, spline_Parabolic, 'k--', 'LineWidth', 1, 'DisplayName', 'Parabolic');
    plot(Alpha_fine, spline_Sinsusoidal, 'k:', 'LineWidth', 1, 'DisplayName', 'Sinusoidal');
    
    
    % Add dashed lines at alpha = 0 and y = 1
    plot([0 0], [0 1], 'k--', 'LineWidth', 0.8); % Vertical dashed line at alpha = 0
    plot([-1,0], [1, 1], 'k--', 'LineWidth', 0.8);

    % Add Marker at the intersection (0, 1)
    
    plot(0, 1, 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w'); 


    hold off;
    legend('show', 'Location', 'best');
    xlabel('\alpha');
    ylabel('f_{1}/f_{1H}', 'Interpreter', 'tex');
    title(['\lambdaH = ' num2str(lambda(j) * H)]);
  
end

toc