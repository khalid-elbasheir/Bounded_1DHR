% This code plots the results in figure 14
clear 
close all
clc
tic

H = 50;
lambda = [2e-04 0.01 0.02 0.04 0.06 0.08 0.1 0.12 0.14 0.16 0.18 0.2];
Alpha = [0.1  0.4  0.6  0.9];
sheetNames = {'Alpha 0.1',  'Alpha 0.4',  'Alpha 0.6',  'Alpha 0.9'};

data = cell(4, 1);

for i = 1:4
    data{i} = readtable('C9BC_RayleighPlots_FrequencyTables.xlsx', 'Sheet', sheetNames{i});
end

figure;


% Generate a finer grid for lambda
lambda_fine = linspace(min(lambda*H), max(lambda*H), 500);

for i = 1:4
    Finhomo_Linear = data{i}{:, 1};
    Finhomo_Parabolic = data{i}{:, 2};
    Finhomo_Sinsusoidal = data{i}{:, 3};
    Finhomo_Exact = data{i}{:, 4};

    subplot(2, 2, i);
    hold on;

    % Spline interpolation for each vector
    spline_Exact = spline(lambda*H, Finhomo_Exact, lambda_fine);
    spline_Linear = spline(lambda*H, Finhomo_Linear, lambda_fine);
    spline_Parabolic = spline(lambda*H, Finhomo_Parabolic, lambda_fine);
    spline_Sinsusoidal = spline(lambda*H, Finhomo_Sinsusoidal, lambda_fine);
    

    % Plotting the spline-interpolated data
    plot(lambda_fine, spline_Exact, 'r-', 'LineWidth', 2, 'DisplayName', 'Exact');
    plot(lambda_fine, spline_Linear, 'k-', 'LineWidth', 1, 'DisplayName', 'Linear');
    plot(lambda_fine, spline_Parabolic, 'K--', 'LineWidth', 1, 'DisplayName', 'Parabolic');
    plot(lambda_fine, spline_Sinsusoidal, 'K:', 'LineWidth', 1, 'DisplayName', 'Sinusoidal');
    
    
    hold off;
    legend('show', 'Location', 'best');
    xlabel('\lambdaH');
    ylabel('f_{1}/f_{1H}', 'Interpreter', 'tex');
    title(['\alpha = ' num2str(Alpha(i))]);
    
end



toc
