clear 
close all
clc
tic;
% Data for Site 1
z1 = [0 1 16 24 34]; 
Vs1 = [170 170 280 400 600];

%Elbasheir Site 1
alpha_1 = 0.9; 
lambda_1 = 0.0116;
Vinf_1 = 1365;
V0_1 = Vinf_1 *(1-alpha_1);

%Rovithis Site 1
V0_RM1 = 197;        
ratio_RM1 = 0.303;
alpha_RM1 = log(1 / ratio_RM1);




% Data for Site 2
z2 = [0 0.42 0.84 1.27 1.7 2.34 2.98 3.62 4.26 5.97 7.68 9.61 11.79 15 18.21 22.51 27.51 29.51 33.51 50 60]; 
Vs2 = [112 112 135 159 165 165 165 165 165 130 130 130 130 184 184 184 257 232 300 300 550]; 

%Elbasheir Site 2
V0_2 = 114.59; 
alpha_2 = 0.7916; 
lambda_2 = 0.0153;

%Rovithis Site 2
V0_RM2 = 134;        
ratio_RM2 = 0.265;
alpha_RM2 = log(1 / ratio_RM2);

% Exponential fit function for Site 1
exp_fit1 = @(z) V0_1 * (1/(1 - alpha_1)) * (1 - alpha_1 * exp(-lambda_1 * z));
exp_fit_RM1 = @(z) V0_RM1 *  exp(alpha_RM1 * z /34);

% Exponential fit function for Site 2
exp_fit2 = @(z) V0_2 * (1/(1 - alpha_2)) * (1 - alpha_2 * exp(-lambda_2 * z));
exp_fit_RM2 = @(z) V0_RM2 *  exp(alpha_RM2 * z /60);


% Generating finer z values for smooth exponential curve plotting
z_fine1 = linspace(min(z1), max(z1), 100);
Vs_fine1 = exp_fit1(z_fine1);
Vs_fine_RM1 = exp_fit_RM1(z_fine1);


z_fine2 = linspace(min(z2), max(z2), 100);
Vs_fine2 = exp_fit2(z_fine2);
Vs_fine_RM2 = exp_fit_RM2(z_fine2);


% Create the figure for Profile 1
figure;
hold on;

% Plot the stair-step (staircase) graph for Profile 1
stairs(Vs1, z1, 'DisplayName', 'Multi-layer Profile');

% Plot the exponential fit for Profile 1
plot(Vs_fine1, z_fine1, '--r', 'DisplayName', 'Exponential Fit');
plot(Vs_fine_RM1, z_fine1, ':k', 'DisplayName', 'Rovithis and Mylonakis 2022');

% Axis formatting for Profile 1
set(gca, 'XAxisLocation', 'top');
set(gca, 'YDir', 'reverse'); % Invert the y-axis
xlabel('V_{s} (m/s)');
ylabel('z (m)');
title('Site 1');
legend('Location', 'best');
hold off;

% Create the figure for Profile 2
figure;
hold on;

% Plot the stair-step (staircase) graph for Profile 2
stairs(Vs2, z2,'DisplayName', 'Multi-layer Profile');

% Plot the exponential fit for Profile 2
plot(Vs_fine2, z_fine2, '--r', 'DisplayName', 'Exponential Fit');
plot(Vs_fine_RM2, z_fine2, ':k', 'DisplayName', 'Rovithis and Mylonakis 2022');

% Axis formatting for Profile 2
set(gca, 'XAxisLocation', 'top');
set(gca, 'YDir', 'reverse'); % Invert the y-axis
xlabel('V_{s} (m/s)');
ylabel('z (m)');
title('Site 2');
legend('Location', 'best');
hold off;

toc;