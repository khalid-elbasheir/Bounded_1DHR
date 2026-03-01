clear; 
close all; 
clc;

% parameters 
H = 50;  
alpha_vec    = [0.9 0.9 0.8 0.5 0.5 -0.5 -0.9 -1 -3];
lambdaH_vec  = [3   1   3   8   3    8    5    3  10];

% Depth vector 
z  = linspace(0, H, 1000);     
zn = z / H;                    % normalized depth

% Preallocate
all_x = [];

figure(1); 
clf; 
hold on; 
box on;

% Loop
for k = 1:numel(alpha_vec)
    alpha    = alpha_vec(k);
    lambdaH  = lambdaH_vec(k);
    lambda   = lambdaH / H;           
  
    Vs_over_Vinf = 1 - alpha .* exp(-lambda .* z);

    % Store
    all_x = [all_x, Vs_over_Vinf];

    % Plot
    plot(Vs_over_Vinf, zn, 'LineWidth', 1, ...
        'DisplayName', sprintf('\\alpha=%.2g, \\lambdaH=%.2g', alpha, lambdaH));
end

% Axis formatting
set(gca, 'XAxisLocation', 'top');
set(gca, 'YDir', 'reverse');

% Labels
xlabel('V_s(z)/V_{\infty}');
ylabel('z/H');


ylim([0 1]);

legend('Location','bestoutside');
hold off;