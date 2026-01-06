%% MULTI LAYER SOLUTION PROFILE 1 (GAZETAS AND DOBRY)
clear
close all
clc
tic;

f = linspace(0, 50, 10E3); % Adjust frequency range
LENGTH = length(f);

%Circular Frequency
omega = (2 * pi) .* f;

% Given data
xi = 0.05;
rho =2.0; % Mg/m^3
H = 60; % m
VZERO =114.59; % m/s109.23
VZERO_COMPLEX = VZERO * (1 + (1i * xi));

lambda_FIT = 0.0153;
alpha_FIT = 0.7916;
A_FIT = 1/(1-alpha_FIT);
VINF_COMPLEX = A_FIT * VZERO_COMPLEX;

% DATA of Rovithis and Mylonakis
V0_RM = 134;     % m/s
V0_COMPLEX_RM = V0_RM * (1 + (1i * xi));
ratio_RM= 0.265; 
alpha_RM= log(1/ratio_RM);

% Preallocate arrays
TRANSFER_ABSOLUTE_EXP = ones(1, length(f));
psi = zeros(1, length(f));
a1 = zeros(1, length(f)); b1 = zeros(1, length(f)); c1 = zeros(1, length(f));
a2 = zeros(1, length(f)); b2 = zeros(1, length(f)); c2 = zeros(1, length(f));
a3 = zeros(1, length(f)); b3 = zeros(1, length(f)); c3 = zeros(1, length(f));
a4 = zeros(1, length(f)); b4 = zeros(1, length(f)); c4 = zeros(1, length(f));
F1baseAsymp = zeros(1, length(f)); F1surfAsymp = zeros(1, length(f));
F2baseAsymp = zeros(1, length(f)); F2surfAsymp = zeros(1, length(f));
PH1baseAsymp = zeros(1, length(f)); PH1surfAsymp = zeros(1, length(f));
PH2baseAsymp = zeros(1, length(f)); PH2surfAsymp = zeros(1, length(f));
D = zeros(1, length(f)); E = zeros(1, length(f)); M = zeros(1, length(f));
TRANSFERF = zeros(1, length(f));
F_RM_Bessel = zeros(1, length(f));

% Exponential Fit Approach
parfor j = 1:length(f)
    psi(j) = (omega(j))^2 / (lambda_FIT^2 * VINF_COMPLEX^2);

    a1(j) = (1 - 2*1i*sqrt(psi(j)) + sqrt(1 - 4*psi(j))) / 2;
    b1(j) = (1 + 2*1i*sqrt(psi(j)) + sqrt(1 - 4*psi(j))) / 2;
    c1(j) = 1 + sqrt(1 - 4*psi(j));

    a2(j) = (1 - 2*1i*sqrt(psi(j)) - sqrt(1 - 4*psi(j))) / 2;
    b2(j) = (1 + 2*1i*sqrt(psi(j)) - sqrt(1 - 4*psi(j))) / 2;
    c2(j) = 1 - sqrt(1 - 4*psi(j));

    a3(j) = (3 - 2*1i*sqrt(psi(j)) + sqrt(1 - 4*psi(j))) / 2;
    b3(j) = (3 + 2*1i*sqrt(psi(j)) + sqrt(1 - 4*psi(j))) / 2;
    c3(j) = 2 + sqrt(1 - 4*psi(j));

    a4(j) = (3 - 2*1i*sqrt(psi(j)) - sqrt(1 - 4*psi(j))) / 2;
    b4(j) = (3 + 2*1i*sqrt(psi(j)) - sqrt(1 - 4*psi(j))) / 2;
    c4(j) = 2 - sqrt(1 - 4*psi(j));

    z1 = -((1 - alpha_FIT) / alpha_FIT);
    z2 = -((exp(lambda_FIT * H) - alpha_FIT) / alpha_FIT);

    F1baseAsymp(j) = hypergeom([a1(j), b1(j)], c1(j), z2);
    F1surfAsymp(j) = hypergeom([a1(j), b1(j)], c1(j), z1);
    F2baseAsymp(j) = hypergeom([a2(j), b2(j)], c2(j), z2);
    F2surfAsymp(j) = hypergeom([a2(j), b2(j)], c2(j), z1);

    PH1baseAsymp(j) = ((-lambda_FIT * exp(lambda_FIT * H)) / (2 * alpha_FIT)) * hypergeom([a3(j), b3(j)], c3(j), z2);
    PH1surfAsymp(j) = ((-lambda_FIT) / (2 * alpha_FIT)) * hypergeom([a3(j), b3(j)], c3(j), z1);
    PH2baseAsymp(j) = ((-lambda_FIT * exp(lambda_FIT * H)) / (2 * alpha_FIT)) * hypergeom([a4(j), b4(j)], c4(j), z2);
    PH2surfAsymp(j) = (-lambda_FIT / (2 * alpha_FIT)) * hypergeom([a4(j), b4(j)], c4(j), z1);

    D(j) = ((lambda_FIT * sqrt(1-4*psi(j))) / (1 - alpha_FIT)) +  (PH1surfAsymp(j)/F1surfAsymp(j)) - (PH2surfAsymp(j) /F2surfAsymp(j));
    E(j) = lambda_FIT * exp(lambda_FIT * H) * ((1 + sqrt(1 - 4*psi(j))) / 2) * ...
           (exp(lambda_FIT * H) - alpha_FIT)^((-1 + sqrt(1 - 4*psi(j))) / 2) * F1baseAsymp(j) + ...
           (exp(lambda_FIT * H) - alpha_FIT)^((1 + sqrt(1 - 4*psi(j))) / 2) * PH1baseAsymp(j);

    M(j) = lambda_FIT * exp(lambda_FIT * H) * ((1 - sqrt(1 - 4*psi(j))) / 2) * ...
           (exp(lambda_FIT * H) - alpha_FIT)^((-1 - sqrt(1 - 4*psi(j))) / 2) * F2baseAsymp(j) + ...
           (exp(lambda_FIT * H) - alpha_FIT)^((1 - sqrt(1 - 4*psi(j))) / 2) * PH2baseAsymp(j);

    TRANSFERF(j) = ((1 - alpha_FIT)^((1 + sqrt(1 - 4*psi(j)))/2)*F1surfAsymp(j) * D(j)) /  ...
                   (E(j) - (1 - alpha_FIT)^sqrt(1-4*psi(j)) * F1surfAsymp(j) / F2surfAsymp(j) * M(j));
    TRANSFER_ABSOLUTE_EXP(j) = abs(TRANSFERF(j));
end

% Rovithis and Mylonakis Approach
parfor j = 1:length(f)
    k0 = omega(j) / V0_COMPLEX_RM;

    A = (k0 * H / alpha_RM) * exp(-alpha_RM);
    B = (k0 * H / alpha_RM);

    % Bessel function terms
    term1 = besselj(1, A) * bessely(0, B);
    term2 = besselj(0, B) * bessely(1, A);

    numerator = (2 / pi);
    inner = A * (term1 - term2);

        F_RM_Bessel(j) = numerator / inner;

end

% Multi-Layer Transfer Matrix Approach
z = [0.42 0.84 1.27 1.7 2.34 2.98 3.62 4.26 5.97 7.68 9.61 11.79 15 18.21 22.51 27.51 29.51 33.51 50 60];
Vs = [112 135 159 165 165 165 165 165 130 130 130 130 184 184 184 257 232 300 300 550];
H_layers = diff([0, z]);

% Calculate alpha for each layer
Vs_complex = zeros(1, length(H_layers));
G = zeros(1, length(H_layers));

for p = 1:length(H_layers)
    Vs_complex(p) = Vs(p) * (1 + (1i * xi));
    G(p) = rho*(Vs_complex(p)^2);    
end

% Vs30 Calculation
depth_limit = 30; % depth of interest in meters (Vs30)
cumulative_depth = 0;
inv_sum_vs = 0;
total_thickness = 0;

for w = 1:length(H_layers)
    if cumulative_depth + H_layers(w) <= depth_limit
        inv_sum_vs = inv_sum_vs + H_layers(w) / Vs(w);
        cumulative_depth = cumulative_depth + H_layers(w);
        total_thickness = total_thickness + H_layers(w);
    else
        remaining_depth = depth_limit - cumulative_depth;
        inv_sum_vs = inv_sum_vs + remaining_depth / Vs(w);
        total_thickness = total_thickness + remaining_depth;
        break;
    end
end

Vs30 = total_thickness / inv_sum_vs; % Harmonic mean for Vs30
Vs30_complex = Vs30 * (1 + (1i * xi));
TRANSFER_VS30 = zeros(1, length(f));

for k = 1:length(f)
    TRANSFER_VS30(k) = abs(1 / cos(omega(k) * 30 / Vs30_complex));
end

% Preallocate arrays and matrices for Multi-Layer Approach
TRANSFER_ABSOLUTE_ML = zeros(1, length(f));
KAPPA = zeros(1, length(f));

TRANSFERF_SCALAR = zeros(1, length(f));

for g = 1:length(f)
    % Initialize Pimatrix for each frequency
    Pimatrix = eye(2);

    % Handling the zero frequency case separately
    if f(g) == 0
        % Set a static response for zero frequency
       TRANSFERF_SCALAR(g) = 1.0; % Static case: no dynamic response
        TRANSFER_ABSOLUTE_ML(g) = 1.0; % Absolute value is also 1 for static case
    else
        for i = 1:length(H_layers)
            % Compute KAPPA for each layer
            KAPPA(g) = (omega(g)/Vs_complex(i));

            % Define DeltaInv for the current layer (no inversion needed)
            DeltaInv = [cos(KAPPA(g)*H_layers(i)), sin(KAPPA(g)*H_layers(i))/(G(i)*KAPPA(g)); 
                               - (G(i) * KAPPA(g)) * sin(KAPPA(g) * H_layers(i)), cos(KAPPA(g) * H_layers(i))];
            
            % Multiply the transfer matrix of this layer to Pimatrix
            Pimatrix = Pimatrix * DeltaInv;  % No pinv or inv
        end
        
        % Calculate the transfer function scalar for non-zero frequencies
        TRANSFERF_SCALAR(g) = Pimatrix(1,1) - (Pimatrix(1,2) * Pimatrix(2,1) / Pimatrix(2,2));
       
        TRANSFER_ABSOLUTE_ML(g) = abs(TRANSFERF_SCALAR(g));
    end
end
 
% Plotting All Solutions for Comparison

figure (1);
hold on;
plot(f, TRANSFER_ABSOLUTE_EXP, 'k--', 'LineWidth', 1); % Exponential Fit Solution
plot(f, abs(F_RM_Bessel), 'k:', 'LineWidth', 1); 
plot(f, TRANSFER_ABSOLUTE_ML, 'k-', 'LineWidth', 1); % Multi-Layer Solution
plot(f, TRANSFER_VS30, 'g-', 'LineWidth', 1); % Vs30 Approximation Solution
set(gca, 'FontSize', 12);
xlabel('f (Hz)');
ylabel('F(\omega)');
legend('Exponential Fit (Eq. 30)','Rovithis and Mylonakis (2022)', 'Multi-layer Profile', 'Hom. Profile Approx. (H=30m, V_{s}=V_{s30})');
xlim([0 16]);
hold off;

max(TRANSFER_ABSOLUTE_ML)
max(TRANSFER_ABSOLUTE_EXP)
max(abs(TRANSFER_VS30))
max(abs(F_RM_Bessel))

% Find the first three maxima for Exponential Fit Transfer Function
[pks_exp, locs_exp] = findpeaks(TRANSFER_ABSOLUTE_EXP, f, 'SortStr', 'descend', 'NPeaks', 3);
disp('Exponential Fit:');
disp(table(pks_exp(:), locs_exp(:), 'VariableNames', {'Amplitude', 'Frequency'}));

% Find the first three maxima for Multi-Layer Transfer Function
[pks_ml, locs_ml] = findpeaks(TRANSFER_ABSOLUTE_ML, f, 'SortStr', 'descend', 'NPeaks', 3);
disp('Multi-Layer Transfer Matrix:');
disp(table(pks_ml(:), locs_ml(:), 'VariableNames', {'Amplitude', 'Frequency'}));

% Find the first three maxima for Vs30 Approximation Transfer Function
[pks_vs30, locs_vs30] = findpeaks(TRANSFER_VS30, f, 'SortStr', 'descend', 'NPeaks', 3);
disp('Vs30 Approximation:');
disp(table(pks_vs30(:), locs_vs30(:), 'VariableNames', {'Amplitude', 'Frequency'}));

% Find the first three maxima for Vs30 Approximation Transfer Function
[pks_RM, locs_RM] = findpeaks(abs(F_RM_Bessel), f, 'SortStr', 'descend', 'NPeaks', 3);
disp('Rovithis and Mylonakis 2022:');
disp(table(pks_RM(:), locs_RM(:), 'VariableNames', {'Amplitude', 'Frequency'}));

toc;