%% Surface Response - Site 1 - ChiChi
clear all
close all
clc
tic;

% Load the input motion
load ChiChi.txt  % Load Northridge earthquake data

LENGTH = 16384;

% Original time step for Northridge
dt = 0.005;  % Original time step for Northridge

BaseMotion= ChiChi;

% Normalize the base motion by its peak value
%BaseMotion = transpose(BaseMotion_RAW / max(abs(BaseMotion_RAW)));  % Normalize to gqs (peak normalization)

% Step 2: Frequency and time vectors for FFT analysis
df = 1 / (LENGTH * dt);  % Frequency increment in Hz

% Set up the time and frequency vectors
t = dt:dt:(LENGTH * dt);  % Time vector based on resampled data
f = df:df:(LENGTH * df);  % Full frequency vector

% Add zero time to time and frequency vectors
%BaseMotion = [0 BaseMotion];
t = [0 t];
f = [0 f];

% Compute FFT of the resampled base motion
BaseFFT = fft(BaseMotion, length(f));  % Compute the FFT for the resampled motion

%Circular Frequency
omega = (2 * pi) .* f;

% Given data
xi = 0.05;
rho =2.0; % Mg/m^3
H = 34; % m

lambda_FIT = 0.0116;
alpha_FIT = 0.9;
VINF= 1365;
VZERO = VINF * (1 - alpha_FIT);
VINF_COMPLEX = VINF * (1 + (1i * xi));
VZERO_COMPLEX = VZERO * (1 + (1i * xi));

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

% Multi-Layer Transfer Matrix Approach
z = [1 16 24 34];
Vs = [170 280 400 600];
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
    TRANSFER_VS30(k) = (1 / cos(omega(k) * 30 / Vs30_complex));
end
TRANSFER_VS30_ABSOLUTE=abs(TRANSFER_VS30);

% Multi-Layer Haskell Thompson Approach
TRANSFER_ABSOLUTE_ML = zeros(1, length(f));
Pimatrix = eye(2);
KAPPA = zeros(1, length(f));
TRANSFERF_SCALAR = zeros(1, length(f));

for g = 1:length(f)
    % Initialize Pimatrix for each frequency
    Pimatrix = eye(2);

    % Handling the zero frequency case separately
    if f(g) == 0
         %Set a static response for zero frequency
       TRANSFERF_SCALAR(g) = 1.0; % Static case: no dynamic response
        TRANSFER_ABSOLUTE_ML(g) = 1.0; % Absolute value is also 1 for static case
    else
        for i = 1:length(H_layers)
            % Compute KAPPA for each layer
            KAPPA(g) = (omega(g)/Vs_complex(i));

            % Define DeltaInv for the current layer (no inversion needed)
            DeltaInv= [cos(KAPPA(g)*H_layers(i)), sin(KAPPA(g)*H_layers(i))/(G(i)*KAPPA(g)); 
                                -(G(i) * KAPPA(g)) * sin(KAPPA(g) * H_layers(i)), cos(KAPPA(g) * H_layers(i))];
            
            % Multiply the transfer matrix of this layer to Pimatrix
            Pimatrix = Pimatrix * DeltaInv;  % No pinv or inv
        end
        
        % Calculate the transfer function scalar for non-zero frequencies
        TRANSFERF_SCALAR(g) = Pimatrix(1,1) - (Pimatrix(1,2) * Pimatrix(2,1) / Pimatrix(2,2));
        TRANSFER_ABSOLUTE_ML(g) = abs(TRANSFERF_SCALAR(g));
    end
end
 
% Time Domain Response
for w=1:length(f)
    A_M_i(w) =TRANSFERF_SCALAR(w) * BaseFFT(w);  
    A_EXP_i(w) =TRANSFERF(w) * BaseFFT(w);
    A_VS30_i(w) =TRANSFER_VS30(w) * BaseFFT(w);
end

    A_M_0=real(2*ifft(A_M_i));
    A_EXP_0=real(2*ifft(A_EXP_i));
    A_VS30_0=real(2*ifft(A_VS30_i));
 
% Saving

save IBRH13ChiChi  A_M_0  A_EXP_0  A_VS30_0  BaseMotion t




max(TRANSFER_ABSOLUTE_ML)
max(TRANSFER_ABSOLUTE_EXP)
max(TRANSFER_VS30_ABSOLUTE)

toc;