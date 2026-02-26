clear 
close all
clc

% Discretized Shear Wave Velocity Profile
V0 = 213;  % m/s (Vs at the surface)
alpha = 0.9;  
lam = 0.02;  % Lambda
H = 35;  % m (Total depth)

% Setting Vinf
V_inf = V0 / (1 - alpha);  % m/s
V_BASE = V_inf * (1 - alpha * exp(-lam * H));

% depth array and corresponding Vs values
depths = [];
Vs = [];

% depth interval for each layer
depth_interval = 3;

% Start at the surface (depth = 0)
current_depth = 0;

while current_depth < H
    % Define the depth for the next layer
    layer_depth = current_depth + depth_interval;  
    
    if layer_depth > H
        layer_depth = H;  
    end
    
    % Calculate Vs 
    Vs_upper = V_inf * (1 - alpha * exp(-lam * current_depth));
    Vs_lower = V_inf * (1 - alpha * exp(-lam * layer_depth));
    
    % Average Vs
    Vs_layer = (Vs_upper + Vs_lower) / 2;
    
    % Add to the arrays
    depths = [depths, current_depth, layer_depth];
    Vs = [Vs, Vs_layer, Vs_layer];  % Same Vs for the entire layer (average value)
    
    % Update the current depth
    current_depth = layer_depth;
end

% calculate the exponential function values
z_values = linspace(0, H, 100);  
Vs_function = V_inf * (1 - alpha * exp(-lam * z_values));

% Normalize the x-axis by Vs/V_inf
Vs_normalized = Vs / V_BASE;
Vs_function_normalized = Vs_function / V_BASE;

% Normalize the y-axis by z/H
depths_normalized = depths / H;
z_values_normalized = z_values / H;

% Multi-Layer Transfer Function (Haskel Thompson)
% Soil and wave parameters
rho = 2; % Mg/m^3
xi = 0.05;
C_m = [sqrt(0.8036)];
V0_complex = V0 * sqrt(1 + (2 * 1i * xi)); 
Vinf = V0 / (1 - alpha);
Vinf_complex = Vinf * sqrt(1 + (2 * 1i * xi));  

% Frequency range
f = linspace(0, 80, 2E3);
f1n = (Vinf_complex / (4 * H)) * (2 / pi) * lam * H * C_m(1);
f_normalized = f ./ f1n;

% Transfer function calculation
TRANSFERF = zeros(1, length(f));

z1 = -((1 - alpha) / alpha);
z2 = -((exp(lam * H) - alpha) / alpha);

% Exponential exact solution 
parfor j = 1:length(f)
        k_infty = (2 * pi * f(j)) / Vinf_complex;
        psi = (k_infty) / (lam);
        sq = sqrt(1 - 4 * psi^2);
        rt = 2 * 1i * psi;
        
        S1 = 0.5 * (1 + sq);
        S2 = 0.5 * (1 - sq);

        a1 = (1 - rt + sq) / 2;
        b1 = (1 + rt + sq) / 2;
        c1 = 1 + sq;
        a2 = (1 - rt - sq) / 2;
        b2 = (1 + rt - sq) / 2;
        c2 = 1 - sq;
        a3 = a1 + 1;
        b3 = b1 + 1;
        c3 = c1 + 1;
        a4 = a2 + 1;
        b4 = b2 + 1;
        c4 = c2 + 1;

        F1base = hypergeom([a1, b1], c1, z2);
        F1surf = hypergeom([a1, b1], c1, z1);
        F2base = hypergeom([a2, b2], c2, z2);
        F2surf = hypergeom([a2, b2], c2, z1);

        PH1base = ((-lam * exp(lam * H)) / (2 * alpha)) * hypergeom([a3, b3], c3, z2);
        PH1surf = ((-lam) / (2 * alpha)) * hypergeom([a3, b3], c3, z1);
        PH2base = ((-lam * exp(lam * H)) / (2 * alpha)) * hypergeom([a4, b4], c4, z2);
        PH2surf = (-lam / (2 * alpha)) * hypergeom([a4, b4], c4, z1);
   
        
        E_SURF = lam * 1 * S1 * ...
            (1 - alpha)^(-S2) * F1surf + ...
            (1 - alpha)^(S1) * PH1surf;

        E_BASE = lam * exp(lam * H) * S1 * ...
            (exp(lam * H) - alpha)^(-S2) * F1base + ...
            (exp(lam * H) - alpha)^(S1) * PH1base;
        
        M_SURF = lam * 1 * (S2) * ...
            (1 - alpha)^(-S1) * F2surf + ...
            (1 - alpha)^(S2) * PH2surf;

        M_BASE = lam * exp(lam * H) * (S2) * ...
            (exp(lam * H) - alpha)^(-S1) * F2base + ...
            (exp(lam * H) - alpha)^(S2) * PH2base;

        CONSTANT = (1 - alpha)^(sq) * (F1surf/F2surf);

        TRANSFERF(j) = (E_SURF - CONSTANT*M_SURF) / ...
        (E_BASE - CONSTANT*M_BASE);
end

TRANSFER_ABSOLUTE_EXP = abs(TRANSFERF);

% Multi-Layer Calculation (Haskel Thompson)
H_layers = diff([0, depths]);

% Calculate alpha for each layer
Vs_complex = zeros(1, length(H_layers));
G = zeros(1, length(H_layers));

for p = 1:length(H_layers)
    Vs_complex(p) = Vs(p) * (1 + (1i * xi));
    G(p) = rho*(Vs_complex(p)^2);
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
        TRANSFERF_SCALAR(g) = 1.0; % Static case
        TRANSFER_ABSOLUTE_ML(g) = 1.0;
    else
        for i = 1:length(H_layers)
            % Compute KAPPA for each layer
            KAPPA(g) = (2 * pi * f(g)) / Vs_complex(i);

            % Define DeltaInv
            DeltaInv = [cos(KAPPA(g) * H_layers(i)), sin(KAPPA(g) * H_layers(i)) / (G(i) * KAPPA(g));
                        - (G(i) * KAPPA(g)) * sin(KAPPA(g) * H_layers(i)), cos(KAPPA(g) * H_layers(i))];

            % Multiply the transfer matrix of this layer to Pimatrix
            Pimatrix = Pimatrix * DeltaInv; 
        end

        % Calculate the transfer function scalar
        TRANSFERF_SCALAR(g) = Pimatrix(1,1) - (Pimatrix(1,2) * Pimatrix(2,1) / Pimatrix(2,2));

        TRANSFER_ABSOLUTE_ML(g) = abs(TRANSFERF_SCALAR(g));
    end
end

% --- Plotting ---
figure;

% --- First Subplot: Normalized Shear Wave Velocity Profile
subplot(1, 2, 1); 
plot(Vs_function_normalized, z_values_normalized, 'k--', 'LineWidth', 2);
hold on;
stairs(Vs_normalized, depths_normalized, 'k-', 'LineWidth', 2);
xlabel('V_{s}(z) / V_{s}(H)');
ylabel('z/H');
legend('Theoretical V_{s} profile', 'Numerical discretization', 'Location', 'Best');
ax1 = gca;
ax1.XAxisLocation = 'top';
set(gca, 'YDir', 'reverse');


% --- Second Subplot: Transfer Function Comparison
subplot(1, 2, 2);
plot(f_normalized, TRANSFER_ABSOLUTE_EXP, 'k--', 'LineWidth', 2);  % Exponential Solution
hold on;
plot(f_normalized, TRANSFER_ABSOLUTE_ML, 'k-', 'LineWidth', 2);  % Multi-Layer Solution
xlabel('f/f_{1}');
ylabel('|F(\omega)|');
legend('Exact solution (Eq. 30)', 'Numerical multi-layer solution', 'Location', 'Best');
xlim([0 6]);

