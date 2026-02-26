clear
clc 
close all
tic;

% Fixed parameters
H = 50;
alpha = 1;

% lambdaH values and corresponding psi sets (only use first freq psi(1)) 
lambdaH_list = [5, 8, 10, 15];
psi_values_list = {
    [sqrt(0.0829), sqrt(0.3367), sqrt(0.5646)];
    [sqrt(0.0366), sqrt(0.2402), sqrt(0.9333)];
    [sqrt(0.0240), sqrt(0.1806), sqrt(0.9333)];
    [sqrt(0.0109), sqrt(0.2271), sqrt(0.9151)];
};

% Discretize depth
z = 1E-5:0.2:H;
z_H = z / H;
nZ = length(z);

% Prepare outputs
STRAIN_RATIOS = zeros(nZ, length(lambdaH_list));
DISPLACEMENT_RATIOS = zeros(nZ, length(lambdaH_list));

% Loop over lambdaH values
for idx = 1:length(lambdaH_list)
    lambdaH = lambdaH_list(idx);
    psi = psi_values_list{idx}(1);  % Only first frequency

    % Precompute constants
    S1 = 0.5 * (1 + sqrt(1 - 4 * psi^2));
    S2 = 0.5 * (1 - sqrt(1 - 4 * psi^2));
    a1 = (1 - 2i*psi + sqrt(1 - 4*psi^2)) / 2;
    b1 = (1 + 2i*psi + sqrt(1 - 4*psi^2)) / 2;
    c1 = 1 + sqrt(1 - 4*psi^2);
    a2 = (1 - 2i*psi - sqrt(1 - 4*psi^2)) / 2;
    b2 = (1 + 2i*psi - sqrt(1 - 4*psi^2)) / 2;
    c2 = 1 - sqrt(1 - 4*psi^2);
    a3 = a1 + 1;
    b3 = b1 + 1;
    c3 = c1 + 1;
    a4 = a2 + 1;
    b4 = b2 + 1;
    c4 = c2 + 1;

    Z2 = -((exp(lambdaH) - alpha) / alpha);
    Z3 = -((1 - alpha) / alpha);

    
    % compute 
    F1SURF = hypergeom([a1, b1], c1, Z3);
    F1Base = hypergeom([a1,b1], c1, Z2);

   
    PH1surf = (-(lambdaH/H)) / (2 * alpha) * hypergeom([a3, b3], c3, Z3);
    PH1Base = ((-(lambdaH/H) * exp(lambdaH)) / (2 * alpha)) * hypergeom([a3,b3], c3, Z2);

    
    E_surf = (lambdaH/H) * S1 * (1 - alpha)^(-S2) * F1SURF + (1 - alpha)^(S1) * PH1surf;
    E_base = (lambdaH/H) * exp(lambdaH) * S1 * (exp(lambdaH) - alpha)^(-S2) * F1Base + (exp(lambdaH) - alpha)^(S1) * PH1Base;

    
    CON1 = 1;
    u0 = CON1 * ( E_surf );
    uH = CON1 * ( E_base );

    Denom_strain = - CON1 * psi^(2) * (lambdaH/H)^(2) * exp(2*lambdaH) * ( (exp(lambdaH) - alpha)^(S1 - 2) * F1Base);

    for i = 1:nZ
        z_i = z(i);
        Z1 = -((exp(lambdaH*z_i/H) - alpha) / alpha);

        % --- Strain Ratio ---
        F1Depth = hypergeom([a1, b1], c1, Z1);
        
        Num_strain = -CON1 * psi^(2) * (lambdaH/H)^(2) * exp(2*lambdaH*z_i/H) * ...
                     ( (exp(lambdaH*z_i/H) - alpha)^(S1 - 2) * F1Depth) ;
        
        STRAIN_RATIOS(i, idx) = Num_strain / Denom_strain;

        % --- Displacement Ratio ---
        PH1Depth = (-(lambdaH/H) * exp(lambdaH*z_i/H)) / (2 * alpha) * hypergeom([a3, b3], c3, Z1);
        
        E_depth = (lambdaH/H) * exp(lambdaH*z_i/H) * S1 * (exp(lambdaH*z_i/H) - alpha)^(-S2) * F1Depth + ...
                  (exp(lambdaH*z_i/H) - alpha)^S1 * PH1Depth;
        
        u_depth = CON1 * ( E_depth );

        DISPLACEMENT_RATIOS(i, idx) = u_depth / uH;
    end
end

% ==================== Plotting ====================
figure;



% Subplot 1: Strain
subplot(1,2,1);
hold on; 
for idx = 1:length(lambdaH_list)
    plot(STRAIN_RATIOS(:, idx), z_H);
end
xlabel('\gamma(z) / \gamma(H)');
ylabel('z/H');
xlim([0 6]);
set(gca, 'YDir', 'reverse');
set(gca, 'XAxisLocation', 'top');
title('Strain Profile');
legend({'\lambdaH = 5', '\lambdaH = 8', '\lambdaH = 10', '\lambdaH = 15'}, 'Location', 'best');

% Subplot 2: Displacement
subplot(1,2,2);
hold on; grid on;
for idx = 1:length(lambdaH_list)
    plot(DISPLACEMENT_RATIOS(:, idx), z_H);
end
xlabel('u(z) / u(H)');
ylabel('z/H');
set(gca, 'YDir', 'reverse');
set(gca, 'XAxisLocation', 'top');
title('Displacement Profile');
legend({'\lambdaH = 5', '\lambdaH = 8', '\lambdaH = 10', '\lambdaH = 15'}, 'Location', 'best');

toc;
