clear all
close all
clc
tic;

% Parameters
xi = 0.0;  % Damping factor
lambda = 0.02;  % Wavelength factor

alpha_values = [0.5];  % Different alpha values to test 

% Core region
ZETA12_dense = linspace(-100, 100, 30000);
ZETA12_dense = ZETA12_dense(ZETA12_dense ~= 0 & ZETA12_dense ~= 1);  % avoid duplication

% Tail regions (log spaced)
neg_tail = -logspace(log10(1e2), log10(1e4), 1000);  % from -100 to -1e4
pos_tail =  logspace(log10(1e2), log10(1e4), 1000);  % from 100 to 1e4

% Key points to explicitly include
key_points = [0, 0.99999];  % exact points

% Combine and sort
ZETA12 = unique([neg_tail, ZETA12_dense, pos_tail, key_points]);
LENGTH = length(ZETA12);



psi_values = [10 5 1 0.5];%(2 * pi * f) ./ (lambda * Vinf_complex);
num_psi = length(psi_values);

% Loop over different alpha values
for a_idx = 1:length(alpha_values)
    alpha = alpha_values(a_idx);
    
        % Preallocate arrays
    %F1 = zeros(num_psi, LENGTH);
    %F2 = zeros(num_psi, LENGTH);
    PHI1 = zeros(num_psi, LENGTH);
    %PHI2 = zeros(num_psi, LENGTH);

    parfor p = 1:num_psi
        psi = psi_values(p);

        % Hypergeometric parameters
        a1 = (1 - 2 * 1i * psi + sqrt(1 - 4 * psi^2)) / 2;
        b1 = (1 + 2 * 1i * psi + sqrt(1 - 4 * psi^2)) / 2;
        c1 = 1 + sqrt(1 - 4 * psi^2);

        a_values = [a1, a1 + 1 - c1];
        b_values = [b1, b1 + 1 - c1];
        c_values = [c1, 2 - c1];

        for o = 1:LENGTH
            EXPlambdaZ = alpha * (1 - ZETA12(o));

            % Compute Hypergeometric functions
            %F = arrayfun(@(a, b, c, z) hypergeom([a, b], c, z), ...
             %            a_values, b_values, c_values, [ZETA12(o), ZETA12(o)]);
            PHI = arrayfun(@(a, b, c, z) ((-lambda * EXPlambdaZ) / (2 * alpha)) * ...
                            hypergeom([a + 1, b + 1], c + 1, z), ...
                            a_values, b_values, c_values, [ZETA12(o), ZETA12(o)]);

            %F1(p, o) = F(1);
            %F2(p, o) = F(2);
            PHI1(p, o) = PHI(1);
            %PHI2(p, o) = PHI(2);
        end
    end


    % Create new masks for 3 specific regions
mask_neg = ZETA12 < 0;                                 % Plot 1
mask_pos_small = ZETA12 >= -1 & ZETA12 <= 1;            % Plot 2: INCLUDE ζ = 0 and 1
mask_pos_large = ZETA12 > 0 & ZETA12 <= 1e12;



% Create a figure with horizontal layout (1 row, 3 columns)
figure('Units', 'normalized', 'Position', [0.05, 0.3, 0.9, 0.4]); % Wider layout

% ==== Plot 1: Negative zeta (1/zeta -> large negative) ====
subplot(1,3,1);
hold on;
legend_entries = cell(1, num_psi);
for p = 1:num_psi
    plot(1./ZETA12(mask_neg), PHI1(p, mask_neg), 'LineWidth', 1.5);
    legend_entries{p} = ['\psi = ', num2str(real(psi_values(p)), '%.2f')];
end
xlabel('1/\zeta (Negative \zeta)');
ylabel('F_1');
title('Plot 1: \zeta < 0');
legend(legend_entries, 'Location', 'best');
set(gca, 'XDir', 'reverse');
set(gca, 'YAxisLocation', 'right');
xlim([-1 0]);
xticks([-1 -0.8 -0.6 -0.4 -0.2 0]);
grid on;
hold off;


% ==== Plot 2: Small positive zeta (0 < \zeta <= 1) ====
subplot(1,3,2);
hold on;
legend_entries = cell(1, num_psi);
for p = 1:num_psi
    plot(ZETA12(mask_pos_small), PHI1(p, mask_pos_small), 'LineWidth', 1.5);
    legend_entries{p} = ['\psi = ', num2str(real(psi_values(p)), '%.2f')];
end
xlabel('\zeta (Small Positive)');
ylabel('F_1');
title('Plot 2: 0 < \zeta \leq 1');
legend(legend_entries, 'Location', 'best');
xlim([-1 1]);
grid on;
hold off;

% ==== Plot 3: Large positive zeta (1/zeta from 0.1 to 1) ====
subplot(1,3,3);
hold on;
legend_entries = cell(1, num_psi);
for p = 1:num_psi
    plot(1./ZETA12(mask_pos_large), PHI1(p, mask_pos_large), 'LineWidth', 1.5);
    legend_entries{p} = ['\psi = ', num2str(real(psi_values(p)), '%.2f')];
end
xlabel('1/\zeta (Large Positive)');
ylabel('F_1');
title('Plot 3: 1 \leq \zeta \leq 10');
legend(legend_entries, 'Location', 'best');
set(gca, 'XDir', 'reverse');              % From 1 to 0
xlim([0 1]);                              % Set desired x limits
grid on;
hold off;



end

toc;
