clear
close all
clc
tic;
 

% Core region
Argument_dense = linspace(-100, 100, 3E4);
Argument_dense = Argument_dense(Argument_dense ~= 0 & Argument_dense ~= 1);  

% Tail regions (log spaced)
neg_tail = -logspace(log10(1e2), log10(1e4), 1000);  
pos_tail =  logspace(log10(1e2), log10(1e4), 1000);  

% Key points
key_points = [-1, 0, 0.99999];  % exact points

% Combine and sort
Argument = unique([neg_tail, Argument_dense, pos_tail, key_points]);
LENGTH = length(Argument);

psi_values = [20 10 5 1 0.5];
num_psi = length(psi_values);


% Preallocate arrays
    PHI2 = zeros(num_psi, LENGTH);

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
            
            % Compute Hypergeometric functions
            alpha = 0.5;
            lambda = 0.02;
            exponent = alpha*(1-Argument(o));
            Multiplier = (- lambda*exponent)/(2*alpha);
            PHI = arrayfun(@(a, b, c, z) Multiplier*...
                            hypergeom([a + 1, b + 1], c + 1, z), ...
                            a_values, b_values, c_values, [Argument(o), Argument(o)]);

            PHI2(p, o) = PHI(2);
        end
    end


    % Create new masks for 3 specific regions
mask_neg = Argument < 0;                                 
mask_pos_small = Argument >= -1 & Argument <= 1;            
mask_pos_large = Argument > 0 & Argument <= 1e12;



% Create figure
figure('Units', 'normalized', 'Position', [0.05, 0.3, 0.9, 0.4]); 

% ==== Plot 1: Negative x ====
subplot(1,3,1);
hold on;
legend_entries = cell(1, num_psi);
for p = 1:num_psi
    plot(1./Argument(mask_neg), PHI2(p, mask_neg), 'LineWidth', 1.5);
    legend_entries{p} = ['\psi = ', num2str(real(psi_values(p)), '%.2f')];
end
xlabel('1/x (Negative x)');
ylabel('P_2');
legend(legend_entries, 'Location', 'best');
set(gca, 'XDir', 'reverse');
set(gca, 'YAxisLocation', 'right');
xlim([-1 0]);
ylim([-0.06 0.05])
xticks([-1 -0.8 -0.6 -0.4 -0.2 0]);
grid on;
hold off;


% ==== Plot 2: Small positive x ====
subplot(1,3,2);
hold on;
legend_entries = cell(1, num_psi);
for p = 1:num_psi
    plot(Argument(mask_pos_small), PHI2(p, mask_pos_small), 'LineWidth', 1.5);
    legend_entries{p} = ['\psi = ', num2str(real(psi_values(p)), '%.2f')];
end
xlabel('x');
ylabel('P_2');
legend(legend_entries, 'Location', 'best');
xlim([-1 1]);
ylim([-0.06 0.05])
grid on;
hold off;

% ==== Plot 3: Large positive x ====
subplot(1,3,3);
hold on;
legend_entries = cell(1, num_psi);
for p = 1:num_psi
    plot(1./Argument(mask_pos_large), PHI2(p, mask_pos_large), 'LineWidth', 1.5);
    legend_entries{p} = ['\psi = ', num2str(real(psi_values(p)), '%.2f')];
end
xlabel('1/x (Large Positive)');
ylabel('P_2');
legend(legend_entries, 'Location', 'best');
set(gca, 'XDir', 'reverse');              
xlim([0 1]);
ylim([-0.06 0.05])
grid on;
hold off;

toc;