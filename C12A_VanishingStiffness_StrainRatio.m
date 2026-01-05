clear
close all
clc
tic
% Alpha is Fixed - 4 types of LambdaH - 3 Frequniceis
H = 50;
Vinf = 200;
alpha = 0.9995;


lambdaH_list = [5, 8, 10, 15];
psi_values_list = {
    [sqrt(0.0829), sqrt(0.3367), sqrt(0.5646)];
    [sqrt(0.0366), sqrt(0.2402), sqrt(0.9333)];
    [sqrt(0.024), sqrt(0.1806), sqrt(0.9333)];
    [sqrt(0.0109), sqrt(0.2271), sqrt(0.9151)];
};

z = 0:0.2:H;
z_H = z ./ H;
nZ = length(z);

figure;

for lambdaIdx = 1:length(lambdaH_list)
    lambdaH = lambdaH_list(lambdaIdx);
    psi = psi_values_list{lambdaIdx};
    nModes = length(psi);

    STRAIN_MATRIX = NaN(nZ, nModes);

    for j = 1:nModes
        psi_j = psi(j);

        S1 = 0.5 * (1 + sqrt(1 - 4 * psi_j^2));
        a1 = (1 - 2*1i*psi_j + sqrt(1 - 4*psi_j^2)) / 2;
        b1 = (1 + 2*1i*psi_j + sqrt(1 - 4*psi_j^2)) / 2;
        c1 = 1 + sqrt(1 - 4*psi_j^2);

        Z2 = -((exp(lambdaH) - alpha) / alpha);

        % Base value
        F1Base = hypergeom([a1, b1], c1, Z2);
        Denominator = - psi_j^2 * (lambdaH/H)^2 * exp(2*lambdaH) * (exp(lambdaH) - alpha)^(S1 - 2) * F1Base;

        % Strain ratio profile
        GAMMA_ratio = NaN(1, nZ);
        for i = 1:nZ
            Zi = -((exp(lambdaH * z(i)/H) - alpha) / alpha);
            F1Depth = hypergeom([a1, b1], c1, Zi);

            Numerator = - psi_j^2 * (lambdaH/H)^2 * exp(2 * lambdaH * z(i)/H) * ...
                        (exp(lambdaH * z(i)/H) - alpha)^(S1 - 2) * F1Depth;

            GAMMA_ratio(i) = Numerator / Denominator;
        end

        STRAIN_MATRIX(:, j) = GAMMA_ratio.';
    end

    % Plot subplot
    subplot(2, 2, lambdaIdx);
    hold on;

    for j = 1:nModes
        plot(STRAIN_MATRIX(:, j), z_H);
    end

    xlabel('\gamma(z) / \gamma(H)');
    ylabel('z/H');
    xlim([-8 8]);
    set(gca, 'YDir', 'reverse');
    set(gca, 'XAxisLocation', 'top');
    title(['\lambdaH = ', num2str(lambdaH)], 'FontWeight', 'bold');
    legend({'f₁', 'f₂', 'f₃'}, 'Location', 'best');

end

toc;
