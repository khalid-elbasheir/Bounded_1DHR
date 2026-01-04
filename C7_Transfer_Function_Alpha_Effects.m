clear; close all; clc; tic;

% Parameters
H = 50; 
V0 = 100;
linew = 2.0;

% Lambda*H Parameters
lambdaH_list = [1, 3, 5, 10];      

% Alpha vector
alpha_h = 1e-5;
alpha_neg_far = [-1 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1];
alpha_pos_far = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
alpha_vec = [alpha_neg_far, alpha_h, alpha_pos_far];

% Damping ratios
zeta_list = [0.01 0.05 0.1 0.2];
colors = lines(numel(zeta_list));


Cm_table = [
    5.5030	5.1506	4.8091	4.4785	4.1588	3.8499	3.5519	3.2647	2.9883	2.7225	2.4674	2.2228	1.9886	1.764587429	1.550487716	1.345865714	1.150003915	0.961596203	0.777864717	0.590853178;
    0.3531 0.3453 0.3376 0.3298 0.3220 0.3142 0.3063 0.2984 0.2904 0.2823 0.2742	0.2659	0.2574	0.248832944	0.239975445	0.230798271	0.221182508	0.210929347	0.199655805	0.186395839;
    0.108200279146990	0.107330163484620	0.106446284023730	0.105547489778634	0.104632466539075	0.103699703022280	0.102747445727425	0.101773641107436	0.1008000000000000	0.0997511749525153	0.098700	0.0976000000000000	0.0964756612661318	0.0952976052736806	0.0940621870212438	0.0927558166879089	0.0913583486255917	0.0898373578345756	0.0881334183045796	0.0861074678174316;
    0.0250532083645810	0.0250196566645292	0.0249853312400655	0.0249501754100588	0.0249141250335655	0.0248771068858086	0.0248390371614201	0.0247998180674474	0.0247593355071023	0.0247174540803098	0.024674	0.0246288074656741	0.0245815957984251	0.0245320600015119	0.0244797846699243	0.0244242015418849	0.0243644923314265	0.0242993846991781	0.0242266516791673	0.0241414209387241
];

% FIGURE
fig = figure('Color','w','Units','normalized','Position',[0.03 0.06 0.94 0.84]);

ax2fig = @(ax,x,y) deal( ...
    ax.Position(1) + (x-ax.XLim(1))/(ax.XLim(2)-ax.XLim(1))*ax.Position(3), ...
    ax.Position(2) + (y-ax.YLim(1))/(ax.YLim(2)-ax.YLim(1))*ax.Position(4) );

for k = 1:numel(lambdaH_list)
    lambdaH = lambdaH_list(k);
    lambda  = lambdaH / H;

    % Common alpha grid used for smoothing and for inset masking
    alpha_smooth = linspace(min(alpha_vec), max(alpha_vec), 400);


    % To compute the inset y-range across all zeta curves, store the smoothed arrays
    AF1_smooth_all = zeros(numel(zeta_list), numel(alpha_smooth));

    % === MAIN SUBPLOT ===
    axMain = subplot(2,2,k); hold(axMain,'on'); grid(axMain,'on');

    for zidx = 1:numel(zeta_list)
        zeta = zeta_list(zidx);
        AF1_norm_vals = zeros(size(alpha_vec));

       
        for i = 1:numel(alpha_vec)
            alpha = alpha_vec(i);

            Vinf        = V0/(1 - alpha);
            Vinf_Damped = Vinf * sqrt(1 + 2i*zeta); 

            Cm = Cm_table(k, i);

            % === Exact psi and hypergeometric build ===
            psi = sqrt(Cm) * sqrt(1 + (2*1i*zeta));

            z1 = -((1 - alpha) / alpha);
            z2 = -((exp(lambda*H) - alpha) / alpha);

            sq = sqrt(1 - 4 * psi^2);
            rt = 2 * 1i * psi;

            S1 = 0.5 * (1 + sq);
            S2 = 0.5 * (1 - sq);

            a1 = (1 - rt + sq)/2;  b1 = (1 + rt + sq)/2;  c1 = 1 + sq;
            a2 = (1 - rt - sq)/2;  b2 = (1 + rt - sq)/2;  c2 = 1 - sq;
            a3 = a1 + 1;           b3 = b1 + 1;           c3 = c1 + 1;
            a4 = a2 + 1;           b4 = b2 + 1;           c4 = c2 + 1;

            F1base = hypergeom([a1, b1], c1, z2);
            F1surf = hypergeom([a1, b1], c1, z1);
            F2base = hypergeom([a2, b2], c2, z2);
            F2surf = hypergeom([a2, b2], c2, z1);

            PH1base = ((-lambda * exp(lambda * H)) / (2 * alpha)) * hypergeom([a3, b3], c3, z2);
            PH1surf = ((-lambda) / (2 * alpha))                   * hypergeom([a3, b3], c3, z1);
            PH2base = ((-lambda * exp(lambda * H)) / (2 * alpha)) * hypergeom([a4, b4], c4, z2);
            PH2surf = (-lambda / (2 * alpha))                      * hypergeom([a4, b4], c4, z1);

            E_SURF = lambda * S1 * (1 - alpha)^(-S2) * F1surf + (1 - alpha)^(S1) * PH1surf;
            E_BASE = lambda * exp(lambda * H) * S1 * (exp(lambda * H) - alpha)^(-S2) * F1base + ...
                     (exp(lambda * H) - alpha)^(S1) * PH1base;

            M_SURF = lambda * S2 * (1 - alpha)^(-S1) * F2surf + (1 - alpha)^(S2) * PH2surf;
            M_BASE = lambda * exp(lambda * H) * S2 * (exp(lambda * H) - alpha)^(-S1) * F2base + ...
                     (exp(lambda * H) - alpha)^(S2) * PH2base;

            CONSTANT = (1 - alpha)^(sq) * (F1surf / F2surf);
            Fsurf = (E_SURF - CONSTANT*M_SURF) / (E_BASE - CONSTANT*M_BASE);

            AF1_norm_vals(i) = abs(Fsurf) / (2/(pi*zeta));
        end

        
        AF1_smooth = spline(alpha_vec, AF1_norm_vals, alpha_smooth);
        AF1_smooth_all(zidx, :) = AF1_smooth;  % stash for inset y-limits

        
        plot(axMain, alpha_smooth, AF1_smooth, 'LineWidth', linew, ...
             'Color', colors(zidx,:), 'DisplayName', sprintf('\\zeta = %g%%', zeta*100));

        

    end
      
    xl = [-1 1]; xlim(axMain, xl);
    xlabel(axMain, '\alpha'); ylabel(axMain, '|F(\omega)| / (\pi\xi/2)');
    title(axMain, sprintf('\\lambdaH = %.3g', lambdaH), 'Interpreter','tex', 'FontWeight','bold');
    set(axMain, 'FontSize', 12);
    legend(axMain, 'show', 'Location', 'best');
    
end

toc;
