clear; close all; clc; tic;

% Model Parameters
H = 50;
V0 = 100;
HOMOG_REF_TF = 12.73;
linew = 2.0;

% alpha values
alpha_set   = [0.5, 0.9, -0.5, -0.9];
lambdaH_vec = [0.01 0.51 0.71 1:1:10];
n_lambda    = numel(lambdaH_vec);
n_alpha     = numel(alpha_set);

% Damping ratios
xi_list = [0.01 0.05 0.1 0.2];
colors    = lines(numel(xi_list));  
nz        = numel(xi_list);


Vinf_vec = V0 ./ (1 - alpha_set);


Cm_table = [
4.03E+03	3.9388	2.3091	1.3459	0.4518	0.2308	0.1394	0.0928	0.0658	0.049	0.0378	0.03	0.0244;
279 1.0574 0.7939  0.590850 0.303840 0.186400 0.12333 0.086107 0.062797 0.047488 0.037011 0.029579 0.024141;
5.22E+04	17.278	8.3249	3.8499	0.7861	0.3142	0.1671	0.1037	0.0708	0.0515	0.0391	0.0308	0.0249;
88483 25.1338 11.6845 5.1506   0.9282   0.3453   0.1767  0.107330 0.072361 0.052235 0.039557 0.031035 0.025020
];

figure('Color','w','Units','normalized','Position',[0.03 0.06 0.94 0.84]);


ax2fig = @(ax,x,y) deal( ...
    ax.Position(1) + (x-ax.XLim(1))/(ax.XLim(2)-ax.XLim(1))*ax.Position(3), ...
    ax.Position(2) + (y-ax.YLim(1))/(ax.YLim(2)-ax.YLim(1))*ax.Position(4) );

for aidx = 1:n_alpha
    alpha = alpha_set(aidx);
    Vinf  = Vinf_vec(aidx);

    ax = subplot(2,2,aidx); hold(ax,'on'); 

    
    lambdaH_smooth = linspace(min(lambdaH_vec), max(lambdaH_vec), 200);
    AF1_smooth_all = zeros(nz, numel(lambdaH_smooth));

    for zidx = 1:nz
        xi = xi_list(zidx);
        Vinf_Damped = Vinf * sqrt(1 + 2i*xi); 

        AF1_norm_vals = zeros(1, n_lambda);

        for i = 1:n_lambda
            lambdaH = lambdaH_vec(i);
            lambda  = lambdaH / H;
            Cm      = Cm_table(aidx, i);

            % ----- analytical solution -----
            psi = sqrt(Cm)* sqrt(1 + (2*1i*xi));
            z1 = -((1 - alpha) / alpha);
            z2 = -((exp(lambda*H) - alpha) / alpha);

            sq = sqrt(1 - 4 * psi^2);
            rt = 2 * 1i * psi;

            S1 = 0.5 * (1 + sq);
            S2 = 0.5 * (1 - sq);

            a1 = (1 - rt + sq) / 2;  b1 = (1 + rt + sq) / 2;  c1 = 1 + sq;
            a2 = (1 - rt - sq) / 2;  b2 = (1 + rt - sq) / 2;  c2 = 1 - sq;
            a3 = a1 + 1;  b3 = b1 + 1;  c3 = c1 + 1;
            a4 = a2 + 1;  b4 = b2 + 1;  c4 = c2 + 1;

            F1base = hypergeom([a1, b1], c1, z2);
            F1surf = hypergeom([a1, b1], c1, z1);
            F2base = hypergeom([a2, b2], c2, z2);
            F2surf = hypergeom([a2, b2], c2, z1);

            PH1base = ((-lambda * exp(lambda * H)) / (2 * alpha)) * hypergeom([a3, b3], c3, z2);
            PH1surf = ((-lambda) / (2 * alpha)) * hypergeom([a3, b3], c3, z1);
            PH2base = ((-lambda * exp(lambda * H)) / (2 * alpha)) * hypergeom([a4, b4], c4, z2);
            PH2surf = (-lambda / (2 * alpha)) * hypergeom([a4, b4], c4, z1);

            E_SURF = lambda * 1 * S1 * (1 - alpha)^(-S2) * F1surf + (1 - alpha)^(S1) * PH1surf;
            E_BASE = lambda * exp(lambda * H) * S1 * (exp(lambda * H) - alpha)^(-S2) * F1base + (exp(lambda * H) - alpha)^(S1) * PH1base;
            M_SURF = lambda * 1 * (S2) * (1 - alpha)^(-S1) * F2surf + (1 - alpha)^(S2) * PH2surf;
            M_BASE = lambda * exp(lambda * H) * (S2) * (exp(lambda * H) - alpha)^(-S1) * F2base + (exp(lambda * H) - alpha)^(S2) * PH2base;

            CONSTANT = (1 - alpha)^(sq) * (F1surf / F2surf);
            Fsurf = (E_SURF - CONSTANT*M_SURF) / (E_BASE - CONSTANT*M_BASE);

            % Normalized amplification
            AF1_norm_vals(i) = abs(Fsurf) / (2/(pi*xi));  
        end
        AF1_norm_vals(1)=1;
        
        AF1_smooth = spline(lambdaH_vec, AF1_norm_vals, lambdaH_smooth);

        AF1_smooth_all(zidx,:) = AF1_smooth;

        % Plot
        plot(lambdaH_vec, AF1_norm_vals, 'ko', 'MarkerFaceColor','k');
        plot(ax, lambdaH_smooth, AF1_smooth, 'Color', colors(zidx,:), ...
            'LineWidth', linew, 'DisplayName', sprintf('\\zeta = %.0f%%', xi*100));
    end

    xlabel(ax,'\lambdaH'); 
    ylabel(ax,'A_{1}/A_{1H}');
    title(ax, sprintf('\\alpha = %.2f', alpha), 'FontWeight','bold');
    xlim(ax, [0 max(lambdaH_vec)]);
    ylim([0.9,1.3])
    set(ax, 'FontSize', 12);
    legend(ax,'show','Location','best');

   %% ----- Build the zoom inset for lambdaH in [1,3] -----
    zoom_x = [1 3];

    % Find y-lims from all curves within the zoom window
    idx_zoom = lambdaH_smooth >= zoom_x(1) & lambdaH_smooth <= zoom_x(2);
    Yzoom_all = AF1_smooth_all(:, idx_zoom);
    y_min = min(Yzoom_all, [], 'all');
    y_max = max(Yzoom_all, [], 'all');
    pad = 0.06*(y_max - y_min + eps);
    zoom_y = [y_min - pad, y_max + pad];

    % Place inset inside the current subplot (top-right quadrant)
    inset_scale = 0.35; % size relative to subplot
    inset_pos   = ax.Position;
    axInset = axes('Position', [inset_pos(1)+0.60*inset_pos(3), ...
                                inset_pos(2)+0.52*inset_pos(4), ...
                                inset_scale*inset_pos(3), ...
                                inset_scale*inset_pos(4)]);
    hold(axInset,'on'); grid(axInset,'on');
    for zidx = 1:nz
        plot(axInset, lambdaH_smooth, AF1_smooth_all(zidx,:), ...
            'Color', colors(zidx,:), 'LineWidth', 1.6);
    end
    xlim(axInset, zoom_x);
    ylim(axInset, zoom_y);
    set(axInset,'FontSize',10,'Box','on');
    
    title(axInset,'Zoom');

    % ----- Draw a red hollow circle over the zoom area on the main plot -----
    % Circle center at x=2 and average y across curves at x=2
    y_at_2 = zeros(nz,1);
    for zidx = 1:nz
        y_at_2(zidx) = interp1(lambdaH_smooth, AF1_smooth_all(zidx,:), 2, 'linear', 'extrap');
    end
    yc = mean(y_at_2, 'omitnan');
    xc = 2;

    % Choose a circle radius as ~12% of the subplot width in normalized coords
    [xcn, ycn] = ax2fig(ax, xc, yc);
    r_norm = 0.06; % radius in figure-normalized units (tweak if you like)
    ellipse_rect = [xcn - r_norm, ycn - r_norm, 2*r_norm, 2*r_norm];

    ann_circle = annotation('ellipse', ellipse_rect, ...
        'Color','r','LineWidth',1.6, 'LineStyle','-','FaceColor','none');

    % ----- Arrow from circle to the inset -----
    inset_center = [axInset.Position(1)+0.5*axInset.Position(3), ...
                    axInset.Position(2)+0.5*axInset.Position(4)];
    % Start point: rightmost point of the circle
    start_pt = [xcn + r_norm, ycn];
    annotation('textarrow', [start_pt(1) inset_center(1)], ...
                           [start_pt(2) inset_center(2)], ...
               'Color','r','LineWidth',1.4, 'HeadLength',8, 'HeadWidth',8);

end

toc;