clear all
close all
clc
c_light = 3e8;
colors = distinguishable_colors(15);
iteration = 1e6;
[d_u_min, d_u_max, d_r_min, d_r_max, N_F, B, P_0, F_c, ...
    I_min, I_max, eta, V_t, I_0, Phi, Theta, q_e, P_LED, I_i, ...
    A_R, h_Delta, B_RF, B_VLC, ~] = initialization();
N_0 = P_0 + 10*log10(B_RF) + N_F ;
N_0 = 10^(N_0/10);
l = h_Delta;
g = 1;
Phi_1_2 = deg2rad(60);
gamma = -1/log2(cos(Phi_1_2));
sigma_2 = q_e*I_i*B_VLC;
beta = 1.6;
theta_2_total_deg = [5, 20];
theta_1_total = deg2rad([10,30]);
theta_2_total = deg2rad(theta_2_total_deg);


radious_total = 0:0.1:8;

%% Define Arrays
pr_out_vlc = zeros(length(theta_2_total), length(radious_total));
pr_out_rf = zeros(length(theta_2_total), length(radious_total));
pr_out_vlc_rf = zeros(length(theta_2_total), length(radious_total));
E_h_sim= zeros(length(theta_2_total), length(radious_total));
I_b = 0.85;
R_th = 1e6;
for cnt_theta = 1: length(theta_2_total)
    theta_2 = theta_2_total(cnt_theta);
    theta_1 = theta_1_total(cnt_theta);

    for cnt_r = 1:length(radious_total)

        radious = radious_total(cnt_r);
        d_u_max = radious;
        d_r_max = radious;
        mu_r = radious;

        A_VLC = I_max - I_b;

        [x1, y1, x2, y2] = randomPointsInCircle_sqrt(d_u_max, iteration);
        d_ur = sqrt((x1-x2).^2 + (y1-y2).^2);
        d_u = raylrnd(mu_r, [1, iteration]);
        d_r = raylrnd(mu_r, [1, iteration]);

        thetha_random_u =  normrnd(theta_1, theta_2, [1, iteration]);

        thetha_random_r =  normrnd(theta_1, theta_2, [1, iteration]);


        B = (exp(1)*(eta*P_LED*A_VLC)^2) / (2*pi*sigma_2);

        h_c_u = (gamma + 1).*A_R.*(h_Delta^gamma).*g.*((h_Delta^2 + d_u.^2).^(-(gamma+2)/2))./(2*pi);

        H_VLC_u = h_c_u.*cos(thetha_random_u);
        R_VLC_u = B_VLC.*log2(1+ B.*H_VLC_u.^2);


        h_c_r = (gamma + 1)*A_R*(h_Delta^gamma)*g*((h_Delta^2 + d_r.^2).^(-(gamma+2)/2))/(2*pi);
        H_VLC_r = h_c_r.*cos(thetha_random_r);
        R_VLC_r = B_VLC.*log2(1+ B.*H_VLC_r.^2);

        G_RF = (4*pi*F_c/c_light)^2.*d_ur.^beta;
        h_RF = h_RF_func(1, 1, a); % Rayleigh fading
        zeta = abs(h_RF).^2 ./ (G_RF.*N_0);
        Beta = eta.*P_LED.*H_VLC_r;
        f =  0.75.*eta.*H_VLC_r.*P_LED.*V_t;
        E_h = (f.*I_b.*log(1+Beta.*I_b./I_0));
        P_h = (E_h);
        R_RF = B_RF .* log2(1+P_h.*zeta);


        pr_out_vlc(cnt_theta, cnt_r) = sum(R_VLC_u < R_th)/iteration;
        pr_out_rf(cnt_theta, cnt_r) = sum(R_RF < R_th)/iteration;
        pr_out_vlc_rf(cnt_theta, cnt_r) = sum(R_VLC_u < R_th)/iteration * (sum(R_VLC_r < R_th)/iteration + sum(R_RF(R_VLC_r >= R_th) < R_th)/iteration);
        E_h_sim(cnt_theta, cnt_r) = mean(E_h);

    end

end

% Plot the results
figure
plot(radious_total, pr_out_vlc_rf(1,:), '->', 'Color', colors(1,:), 'MarkerIndices',1:10:length(radious_total))
hold on
grid on
plot(radious_total, pr_out_vlc_rf(2,:), '-o', 'Color', colors(7,:), 'MarkerIndices',1:10:length(radious_total))
plot(radious_total, pr_out_vlc(1,:), '-.<', 'Color', colors(2,:), 'MarkerIndices',1:10:length(radious_total))
plot(radious_total, pr_out_vlc(2,:), '-.+', 'Color', colors(4,:), 'MarkerIndices',1:10:length(radious_total))

legend({ '$P^{\mathcal{O}}_{\rm H-VR}$ (simulation-$\theta \sim \mathcal{N}(10^{\circ}, 5^{\circ})$)', ...
    '$P^{\mathcal{O}}_{\rm H-VR}$ (simulation-$\theta \sim \mathcal{N}(30^{\circ}, 20^{\circ})$)', ...
    '$P_{\rm{direct}}^{\mathcal{O}}$ (analytic-$\theta \sim \mathcal{N}(10^{\circ}, 5^{\circ})$)', ...
    '$P_{\rm{direct}}^{\mathcal{O}}$ (analytic-$\theta \sim \mathcal{N}(30^{\circ}, 20^{\circ})$)'},'Interpreter','latex')
xlabel({'Rayleigh scale parameter $\sigma_{\rm{R}}$'},'Interpreter','latex')
ylabel('Outage probability')
ylim([0 1])
ax=gca;
ax.GridAlpha=1;
ax.GridLineStyle=':';
ax.MinorGridAlpha=1;
pause(0.1)
set(gca,'FontSize',12)
set(gca, 'FontName', 'Times New Roman')
set(findall(gca, 'Type', 'Line'),'LineWidth',3)
set(gcf,'windowstyle','normal')
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
