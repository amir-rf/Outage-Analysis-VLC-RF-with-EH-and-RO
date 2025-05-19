clear all
close all
clc


c_light = 3e8; % Speed of light in m/s
colors = distinguishable_colors(15); % Generate distinguishable plot colors
iteration = 1e6; % Number of Monte Carlo iterations

% Load simulation and system parameters
[d_u_min, d_u_max, d_r_min, d_r_max, N_F, B, P_0, F_c, ...
    I_min, I_max, eta, V_t, I_0, Phi, Theta, q_e, P_LED, I_i, ...
    A_R, h_Delta, B_RF, B_VLC, R_th] = initialization();

% Convert noise power from dBm to linear scale
N_0 = P_0 + 10*log10(B_RF) + N_F ;
N_0 = 10^(N_0/10);

% VLC channel and simulation setup
l = h_Delta;
g = 1;
Phi_1_2 = deg2rad(60);
gamma = -1/log2(cos(Phi_1_2));
sigma_2 = q_e*I_i*B_VLC;
beta = 1.6;
theta_1 = deg2rad(10);
theta_2 = deg2rad(40);


radious = 7;
I_b_total = [0.55, 0.95];
threshold_total = 1e6:1e5:1e8;
for cnt_th = 1: length(threshold_total)
    R_th = threshold_total(cnt_th);
    for cnt_i_b = 1:length(I_b_total)

        I_b = I_b_total(cnt_i_b);
        d_u_max = radious;
        d_r_max = radious;

        A_VLC = I_max - I_b;

        [x1, y1, x2, y2] = randomPointsInCircle_sqrt(d_u_max, iteration);
        d_ur = sqrt((x1-x2).^2 + (y1-y2).^2);
        [x1, y1, x2, y2] = randomPointsInCircle(d_u_max, iteration);
        d_u = sqrt(x1.^2 + y1.^2);
        d_r = sqrt(x2.^2 + y2.^2);

        thetha_random_u =  unifrnd(theta_1, theta_2, [1, iteration]);

        thetha_random_r =  unifrnd(theta_1, theta_2, [1, iteration]);


        B = (exp(1)*(eta*P_LED*A_VLC)^2) / (2*pi*sigma_2);

        h_c_u = (gamma + 1).*A_R.*(h_Delta^gamma).*g.*((h_Delta^2 + d_u.^2).^(-(gamma+2)/2))./(2*pi);

        H_VLC_u = h_c_u.*cos(thetha_random_u);
        R_VLC_u = B_VLC.*log2(1+ B.*H_VLC_u.^2);


        h_c_r = (gamma + 1)*A_R*(h_Delta^gamma)*g*((h_Delta^2 + d_r.^2).^(-(gamma+2)/2))/(2*pi);
        H_VLC_r = h_c_r.*cos(thetha_random_r);
        R_VLC_r = B_VLC.*log2(1+ B.*H_VLC_r.^2);

        G_RF = (4*pi*F_c/c_light)^2.*d_ur.^beta;
        h_RF = h_RF_func(1, 1, iteration); % Rayleigh fading
        zeta = abs(h_RF).^2 ./ (G_RF.*N_0);
        Beta = eta.*P_LED.*H_VLC_r;
        f =  0.75.*eta.*H_VLC_r.*P_LED.*V_t;
        E_h = (f.*I_b.*log(1+Beta.*I_b./I_0));
        P_h = (E_h);
        R_RF = B_RF .* log2(1+P_h.*zeta);


        pr_out_vlc(cnt_th, cnt_i_b) = sum(R_VLC_u < R_th)/iteration;
        pr_out_rf(cnt_th, cnt_i_b) = sum(R_RF < R_th)/iteration;
        pr_out_vlc_rf(cnt_th, cnt_i_b) = sum(R_VLC_u < R_th)/iteration * (sum(R_VLC_r < R_th)/iteration + sum(R_RF(R_VLC_r >= R_th) < R_th)/iteration);
        E_h_sim(cnt_th, cnt_i_b) = mean(E_h);

    end

end


% Theoretical analysis follows in the next loop
for cnt_th = 1: length(threshold_total)
    R_th = threshold_total(cnt_th);

    h_c_prime2 = ((gamma + 1)*A_R*(h_Delta^gamma)*g/(2*pi))^2;
    d_u_max = (radious);
    d_r_max = radious;
    func_d = @(d) (1./(d_u_max - d_u_min)) .* double(d >= d_u_min) .* double(d <= d_u_max);
    func_h_d2 = @(y) (y.^(-(gamma+3)./(gamma+2))  .* (y.^(-1/(gamma+2)) - h_Delta^2).^(-0.5) .* func_d((y.^(-1./(gamma+2)) - h_Delta^2).^0.5 )).*double(y >= 0).*double(y <= h_Delta^(-2*(gamma+2)));
    c_d = 1./integral(func_h_d2, 0, h_Delta^(-2*(gamma+2)));
    Func_theta = @(z) ((z - theta_1)./(theta_2-theta_1)).*double(z >= theta_1).*double(z <= theta_2) + double(z > theta_2);
    func_theta = @(z) (1./(theta_2-theta_1)).*double(z>=theta_1).*double(z<=theta_2);
    func_delta = @(a, b) (Func_theta(b) - Func_theta(a)).*double(a <= b);
    Func_h_theta2 = @(t) (func_delta(0.5.*acos(2.*t-1), Theta) + 1 - Func_theta(Theta)).*double(t>=0).*double(t<=1) + double(t >1);
    func_h_theta2 = @(t) (1./sqrt( 4.*t.*(1-t) )  ).*func_theta(0.5.*acos(2.*t-1)).*double((cos(Theta))^2 <= t).*double(t <1);

    for cnt_i_b = 1:length(I_b_total)

        I_b = I_b_total(cnt_i_b);
        A_VLC = I_max - I_b;
        x = (2^(R_th/B_VLC) - 1)*(2*pi*sigma_2)/(exp(1)*(eta*P_LED*A_VLC)^2);
        xh = (2^(R_th/B_VLC) - 1)*(2*pi*sigma_2)/(exp(1)*(eta*P_LED*A_VLC)^2)/h_c_prime2;
        Func_vlc_h2 = @(yy) (1) .* ...
            func_h_d2(yy) .*     (Func_h_theta2(xh./(yy)));

        pr_out_anl(cnt_th, cnt_i_b) =  1 - ( ((Theta - theta_1)/(theta_2 - theta_1)) .* double(Theta < theta_2) + double(Theta > theta_2))  + ...
            (c_d) * integral(Func_vlc_h2,  0, h_Delta^(-2*(gamma+2)));

        func_vlc_h2 = @(yyy, xxx)   (1/h_c_prime2).*(1./yyy).*func_h_d2(yyy).*func_h_theta2(  xxx./(h_c_prime2.*yyy)   );
        x_min = @(yyy)  (cos(Theta))^2.*h_c_prime2.*yyy;
        x_max = @(yyy) h_c_prime2.*yyy;
        c_h = 1/integral2(func_vlc_h2, 0, h_Delta^(-2*(gamma+2)), x_min, x_max);

        func_vlc_h2_prime = @(yyyy, xxxx) ((xxxx./(yyyy)).*func_h_d2(yyyy).*func_h_theta2(xxxx./(yyyy.*h_c_prime2)));
        x_min = @(yyyy)  (cos(Theta))^2.*h_c_prime2.*yyyy;
        x_max = @(yyyy) h_c_prime2.*yyyy;
        f =  0.75.*eta.*P_LED.*I_b.*V_t.*eta.*P_LED.*I_b./(I_0);
        E_h_anl(cnt_th, cnt_i_b) = f.*c_h.*integral2(func_vlc_h2_prime ,  0, h_Delta^(-2*(gamma+2)), x_min, x_max)/h_c_prime2;
        f =  eta.*P_LED.*I_b;
        h_VLC_min = (gamma + 1).*A_R.*(h_Delta^gamma).*g.*((h_Delta^2 + d_r_max.^2).^(-(gamma+2)/2))*cos(theta_2)./(2*pi);
        h_VLC_max = (gamma + 1).*A_R.*(h_Delta^gamma).*g.*((h_Delta^2).^(-(gamma+2)/2))*cos(theta_1)./(2*pi);
        E_h_min = 0.75*V_t*h_VLC_min*f*log(1+f*h_VLC_min/I_0);
        E_h_max = 0.75*V_t*h_VLC_max*f*log(1+f*h_VLC_max/I_0);
        alpha=0.95;
        E_h_alpha = alpha*E_h_min + (1-alpha)*E_h_max;
        radius = radious;
        f_1 = @(t) double(t.^(-1/beta) > 0) .* double(t.^(-1/beta) < 2*radius ) .* 4./(pi*radius^2) .* (pi.*(t.^(-1/beta)).^2 ./4 - (2.*(t.^(-1/beta)).^2 - 4*radius.^2) ./4 .* asin((t.^(-1/beta))./(2.*radius)) - (t.^(-1/beta))./4 .* sqrt(4*radius^2 - (t.^(-1/beta)).^2) );
        f_2 = @(t) double(t.^(-1/beta) > 0) .* double(t.^(-1/beta) < 2*radius ) .*2./(pi*radius) .* (-(t.^(-1/beta)) .* sqrt((1-(t.^(-1/beta)).^2./(4*radius^2)).^3  )  + (t.^(-1/beta))./2 .*sqrt(1- (t.^(-1/beta)).^2./(4.*radius^2)) - radius .* asin(-(t.^(-1/beta))./(2*radius))  );

        xf = (2^(R_th/B_RF) - 1)*(N_0*(4*pi*F_c)^2)/(E_h_sim(cnt_th, cnt_i_b)*c_light^2);
        sigma_h2 = 0.5;
        al_func = @(y) (1 - (  f_1(xf./y) - f_2(xf./y) )    ) .* exp(-y./(2*sigma_h2))  ;
        int_val(cnt_th, cnt_i_b) =  1/(2*sigma_h2) .* integral(al_func, 0, xf*(2*radius)^beta);
        xf_max = (2^(R_th/B_RF) - 1)*(N_0*(4*pi*F_c)^2)/(E_h_min*c_light^2);
        sigma_h2 = 0.5;
        al_func = @(y) (1 - (  f_1(xf_max./y) - f_2(xf_max./y) )    ) .* exp(-y./(2*sigma_h2))  ;
        int_val_max(cnt_th, cnt_i_b) =  1/(2*sigma_h2) .* integral(al_func, 0, xf_max*(2*radius)^beta);

        xf_min = (2^(R_th/B_RF) - 1)*(N_0*(4*pi*F_c)^2)/(E_h_max*c_light^2);
        al_func = @(y) (1 - (  f_1(xf_min./y) - f_2(xf_min./y) )    ) .* exp(-y./(2*sigma_h2))  ;
        int_val_min(cnt_th, cnt_i_b) =  1/(2*sigma_h2) .* integral(al_func, 0, xf_min*(2*radius)^beta);

        xf_alpha = (2^(R_th/B_RF) - 1)*(N_0*(4*pi*F_c)^2)/(E_h_alpha*c_light^2);
        al_func = @(y) (1 - (  f_1(xf_min./y) - f_2(xf_min./y) )    ) .* exp(-y./(2*sigma_h2))  ;
        int_val_alpha(cnt_th, cnt_i_b) =  1/(2*sigma_h2) .* integral(al_func, 0, xf_alpha*(2*radius)^beta);

        Pr_out_last_theo_min(cnt_th, cnt_i_b) = pr_out_anl(cnt_th, cnt_i_b) * (pr_out_anl(cnt_th, cnt_i_b) + (1-pr_out_anl(cnt_th, cnt_i_b))*int_val_min(cnt_th, cnt_i_b));

        Pr_out_last_theo_max(cnt_th, cnt_i_b) = pr_out_anl(cnt_th, cnt_i_b) * (pr_out_anl(cnt_th, cnt_i_b) + (1-pr_out_anl(cnt_th, cnt_i_b))*int_val_max(cnt_th, cnt_i_b));

    end
end

% Plot the results
figure
plot(threshold_total, pr_out_vlc_rf(:,1), 'o', 'Color', colors(2,:), 'MarkerIndices',1:20:length(threshold_total))

hold on
grid on
plot(threshold_total, Pr_out_last_theo_max(:,1), '+--', 'Color', colors(1,:), 'MarkerIndices',1:20:length(threshold_total))
plot(threshold_total, Pr_out_last_theo_min(:,1), '-.', 'Color', colors(7,:))
plot(threshold_total, pr_out_vlc(:,1), '-', 'Color', colors(4,:))

plot(threshold_total, Pr_out_last_theo_max(:,2), '+--', 'Color', colors(1,:), 'MarkerIndices',1:20:length(threshold_total))

plot(threshold_total, pr_out_vlc_rf(:,2), 'o', 'Color', colors(2,:), 'MarkerIndices',1:20:length(threshold_total))
plot(threshold_total, Pr_out_last_theo_min(:,2), '-.', 'Color', colors(7,:))
plot(threshold_total, pr_out_vlc(:,2), '-', 'Color', colors(4,:))
legend({'$P^{\mathcal{O}}_{\rm H-VR}$ (simulation)', '$P^{\mathcal{O}}_{\rm H-VR, UB}$ (analytic)', '$P^{\mathcal{O}}_{\rm H-VR, LB}$ (analytic)', '$P^{\mathcal{O}}_{\rm{direct}}$ (analytic)'},'Interpreter','latex')
xlabel({'$R_{\textrm{th}}$ (bps)'},'Interpreter','latex')
ylabel('Outage probability')
xlim([threshold_total(1) threshold_total(end)])
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