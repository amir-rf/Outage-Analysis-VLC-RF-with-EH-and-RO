function [d_u_min, d_u_max, d_r_min, d_r_max, N_F, B, P_0, F_c, ...
    I_min, I_max, eta, V_t, I_0, Phi, Theta, q_e, P_LED, I_i, ...
    A_R, h_Delta, B_RF, B_VLC, R_th] = initialization()


%% Initialization
%User distance
d_u_min = 0;
d_u_max = 6;
%Relay distance
d_r_min = 0;
d_r_max = 6;
%Noise figure
N_F = 9; %dB
% Signal bandwidth
B = 10e6; % 10 MHz
% Thermal noise
P_0 = -174; %dBm/Hz
P_0 = P_0 - 30; %dBm to dB
% Frequency
F_c = 2.4e9; %GHz
% F_c = 5e9; %GHz
% Minimum DC bias
I_min = 100e-3;
% Maximum DC bias
I_max = 1;
% Photo-detector responsivity A/W
eta = 0.4;
% Thermal voltage
V_t = 25e-3;
% Dark saturation current
I_0 = 1e-7;
% Half FOV
Phi = deg2rad(60);
% Half-power beamwidth
Theta = deg2rad(60);
% Electron charge
q_e = 1.6e-19;
% Induced current
I_i = 5840*1e-6;
% PD detection area
A_R =1e-4;
% AP relative height
h_Delta = 2;


B_RF = B;
B_VLC = B;
P_LED = 1.5; %W/A
% Rate of threshold
R_th = 1e6;


end