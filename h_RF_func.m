function h = h_RF_func(type, K, num)
% h_RF_func generates Rayleigh or Rician fading samples
% type = 1: Rayleigh fading
% type = 2: Rician fading with factor K
% num: number of complex fading samples

if type == 1
    % Rayleigh fading: zero-mean complex Gaussian (variance 1)
    h = (randn(1, num) + 1i * randn(1, num)) / sqrt(2);
    
elseif type == 2
    % Rician fading with K-factor
    % g1 = mean (LOS component), g2 = std. dev. of NLOS component
    g1 = sqrt(K / (K + 1));
    g2 = sqrt(1 / (2 * (K + 1)));
    
    h = (g2 * randn(1, num) + g1) + 1i * (g2 * randn(1, num));
    
else
    error('Invalid type. Use 1 for Rayleigh or 2 for Rician.');
end
end