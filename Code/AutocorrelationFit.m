%
% Dense Flow reConstruction and Correlation (DFCC)
% ----------------------------------------------------------------------- %
%
% Reference to the publication:
%   Shaban, H.A.; Barth, R.; Bystricky, K. title. Nucleic Acids Research, volume(issue):pages1-pages2, 2018
%
% developed at:  
%       Laboratoire de Biologie Moléculaire Eucaryote (LBME), 
%       Centre de Biologie Intégrative (CBI), CNRS; 
%       University of Toulouse, UPS; 31062 
%       Toulouse; France
%
% ----------------------------------------------------------------------- %

function [xi, nu] = AutocorrelationFit(lags, Correlation)
% AUTOCORRELATIONFIT: fits the empirical autocorrelaiton to the
% Whittle-Matern model using MATLAB's lsqcurvefit and returns fit parameter
% which standard deviation
%
%   INPUT
%   lags:         vector containing spatial lags
%   Correlation:  cell-vector containing spatial correlation over space lag
%                 for every time lag over time. Each entry may contain
%                 several curves, depending on the space lag. The function
%                 uses the average of the provided curves for each space
%                 lag
%
%   OUTPUT
%   xi:           array of size [lags, 2], where xi(lag,1) contains the
%                 correlation length for given lag and xi(lag,2) contains 
%                 the corresponding standard deviation
%   nu:           array of size [lags, 2], where nu(lag,1) contains the
%                 smoothness parameter for given lag and nu(lag,2) contains 
%                 the corresponding standard deviation
%
% ----------------------------------------------------------------------- %

% lags to consider
indices = 1:length(Correlation);

% initialize
xi    = zeros(length(indices), 2);
nu    = zeros(length(indices), 2);

% Whittle-Matern model
params.fun = @(b,x) b(1) .* 2^(1-b(3))./gamma(b(3)) .* (x./b(2)).^b(3) .* besselk(b(3), x./b(2));
params.start = [1    3     3];
params.ub    = [2000 10000 20];
params.lb    = [eps  0     0];

% shift 0-lag as function is not defined at r=0
if lags(1) == 0
    lags(1) = eps;
end

% fit for every time lag
for lag = 1:length(indices)
    
    y = mean(Correlation{lag}, 1); % columns are different observations
    
    % fit
    [coeffs, resnorm, ~, ~, ~, ~, J] = lsqcurvefit(params.fun, params.start, lags, y, params.lb, params.ub, optimset('Display','off'));
    % Jacobian
    J = full(J);
    % covariance
    COVB = resnorm/(length(y)-1)*(J'*J)^(-1);
    % variance are diagonal elements
    coeffs_std = sqrt(diag(COVB));
        
    % extract fit parameter and standard deviation
    xi(lag,1) = coeffs(2);
    xi(lag,2) = coeffs_std(2);
    nu(lag,1) = coeffs(3);
    nu(lag,2) = coeffs_std(3);
end