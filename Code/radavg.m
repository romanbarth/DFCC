%
% Dense Flow reConstruction and Correlation (DFCC)
% ----------------------------------------------------------------------- %
%
% 
% Reference to the publication:
%   Haitham A Shaban, Roman Barth, Kerstin Bystricky; Formation of correlated 
%   chromatin domains at nanoscale dynamic resolution during transcription, 
%   Nucleic Acids Research, gky269, https://doi.org/10.1093/nar/gky269
%
% developed at:  
%       Laboratoire de Biologie Moléculaire Eucaryote (LBME), 
%       Centre de Biologie Intégrative (CBI), CNRS; 
%       University of Toulouse, UPS; 31062 
%       Toulouse; France
%
% ----------------------------------------------------------------------- %

function [Radavg, lags] = radavg(C, pixelsize)
% RADAVG calculates the radial average of array C
%
%   INPUT
%   C:  input array
%
%   OUTPUT
%   Radavg:    Radial average vector
%   lags:      vector containing space lags
%
% ----------------------------------------------------------------------- %

% Compute radially average of power spectrum
% Make Cartesian grid
[X, Y] = meshgrid(-size(C,2)/2:size(C,2)/2-1, -size(C,1)/2:size(C,1)/2-1);
% Convert to polar coordinate axes
[~, rho] = cart2pol(X, Y);
rho = round(rho);
lags = 0:max(rho);
Radavg = nan(1, length(lags)-1);

xc = (size(C,1)+1)/2;
yc = (size(C,2)+1)/2;

for r = 2:length(lags)

    if min(xc-lags(r))<1 || min(yc-lags(r))<1
        tempC = C;
        tempRho = rho;
    else
        tempC = C(xc-lags(r):xc+lags(r), yc-lags(r):yc+lags(r));
        tempRho = rho(xc-lags(r):xc+lags(r), yc-lags(r):yc+lags(r));
    end
    Radavg(r) = nanmean(tempC(( (tempRho <= lags(r)) & (tempRho > lags(r-1)) )));
        
    % if the number is NaN, break. We reached the border of the circle
    if isnan(Radavg(r))
        break;
    end
end
Radavg(isnan(Radavg)) = [];

Radavg(1) = C((size(C,1)+1)/2, (size(C,2)+1)/2);
lags = lags(1:length(Radavg))*pixelsize;
