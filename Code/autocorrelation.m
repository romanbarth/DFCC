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

function [R, lags] = autocorrelation(dir_mag, pixelsize, mask, xp, yp)
% AUTOCORRELATION: computes the autocorrelation for direction or magnitude
% of the flow fields at all available time lags. Takes trajectories and 
% calculates the angle at different time lags as complex number. Afterwards
% the spatial correlation is calculated and radially averaged
%
%   INPUT
%   dir_mag:   string to indicate if autocorrelation of direction or
%              magnitude shall be computed
%                    'dir': direction
%                    'mag': magnitude
%   pixelsize: pixelsize in micrometer
%   mask:      logical mask containing zeros outside nucleus and inside 
%              nucleoli and ones otherwise
%   xp:        3D array with size of u and length = length(u)+1. The array 
%              consists of the x-position of every pixel for time t
%   yp:        3D array with size of u and length = length(u)+1. The array 
%              consists of the y-position of every pixel for time t.
%
%   OUTPUT
%   R:         cell-vector containing spatial correlation over space lag
%              for every time lag over time. Each entry may contain
%              several curves, depending on the space lag. The function
%              uses the average of the provided curves for each space
%              lag
%   lags:      time lags at which the correlation was evaluated
%
% ----------------------------------------------------------------------- %

% initialize
[X, Y] = meshgrid(-(size(xp,2)*2-1)/2:(size(xp,2)*2-1)/2-1, ...
    -(size(xp,1)*2-1)/2:(size(xp,1)*2-1)/2-1); % Make Cartesian grid
[~, rho] = cart2pol(X, Y);   % Convert to polar coordinate axes
lags = 0:max(rho(:));

R = cell(1, size(xp,3)-1);

% find inner circle of rescaled mask
mask = double(mask);
maskc = innerCircle(round(imresize(mask, size(mask)*2-1)));
% set outer pixels to NaN. These values are not taken into account in
% averaging, because we use nanmean and cut NaN values at the end.
maskc(maskc==0) = NaN;

% loop through lags
for lag = 1:size(xp,3)-1       
    
    if strcmp(dir_mag, 'mag')
        disp(['Autocorrelation in magnitude: ', ...
            num2str(lag/(size(xp,3)-1)*100), '%'])
        % calculate magnitude
        arg = ( (yp(:,:,1+lag:end)-yp(:,:,1:end-lag)).^2 + ...
            (xp(:,:,1+lag:end)-xp(:,:,1:end-lag)).^2 ).^(1/2);
        
    elseif strcmp(dir_mag, 'dir')
        disp(['Autocorrelation in direction: ', ...
            num2str(lag/(size(xp,3)-1)*100), '%'])
        % calcualte angle in radians
        arg = atan2(yp(:,:,1+lag:end)-yp(:,:,1:end-lag), ...
            xp(:,:,1+lag:end)-xp(:,:,1:end-lag));
    end
    
    C = zeros(size(arg,1)*2-1, size(arg,2)*2-1, size(arg,3));
    for k = 1:size(arg,3)
        
        if strcmp(dir_mag, 'mag')
            % apply mask
            z = arg(:,:,k) .* mask;
            
            % substract the mean
            z(z==0) = NaN;
            zm = nanmean(z(:));            
            z = z-zm;            
            z(isnan(z)) = 0;
            
        elseif strcmp(dir_mag, 'dir')
            % translate to complex plane
            z = exp(1i*arg(:,:,k)) .* mask;
        end
        
        % calculate cross-correlation and normalize by the energy of the 
        % argument
        C(:,:,k) = xcorr2(z, z)/sum(abs(z(:)).^2); % not normalized
        
        % crop correlation function by rescaled version of circle shaped
        % mask
        C(:,:,k) = C(:,:,k).*maskc;
        
        % radial average
        [out_radavg, ~] = radavg(real(C(:,:,k)), pixelsize);
        if k > 1
            if length(out_radavg) >= length(R{lag}(k-1,:))
                R{lag}(k,:) = out_radavg(1:length(R{lag}(k-1,:)));
            elseif length(out_radavg) < length(R{lag}(k-1,:))
                R{lag} = R{lag}(:,1:length(out_radavg));
                R{lag}(k,:) = out_radavg;
            end
        else
            R{lag}(k,:) = out_radavg;
        end              
        
    end
end

% crop to standard length: find smallest length
lengthR = cellfun(@(a) size(a,2), R, 'UniformOutput', true);
lengthMin = min([lengthR length(lags)]);
% crop
R = cellfun(@(a) a(:,1:lengthMin), R, 'UniformOutput', false);
lags = lags(1:lengthMin);
