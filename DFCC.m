%
% Interface to apply Dense Flow reConstruction and Correlation (DFCC)
% ----------------------------------------------------------------------- %
%
% Reference to the publication:
%   Shaban, H.A.; Barth, R.; Bystricky, K. title. Nucleic Acids Research, volume(issue):pages1-pages2, 2018
%
% Demo for example data (U2OS cells expressing H2B-GFP of size 348 pixels
% x 338 pixels x 50 frames)
% The data is a part of the entire dataset of size 1024 pixels x 1024
% pixels x 150 frames (88 nm pixel size and 200 ms per frame), the same as
% shown in Figure 2, Figure 3 and Figure 5 of the manuscript.
%
% The example data takes ~1h 40 min to process on following system:
%       Intel(R) Core(TM) i7-6700IQ CPU @ 2.60 GHz 2.60 Ghz
%       16.0 Gb RAM
%
% developed at:
%       Laboratoire de Biologie Moléculaire Eucaryote (LBME),
%       Centre de Biologie Intégrative (CBI), CNRS;
%       University of Toulouse, UPS; 31062
%       Toulouse; France
%
% ----------------------------------------------------------------------- %

close all
clearvars

%% Set parameters
% give the path to the file to analyse
filepath = pwd;
if ispc
    directory = [filepath '\U2OS_H2BGFP_example_data.tif'];
else
    directory = [filepath '/U2OS_H2BGFP_example_data.tif'];
end

% pixel size in micrometer
pixelsize= 0.088;

% acquisition time per frame in seconds
dT = 0.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             USER INPUT                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Search path
initialize % add all required paths to the Matlab path

%% save date and time for later
cstart = datestr(now); % get start time
datee = datestr(now);
datee(datee==':') = '-';

%% load images, set mask and filter noise (blurring)
[im, numIm, mask, mask_nuc] = Preprocessing(directory);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            START SCRIPT                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Starting calculations...')

%% calculate velocity field
[x, y, u, v] = OpticalFlow(im, mask);
% apply nucleolus mask and convert displacement to SI units
u = cellfun(@(x) x .* mask_nuc * pixelsize, u, 'uniformoutput', false);
v = cellfun(@(x) x .* mask_nuc * pixelsize, v, 'uniformoutput', false);

%% integrate velocity field
disp('Integration...')
[xp, yp] = IntegrateFlowField(x, y, u, v, mask_nuc);

%% compute spatial autocorrelation for direction and magnitude and fit
[R,     lags] = autocorrelation('dir', pixelsize, mask.msk, xp, yp);
[R_mag, ~   ] = autocorrelation('mag', pixelsize, mask.msk, xp, yp);

% fit spatial autocorrelation
[xi, nu]         = AutocorrelationFit(lags, R);
[xi_mag, nu_mag] = AutocorrelationFit(lags, R_mag);


%% save data
if ispc
    save([filepath '\DFCC_analysis_' convertDate(datee) '.mat'])
else
    save([filepath '/DFCC_analysis_' convertDate(datee) '.mat'])
end

%% plot resulting parameters with standard deviation
f = PlotParameters(xi, nu, xi_mag, nu_mag, dT);

cend = datestr(now);
disp('Analysis finished!')
disp(['Started:  ' cstart])
disp(['Finished: ' cend])
% ----------------------------------------------------------------------- %