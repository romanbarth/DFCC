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

function f = PlotParameters(xi_dir, nu_dir, xi_mag, nu_mag, dT)
% PLOTPARAMETERS: plots the correlation length xi and smoothness parameter
% nu as computed by AutocorrelationFit.m over time lag
%
%   INPUT
%   xi_dir: correlation length for the direction in the format as given by 
%           AutocorrelationFit.m
%   nu_dir: smoothness parameter for the direction in the format as given by 
%           AutocorrelationFit.m
%   xi_mag: correlation length for the magnitude in the format as given by 
%           AutocorrelationFit.m
%   nu_mag: smoothness parameter for the magnitude in the format as given by 
%           AutocorrelationFit.m
%
%   OUTPUT
%   f:      handle to the figure
%
% ----------------------------------------------------------------------- %

% settings
linewidth = 2;
linecolor = ['r', 'g'];

% create time lags
timelag = (1:size(xi_dir,1))*dT;

% create figures with 2 subplots
f = figure('unit', 'normalized', 'position', [0.01 0.4 0.7 0.4]); hold on
set(f, 'Color', 'w')

% correlation length
subplot(1,2,1), hold on
h(1) = errorbar(timelag, xi_dir(:,1), xi_dir(:,2));
h(2) = errorbar(timelag, xi_mag(:,1), xi_mag(:,2));
set(gca, 'FontName', 'Arial', 'FontSize', 18)

for ih = 1:length(h)
    set(h(ih), 'color', linecolor(ih), 'linewidth', linewidth)
end
legend(h, {'Direction', 'Magnitude'})
xlabel('Time lag \Deltat [s]')
ylabel('Correlation length \xi [\mum]')

% smoothness parameter
subplot(1,2,2), hold on
h(1) = errorbar(timelag, nu_dir(:,1), nu_dir(:,2));
h(2) = errorbar(timelag, nu_mag(:,1), nu_mag(:,2));
set(gca, 'FontName', 'Arial', 'FontSize', 18)

for ih = 1:length(h)
    set(h(ih), 'color', linecolor(ih), 'linewidth', linewidth)
end
xlabel('Time lag \Deltat [s]')
ylabel('Smoothness parameter \nu')
