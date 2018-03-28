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

function initialize
% INITIALIZE adds required subfolders to matlab path

disp('Initialization...')

here = mfilename('fullpath');
[path, ~, ~] = fileparts(here);
addpath(genpath(path));


end
