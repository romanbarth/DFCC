%
% Dense Flow reConstruction and Correlation (DFCC)
% ----------------------------------------------------------------------- %
%
% Reference to the publication:
%   Haitham A Shaban, Roman Barth, Kerstin Bystricky; Formation of correlated 
%   chromatin domains at nanoscale dynamic resolution during transcription, 
%   Nucleic Acids Research, gky269, https://doi.org/10.1093/nar/gky269
%
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
