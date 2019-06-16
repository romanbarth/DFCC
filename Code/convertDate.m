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

function out = convertDate(in)
% CONVERTDATE converts date into a form without dashes and spaces such that
% names including the date are readable platform- and program-independent.
%
%   INPUT
%   in:      date to be processed as string
%
%   OUTPUT
%   out: modified date given as string
%
% ----------------------------------------------------------------------- %

out = strrep(in, ' ', '_');
out = strrep(out, '-', '');
out = [out(1:12) 'h' out(13:14) 'm'];