%
% Dense Flow reConstruction and Correlation (DFCC)
% ----------------------------------------------------------------------- %
%
% Reference to the publication:
%   Shaban, H.A.; Barth, R.; Bystricky, K. title. Nucleic Acids Research, volume(issue):pages1-pages2, 2018
%
% developed at:  
%       Laboratoire de Biologie Mol�culaire Eucaryote (LBME), 
%       Centre de Biologie Int�grative (CBI), CNRS; 
%       University of Toulouse, UPS; 31062 
%       Toulouse; France
%
% ----------------------------------------------------------------------- %

function im = readTiffStack(directory)
% READTIFFSTACK reads a stack of tiff images
%
%   INPUT
%   directory:  string providing the full path, filename and extension of
%               the file to read in
%
%   OUTPUT
%   im:         cell-vector containing the time series frames. Each cell 
%               entry is a 2D array
%
% ----------------------------------------------------------------------- %

warning('off')
InfoImage = imfinfo(directory);
numIm = length(InfoImage);
im = cell(1, numIm);

obj = Tiff(directory, 'r');
for i = 1:numIm
    
    obj.setDirectory(i);
    im{i} = obj.read();
    
end
obj.close();