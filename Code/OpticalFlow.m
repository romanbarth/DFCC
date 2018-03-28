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

function [x,y,u,v] = OpticalFlow(im, mask)
% OPTICALFLOW: computes Optical Flow for given image series using a Optical
% Flow algorithm by Sun et al:
%
%   Sun,D., Roth,S. and Black,M.J. (2014) A quantitative analysis of 
%   current practices in optical flow estimation and the principles behind 
%   them. Int. J. Comput. Vis., 106, 115–137
%
%   INPUT
%   im:   cell-vector containing the time series frames. Each cell entry is
%         a 2D array
%   mask: structure with fields .msk, .idx and .idy as given by maskROI.m
%
%   OUTPUT
%   x:    array of size of the frames in im with x-indices of pixels
%   y:    array of size of the frames in im with y-indices of pixels
%   u:    cell vector with length length(im)-1 with the x-component of the
%         computed flow field 
%   v:    cell vector with length length(im)-1 with the y-component of the
%         computed flow field 
%
% ----------------------------------------------------------------------- %

% initialize
numIm = length(im);
u = cell(1, numIm-1);
v = u;

% create mesh on which flow field is evaluated, i.e. the centers of each image pixel
[x,y] = meshgrid(1:size(im{1},2), 1:size(im{1},1));

% Optical Flow
disp('Optical Flow..')
for j = 1:numIm-1
    
    % apply Optical Flow to every pair of subsequent frames
    uv = estimate_flow_interface(im{j}, im{j+1}, mask.msk, 'classic+nl-fast');
    u{j} = double(uv(:,:,1));
    v{j} = double(uv(:,:,2));
    
    % apply mask (including nucleoli to flow field
    u{j} = u{j}.*mask.msk;
    v{j} = v{j}.*mask.msk;
    
    % filter to suppress artifacts at borders
    [u{j}, v{j}] = FilterFlowField(u{j}, v{j}, 6);
    
    % interpolate resulting NaNs
    if sum(isnan(u{j}(:))) > 0
        u{j} = inpaint_nans(u{j}, 4);
    end
    if sum(isnan(v{j}(:))) > 0
        v{j} = inpaint_nans(v{j}, 4);
    end
    
    disp(['Optical Flow: ', num2str(round(j/(numIm-1)*100)) '%'])
end
