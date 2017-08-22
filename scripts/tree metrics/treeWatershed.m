function label = treeWatershed(chm, varargin)
% TREEWATERSHED - extract Individual Tree Crowns (ITC) label from a raster Canopy Height Model (CHM) 
% using (optional marker-controlled) watershed segmentation.
%
% LABEL = TREEWATERSHED produces the label matrix LABEL containing individual tree crown labels 
% for each pixel using the approach described in Kwak et al. (2007) [1]. The function takes as input the raster Canopy
% Height Model CHM and the image coordinates MARKERS of individual tree tops.
%
% [1] D.-A. Kwak, W.-K. Lee, J.-H. Lee, G. S. Biging, and P. Gong, 
% "Detection of individual trees and estimation of tree height using LiDAR data,"
% J For Res, vol. 12, no. 6, pp. 425–434, Oct. 2007.
%
% Syntax:
%    label = treeWatershed(chm, ...)
%
% Inputs:
%    chm - RxC numeric matrix, raster Canopy Height Model (CHM) with R rows and C columns
%    
%    markers (optional, default: []) - Mx3 numeric matrix, images coordinates (col, row) and height 
%    values of markers used in the watershed algorithm
%
%    minHeight (optional, default: 1) - numeric value, minimum canopy height, all values 
%    below this threshold are set to zero.
%
%    verbose (optional, default: true) - boolean value, verbosiy switch
%
%    fig (optional, default: true) - boolean value, switch to plot figures
%
% Outputs:
%    label - RxC numeric matrix, individual tree crown labels
%
% Example:
%    info = geotiffinfo('..\data\measurements\raster\chm\so_2014_woodland_pasture.tif');
%    [chm, ~] = geotiffread('..\data\measurements\raster\chm\so_2014_woodland_pasture.tif');
%    
%    [crh, xyh] = canopyPeaks(chm, ...
%        info.RefMatrix, ...
%        'adjacencyRange', 1.5, ...
%        'method', 'allometricRadius', ...
%        'allometry', @(h) 1 + 0.5*log(max(h,1)), ...
%        'fig', true);
%
%    label = treeWatershed(chm, ...
%        'markers', crh, ...
%        'minHeight', 1, ...
%        'fig', true, ...
%        'verbose', true);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2016b
%
% See also: canopyPeaks.m, canopyCover.m
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory (LASIG)
% Website: http://mparkan.github.io/Digital-Forestry-Toolbox/
% Last revision: August 22, 2017
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund (WHFF, OFEV), project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'chm', @isnumeric);
addParameter(arg, 'markers', [], @isnumeric);
addParameter(arg, 'minHeight', 1, @(x) isnumeric(x) && (numel(x) == 1));
addParameter(arg, 'fig', true, @(x) islogical(x) && (numel(x) == 1));
addParameter(arg, 'verbose', true, @(x) islogical(x) && (numel(x) == 1));

parse(arg, chm, varargin{:});


%% compute gradient magnitude

if arg.Results.verbose
    
    fprintf('computing gradient magnitude...');
    
end

chm(chm <= arg.Results.minHeight) = 0;

hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(chm), hy, 'replicate');
Ix = imfilter(double(chm), hx, 'replicate');
gradient_magnitude = sqrt(Ix.^2 + Iy.^2);

% gradient mask
gradient_mask = false(size(chm));
gradient_mask(gradient_magnitude > 1.2) = 1;

% fill holes in gradient mask
connected_components = bwconncomp(~gradient_mask);
numPixels = cellfun(@numel, connected_components.PixelIdxList);
gradient_mask(cell2mat(connected_components.PixelIdxList(numPixels < 10)')) = 1; 

if arg.Results.verbose
    
    fprintf('done!\n');
    
end

   
%% compute watershed transform

if arg.Results.verbose
    
    fprintf('computing watershed transform...');
    
end

if isempty(arg.Results.markers) 
    
    % basic watershed transform
    label = watershed(-chm, 8);
    
else 
    
    % marker controled watershed transform
    ind_peaks = sub2ind(size(chm), arg.Results.markers(:,2), arg.Results.markers(:,1));
    marker = false(size(chm));
    marker(ind_peaks) = 1;
    I = imimposemin(-chm, marker, 8);
    label = watershed(I, 8);
    
end

if arg.Results.verbose
    
    fprintf('done!\n');
    
end

%% filter non vegetation pixels

label(~gradient_mask) = 0;


%% plot results

if arg.Results.fig
    
    % assign distinct colors to adjacent clusters 
    [nrows, ncols] = size(label);
    [xgrid, ygrid] = meshgrid(1:ncols, 1:nrows);
    buffer = min(10, 3*sqrt(median(accumarray(label(:)+1, label(:), [], @numel))));

    [color, ~, cmap] = clusterColor([xgrid(:), ygrid(:)], ...
        label(:), ...
        'adjacency', '2d', ...
        'buffer', buffer, ...
        'colormap', 'cmap12', ...
        'unlabelledColor', [0.1, 0.1, 0.1], ...
        'fig', false, ...
        'verbose', false);
    
    figure
    imagesc(reshape(color, nrows, ncols))
    axis equal tight
    xlabel('col');
    ylabel('row');
    title('Individual tree crown labels')
    colormap(cmap);
    
end
