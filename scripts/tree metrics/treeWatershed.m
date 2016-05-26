function labels = treeWatershed(chm, markers, varargin)
% TREEWATERSHED - extract Individual Tree Crowns (ITC) labels from a raster Canopy Height Model (CHM) 
% using marker-controlled watershed segmentation.
%
% LABELS = TREEWATERSHED produces the label matrix LABELS containing individual tree crown labels 
% for each pixel using the approach described in Kwak et al. (2007) [1]. The function takes as input the raster Canopy
% Height Model CHM and the image coordinates MARKERS of individual tree tops.
%
% [1] D.-A. Kwak, W.-K. Lee, J.-H. Lee, G. S. Biging, and P. Gong, 
% "Detection of individual trees and estimation of tree height using LiDAR data,"
% J For Res, vol. 12, no. 6, pp. 425–434, Oct. 2007.
%
% Syntax:  labels = treeWatershed(chm, markers, ...)
%
% Inputs:
%    chm - RxC numeric matrix, raster Canopy Height Model (CHM) with R rows and C columns
%
%    markers - Mx3 numeric matrix, images coordinates (col, row) and height values of tree tops
%
%    minHeight (optional, default: 1) - numeric value, minimum canopy height, all values 
%    below this threshold are set to zero.
%
%    verbose (optional, default: true) - boolean value, verbosiy switch
%
%    fig (optional, default: true) - boolean value, switch to plot figures
%
% Outputs:
%    labels - RxC numeric matrix, individual tree crown labels
%
% Example:
%    info = geotiffinfo('..\data\measurements\raster\chm\so_2014_woodland_pasture.tif');
%    [chm, ~] = geotiffread('..\data\measurements\raster\chm\so_2014_woodland_pasture.tif');
%    
%    [crh, xyh] = treePeaks(chm, info.RefMatrix, ...
%        'method', 'fixedRadius', ...
%        'windowRadius', 3, ...
%        'fig', true);
%
%    labels = treeWatershed(chm, crh, ...
%        'minHeight', 1, ...
%        'fig', true);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2016a
%
% See also: treePeaks.m, canopyCover.m
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory
% Website: http://lasig.epfl.ch/
% Last revision: May 26, 2016
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund (WHFF, OFEV), project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'chm', @isnumeric);
addRequired(arg, 'markers', @isnumeric);
addParamValue(arg, 'minHeight', 1, @(x) isnumeric(x) && (numel(x) == 1));
addParamValue(arg, 'fig', true, @(x) islogical(x) && (numel(x) == 1));
addParamValue(arg, 'verbose', true, @(x) islogical(x) && (numel(x) == 1));

parse(arg, chm, markers, varargin{:});


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

ind_peaks = sub2ind(size(chm), markers(:,2), markers(:,1));
marker = false(size(chm));
marker(ind_peaks) = 1;
labels = watershed(bwdist(marker), 8);

if arg.Results.verbose
    
    fprintf('done!\n');
    
end


%% match watershed labels with tree top indices

if arg.Results.verbose
    
    fprintf('matching watershed labels with marker indices...');
    
end

peak_labels = [0; labels(ind_peaks)];

[~,idx] = ismember(labels, 0:size(markers,1));
labels = peak_labels(idx);

if arg.Results.verbose
    
    fprintf('done!\n');
    
end


%% filter non forest pixels

labels(~gradient_mask) = 0;

if arg.Results.fig
    
    figure
    imagesc(labels)
    axis equal tight
    xlabel('col');
    ylabel('row');
    title('Individual tree crown labels')
    
end

