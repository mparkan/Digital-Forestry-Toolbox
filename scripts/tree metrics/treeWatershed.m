function varargout = treeWatershed(chm, varargin)
% TREEWATERSHED - extract Individual Tree Crowns (ITC) label from a raster Canopy Height Model (CHM)
% using (optional marker-controlled) watershed segmentation.
%
% [LABEL, ...] = TREEWATERSHED produces the label matrix LABEL containing individual tree crown labels
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
%    minHeight (optional, default: 1) - numeric value, minimum canopy height, all CHM values
%    below this threshold are set to zero
%
%    mask (optional, default: []) - RxC boolean matrix, labels where mask
%    is false are set to 0
%
%    removeBorder (optional, default: false) - boolean value, if set to
%    true, segments that touch the border of the CHM are removed
%
%    verbose (optional, default: true) - boolean value, verbosity switch
%
%    fig (optional, default: true) - boolean value, switch to plot figures
%
% Outputs:
%    label - RxC numeric matrix, individual tree crown labels
%    colors - Nx1 numeric matrix, topological color indices
%
% Example:
%    [chm, refmat, ~] = geotiffread('..\data\measurements\raster\chm\so_2014_woodland_pasture.tif');
%
%    [crh, xyh] = canopyPeaks(double(chm), ...
%        refmat, ...
%        'minTreeHeight', 2, ...
%        'searchRadius', @(h) 1, ...
%        'mergeRadius', @(h) 0.28 * h^0.59, ...
%        'fig', true, ...
%        'verbose', true);
%
%    label = treeWatershed(chm, ...
%        'markers', crh, ...
%        'minHeight', 1, ...
%        'fig', true, ...
%        'verbose', true);
%
% Other m-files required: topoColor.m
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2017b, GNU Octave 4.4.1 (configured for "x86_64-w64-mingw32")
%
% See also: rasterize.m
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory (LASIG)
% Website: http://mparkan.github.io/Digital-Forestry-Toolbox/
% Last revision: January 23, 2019
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund (WHFF, OFEV), project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'chm', @isnumeric);
addParameter(arg, 'markers', [], @isnumeric);
addParameter(arg, 'minHeight', 1, @(x) isnumeric(x) && (numel(x) == 1));
addParameter(arg, 'mask', [], @(x) islogical(x) && any(x(:)) && all(size(x) == size(chm)));
addParameter(arg, 'removeBorder', false, @(x) islogical(x) && (numel(x) == 1));
addParameter(arg, 'fig', true, @(x) islogical(x) && (numel(x) == 1));
addParameter(arg, 'verbose', true, @(x) islogical(x) && (numel(x) == 1));

parse(arg, chm, varargin{:});
nargoutchk(1, 2);


%% apply height threshold

chm = double(chm);
chm(chm <= arg.Results.minHeight) = 0;


%% compute watershed transform

if arg.Results.verbose
    
    fprintf('computing watershed transform...');
    pause(0.01);
    
end

[nrows, ncols] = size(chm);
I = -chm;

% impose regional minima in CHM (marker controled watershed)
if ~isempty(arg.Results.markers)
    
    ind_peaks = sub2ind(size(chm), arg.Results.markers(:,2), arg.Results.markers(:,1));
    marker = false(size(chm));
    marker(ind_peaks) = 1;

    I = imimposemin(I, marker);

end

% apply watershed algorithm
label = single(watershed(I, 8));
label(chm == 0) = 0;

% relabel connected components
label = bwlabel(label ~= 0, 8);
label(~ismember(label, label(marker))) = 0;
% label(label ~= 0) = grp2idx(label(label ~= 0));
[label(label ~= 0), ~] = grp2idx(label(label ~= 0));

if arg.Results.verbose
    
    fprintf('done!\n');
    
end


%% remove segments that touch the border of the CHM

if arg.Results.removeBorder
    
    B = padarray(false(nrows-2, ncols-2), [1 1], true);
    label(ismember(label, label(B))) = 0;
    
end


%% apply mask

if ~isempty(arg.Results.mask)

    label(~arg.Results.mask) = 0;
    
end


%% remove nan values

label(isnan(chm)) = 0;

% relabel connected components
label = bwlabel(label ~= 0, 8);


%% remove seam lines

[d_nn, idxn_nn] = bwdist(label ~= 0);
idxl_dist = ((d_nn >= 1) & (d_nn < 2));
label(idxl_dist) = label(idxn_nn(idxl_dist));


%% reassign labels

[label(label ~= 0), ~] = grp2idx(label(label ~= 0));
varargout{1} = label;


%% plot results

if arg.Results.fig || (nargout == 2)
    
    if arg.Results.verbose
        
        fprintf('topological coloring...');
        pause(0.01);
        
    end
    
    [nrows, ncols] = size(label);
    [cols, rows] = meshgrid(1:ncols, 1:nrows);
    
    [idxn_color, ~, cmap] = topoColor([cols(:), rows(:)], ...
        label(:), ...
        'adjacency', '2d', ...
        'buffer', 3, ...
        'colormap', 'cmap12', ...
        'unlabelledColor', [0, 0, 0], ...
        'fig', false, ...
        'verbose', false);
    
    colors = reshape(idxn_color, nrows, ncols);
    
    if nargout == 2
        
        varargout{2} = colors;
        
    end
    
    if arg.Results.verbose
        
        fprintf('done!\n');
        
    end
    
end

if arg.Results.fig
    
    figure
    imagesc(colors, [0, size(cmap,1)])
    axis equal tight
    xlabel('col');
    ylabel('row');
    colormap([0, 0, 0; cmap])
    
end
