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
%    minHeight (optional, default: 1) - numeric value, minimum canopy height, all values
%    below this threshold are set to zero.
%
%    mask (optional, default: []) - RxC boolean matrix, any segment that is
%    not located entirely within the mask is excluded (i.e. label set to
%    0). By default all segments are included.
%
%    verbose (optional, default: true) - boolean value, verbosiy switch
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
%        'method', 'allometricRadius', ...
%        'allometry', @(h) 0.5 + 0.25*log(max(h,1)), ....
%        'fig', true, ...
%        'verbose', true);
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
% Compatibility: tested on Matlab R2017b, GNU Octave 4.2.1 (configured for "x86_64-w64-mingw32")
%
% See also: rasterize.m
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory (LASIG)
% Website: http://mparkan.github.io/Digital-Forestry-Toolbox/
% Last revision: March 9, 2018
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund (WHFF, OFEV), project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'chm', @isnumeric);
addParameter(arg, 'markers', [], @isnumeric);
addParameter(arg, 'minHeight', 1, @(x) isnumeric(x) && (numel(x) == 1));
addParameter(arg, 'mask', [], @(x) islogical(x) && any(x(:)) && all(size(x) == size(chm)));
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

I = -chm;

% impose regional minima in CHM (marker controled watershed)
if ~isempty(arg.Results.markers)
    
    ind_peaks = sub2ind(size(chm), arg.Results.markers(:,2), arg.Results.markers(:,1));
    marker = false(size(chm));
    marker(ind_peaks) = 1;
    
    I_range = max(I(:)) - min(I(:));
    
    if I_range == 0
        
        h0 = 0.1;
        
    else
        
        h0 = 0.0001 * I_range;
        
    end
    
    E = inf(size(chm));
    E(marker) = -Inf;
    
    J = imreconstruct(imcomplement(E), imcomplement(min(I + h0, E)), 8);
    I = imcomplement(J);
    
end

% apply watershed algorithm
label = single(watershed(I, 8));

if arg.Results.verbose
    
    fprintf('done!\n');
    
end


%% apply mask

label(chm < 1) = 0;

%% remove segments that intersect the inverted mask

if ~isempty(arg.Results.mask)

    idxl_mask = ismember(label, label(~arg.Results.mask));
    label(idxl_mask) = 0;
    
end

%% remove seam lines

[d_nn, idxn_nn] = bwdist(label ~= 0);
idxl_dist = ((d_nn >= 1) & (d_nn < 2));
label(idxl_dist) = label(idxn_nn(idxl_dist));


%% remove nan values

label(isnan(chm)) = 0;


%% reassign labels

[label(label ~= 0), ~] = grp2idx(label(label ~= 0));
varargout{1} = label;


%% plot results

if arg.Results.fig || (nargout == 2)
    
    if arg.Results.verbose
        
        fprintf('topological coloring...');
        pause(0.01);
        
    end
    
    % topological coloring
    cmap = [0, 0, 0;
        166,206,227;
        31,120,180;
        178,223,138;
        51,160,44;
        251,154,153;
        227,26,28;
        253,191,111;
        255,127,0;
        202,178,214;
        106,61,154;
        255,255,153;
        177,89,40] ./ 255;
    
    N = max(label(:));
    colors = zeros(size(label));
    n_colors = 100;
    
    for k = 1:N
        
        idxn_adj = imdilate(label == k, true(4));
        
        % list available colors
        idxl_color_pool = ~ismember(1:n_colors, colors(idxn_adj));
        
        % assign first available color
        colors(idxn_adj) = find(idxl_color_pool,1);
        
    end
    
    colors = min(colors, size(cmap,1));
    colors(label == 0) = 0;
    
    if nargout == 2
        
        varargout{2} = colors;
        
    end
    
    if arg.Results.verbose
        
        fprintf('done!\n');
        
    end
    
end

if arg.Results.fig
    
    figure
    imagesc(colors+1)
    axis equal tight
    xlabel('col');
    ylabel('row');
    colormap(cmap);
    
end
