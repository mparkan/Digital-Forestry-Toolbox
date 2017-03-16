function [crh, xyh] = canopyPeaks(chm, refmat, varargin)
% CANOPYPEAKS - find local maxima coordinates in a raster Canopy Height Model (CHM).
%
% [CRH, XYH] = CANOPYPEAKS(CHM, REFMAT, METHOD, ...) computes candidate
% tree top coordinates in the image (CRH) and map (XYH) reference systems
% using the specified METHOD. Three methods are implemented:
% 1. Fixed size circular window convolution. The optional parameter 'windowRadius' can be 
% specified to set the window radius.
% 2. Variable size circular window (based on allometry) convolution, see
% Popescu et al. (2004) [1] and Chen et al. (2006) [2].
% 3. H-maxima transform, see Kwak et al. (2007) [3]
%
% [1] S. C. Popescu and R. H. Wynne, "Seeing the trees in the forest: Using lidar and 
% multispectral data fusion with local filtering and variable window size for estimating tree height,"
% Photogrammetric Engineering and Remote Sensing, vol. 70, no. 5, pp. 589–604, 2004.
%
% [2] Q. Chen, D. Baldocchi, P. Gong, and M. Kelly, "Isolating Individual Trees in a 
% Savanna Woodland Using Small Footprint Lidar Data," Photogrammetric Engineering & 
% Remote Sensing, vol. 72, no. 8, pp. 923–932, Aug. 2006.
%
% [3] D.-A. Kwak, W.-K. Lee, J.-H. Lee, G. S. Biging, and P. Gong, 
% "Detection of individual trees and estimation of tree height using LiDAR data,"
% J For Res, vol. 12, no. 6, pp. 425–434, Oct. 2007.
%
% Syntax:  [crh, xyh] = canopyPeaks(chm, refmat, method, ...)
%
% Inputs:
%    chm - RxC numeric matrix, raster Canopy Height Model (CHM) with R rows
%    and C columns
%
%    refmat - 3x2 numeric matrix, spatial referencing matrix such that [map_x map_y] = [row col 1] * refmat
%
%    adjacencyRange (optional, default: 1.5) - numeric value, minimum
%    3D distance separating peaks. Peaks separated by a smaller distance
%    are merged and the new peaks are the centroids of all adjacent peaks.
%
%    method (optional, default: 'fixedRadius') - string, peak detection method: 'fixedRadius', 'allometricRadius' or 'hMaxima'.
%
%    minPeakHeight (optional, default: 2) - numeric value, minimum tree top height
%
%    windowRadius (optional, default: 4) - numeric value, fixed circular window radius in map units
%    used to detect local maxima when method is set to 'fixedRadius'.
%
%    allometry (optional, default: @(h) 1 + 0.5*log(max(h,1))) - anonymous function handle of the form @(h) = ...,
%    specifying the allometric relation between tree height and crown
%    diameter (not radius) in map units when method is set to 'allometric radius'.
%     examples:
%     @(h) 3.09632 + 0.00895 * h.^2; % deciduous forest (Popescu et al, 2004)
%     @(h) 3.75105 - 0.17919 * h + 0.01241 * h.^2; % coniferous forests (Popescu et al, 2004)
%     @(h) 2.51503 + 0.00901 * h.^2; % mixed forests (Popescu et al, 2004)
%     @(h) 1.7425 * h.^0.5566; % mixed forests (Chen et al., 2006)
%     @(h) 1.2 + 0.16 * h; % mixed forests (Pitkänen et al., 2004)
%
%    minHeightDifference (optional, default: 0.1) - numeric value, threshold
%    height difference below which the H-maxima transform suppresses all local maxima (only when method is set to 'hMaxima')
%
%    fig (optional, default: false) - boolean value, switch to plot figures
%
% Outputs:
%    crh - Mx3 numeric matrix, images coordinates (col, row) and height values of tree tops
%
%    xyh - Mx3 numeric matrix, map coordinates (x, y) and height values of tree tops
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
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2016b
%
% See also: treeWatershed.m
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory (LASIG)
% Website: http://mparkan.github.io/Digital-Forestry-Toolbox/
% Last revision: March 15, 2017
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund (WHFF, OFEV), project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'chm', @isnumeric);
addRequired(arg, 'refmat', @isnumeric);
addParameter(arg, 'method', 'fixedRadius', @(x) any(validatestring(x, {'fixedRadius', 'allometricRadius', 'hMaxima'})));
addParameter(arg, 'adjacencyRange', 1.5, @(x) isnumeric(x) && x > 0 && (numel(x) == 1));
addParameter(arg, 'minPeakHeight', 2, @(x) isnumeric(x) && (numel(x) == 1));
addParameter(arg, 'windowRadius', 4, @(x) isnumeric(x) && (numel(x) == 1));
addParameter(arg, 'allometry', @(h) 1 + 0.5*log(max(h,1)), @(x) strfind(func2str(x),'@(h)') == 1);
addParameter(arg, 'minHeightDifference', 0.1, @isnumeric);
addParameter(arg, 'fig', true, @(x) islogical(x) && (numel(x) == 1));

parse(arg, chm, refmat, varargin{:});


%% find tree tops

gridResolution = abs(refmat(1,2));

switch arg.Results.method
    
    case 'fixedRadius'
        
        SE = strel('disk', round(arg.Results.windowRadius / gridResolution)).getnhood();
        SE(ceil(size(SE,1)/2), ceil(size(SE,2)/2)) = 0;
        idx_lm = chm > imdilate(chm, SE);
        val_lm = chm(idx_lm);
        [row_lm, col_lm] = ind2sub(size(idx_lm), find(idx_lm));
        
    case 'allometricRadius'
        
        % determine convolution window size for each pixel
        crown_radius = arg.Results.allometry(chm) ./ 2;
        window_radius = max(round(crown_radius ./ gridResolution),1);
        unique_window_radius = unique(window_radius);

        n_radius = length(unique_window_radius);
        
        idx_lm = cell(n_radius,1);
   
        for j = 1:n_radius
            
            SE = strel('disk', unique_window_radius(j)).getnhood();
            SE(ceil(size(SE,1)/2), ceil(size(SE,2)/2)) = 0;
            
            idx_lm{j,1} = find((chm > imdilate(chm, SE)) & (window_radius == unique_window_radius(j)));
            
        end
        
        [row_lm, col_lm] = ind2sub(size(window_radius), cell2mat(idx_lm));
        val_lm = chm(sub2ind(size(chm), row_lm, col_lm));
        
    case 'hMaxima'    
        
        % The H-maxima transform suppresses all maxima in the intensity image I whose height is less than h
        idx_lm = imextendedmax(chm, arg.Results.minHeightDifference, 8);
        val_lm = chm(idx_lm);
        [row_lm, col_lm] = ind2sub(size(idx_lm), find(idx_lm));
        
end

% crh = [col_lm, row_lm, val_lm];
[nrows, ncols] = size(chm);

t = [row_lm, col_lm] * refmat(1:2,:);
x = t(:,1) + refmat(3,1);
y = t(:,2) + refmat(3,2);
h = val_lm;
xyh = [x,y,h];


%% filter tree tops below height threshold

idxl_height_filter = (val_lm >= arg.Results.minPeakHeight);

% crh = crh(idxl_height_filter,:);
xyh = xyh(idxl_height_filter,:);


%% merge adjacent treetops 

[idxn_adj, d_adj] = rangesearch(xyh(:,1:3), ...
    xyh(:,1:3), ...
    arg.Results.adjacencyRange);

n_adj = cellfun(@numel, idxn_adj);
idxl_pair = (n_adj > 1);

pairs = cell2mat(cellfun(@(x) nchoosek(x,2), idxn_adj(idxl_pair), 'UniformOutput', false));

% remove duplicate edges
[~, idxl_unique_edges, ~] = unique(sort(pairs, 2), 'rows');

pairs = pairs(idxl_unique_edges,:);

% build  graph
EdgeTable = table(pairs, ...
    'VariableNames',{'EndNodes'});

NodeTable = table((1:size(xyh,1))', ...
    'VariableNames',{'ID'});

G = graph(EdgeTable, NodeTable);

% find connected components
idxn_cluster = conncomp(G, ...
    'OutputForm', 'vector')';

% compute cluster centroids
xyz_cl(:,1) = accumarray(idxn_cluster, xyh(:,1), [], @mean);
xyz_cl(:,2) = accumarray(idxn_cluster, xyh(:,2), [], @mean);
xyz_cl(:,3) = accumarray(idxn_cluster, xyh(:,3), [], @mean);

rc_cl = [xyz_cl(:,1) - refmat(3,1), xyz_cl(:,2) - refmat(3,2)] / refmat(1:2,:);
rc_cl = round(rc_cl);
rc_cl(rc_cl == 0) = 1;
rc_cl(rc_cl(:,1) > nrows,1) = nrows;
rc_cl(rc_cl(:,2) > ncols,2) = nrows;
crh_cl = [rc_cl(:,[2 1]), xyz_cl(:,3)];
crh = crh_cl;
xyh = xyz_cl;


%% plot individual tree tops

if arg.Results.fig
    
    figure
    imagesc(chm);
    colormap('gray');
    colorbar;
    hold on
    plot(crh(:,1), crh(:,2), 'rx');
    axis equal tight
    caxis(quantile(chm(:), [0.01, 0.99]))
    xlabel('col');
    ylabel('row');
    
end