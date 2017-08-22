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
%    adjacency (optional, default: @(h) min(0.5 + 0.5*log(max(h,1)),4) - anonymous function handle of the form @(h) = ...,
%    which defines the minimum 3D distance separating peaks. Peaks separated by a smaller distance
%    are grouped and only the highest peak within each group is retained.
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
%    
%    [chm, refmat, ~] = geotiffread('..\data\measurements\raster\chm\so_2014_woodland_pasture.tif');
%    
%    [crh, xyh] = canopyPeaks(double(chm), ...
%        refmat, ...
%        'method', 'allometricRadius', ...
%        'allometry', @(h) 1 + 0.5*log(max(h,1)), ...
%        'adjacency' @(h) min(0.5 + 0.5*log(max(h,1)),4), ...
%        'fig', true, ...
%        'verbose', true);
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
% Last revision: August 22, 2017
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund (WHFF, OFEV), project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'chm', @isnumeric);
addRequired(arg, 'refmat', @isnumeric);
addParameter(arg, 'method', 'fixedRadius', @(x) any(validatestring(x, {'fixedRadius', 'allometricRadius', 'hMaxima'})));
addParameter(arg, 'minPeakHeight', 2, @(x) isnumeric(x) && (numel(x) == 1));
addParameter(arg, 'windowRadius', 4, @(x) isnumeric(x) && (numel(x) == 1));
addParameter(arg, 'allometry', @(h) 1 + 0.5*log(max(h,1)), @(x) strfind(func2str(x),'@(h)') == 1);
addParameter(arg, 'adjacency', @(h) min(0.5 + 0.5*log(max(h,1)),4), @(x) strfind(func2str(x),'@(h)') == 1);
addParameter(arg, 'minHeightDifference', 0.1, @isnumeric);
addParameter(arg, 'fig', true, @(x) islogical(x) && (numel(x) == 1));
addParameter(arg, 'verbose', false, @(x) islogical(x) && (numel(x) == 1));

parse(arg, chm, refmat, varargin{:});


%% find tree tops

if arg.Results.verbose
   
    fprintf('detecting peaks...');
    
end

gridResolution = abs(refmat(1,2));
chm(chm < 0) = 0;
chm = double(chm);

switch arg.Results.method
    
    case 'fixedRadius'
        
        SE = strel('disk', round(arg.Results.windowRadius / gridResolution)).getnhood();
        SE(ceil(size(SE,1)/2), ceil(size(SE,2)/2)) = 0; % set central convolution window value to zero
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

% transform image to map coordinates
% crh = [col_lm, row_lm, val_lm];
[nrows, ncols] = size(chm);
xy = [row_lm, col_lm, ones(size(row_lm))] * refmat;
xyh = [xy, val_lm];

if arg.Results.verbose
    
    fprintf('done!\n');
    
end

%% filter tree tops below height threshold

if arg.Results.verbose
   
    fprintf('filtering peaks...');
    
end

idxl_height_filter = (val_lm >= arg.Results.minPeakHeight);

% crh = crh(idxl_height_filter,:);
xyh = xyh(idxl_height_filter,:);

% sort peaks by decreasing height
[~, idxn_sort] = sort(xyh(:,3), 1, 'descend');
xyh = xyh(idxn_sort,:);

if arg.Results.verbose
    
    fprintf('done!\n');
    
end


%% merge adjacent treetops 

if arg.Results.verbose
   
    fprintf('merging peaks...');
    
end

r_adj = arg.Results.adjacency(xyh(:,3));
r_max = max(r_adj);

[idxn_adj, d_adj] = rangesearch(xyh(:,1:3), ...
    xyh(:,1:3), ...
    r_max);

for j = 1:length(idxn_adj)
    
    idxl_range = d_adj{j} <= r_adj(j);
    idxn_adj{j} = idxn_adj{j}(idxl_range);
    d_adj{j} = d_adj{j}(idxl_range);
    
end

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

% compute highest peak in cluster
idxn_peak = accumarray(idxn_cluster, G.Nodes.ID, [], @(x) x(1));
xyz_cl = xyh(idxn_peak,:); % xyh_peak
xyh = xyz_cl;
    
% transform map to image coordinates
rc_cl = [xyz_cl(:,1) - refmat(3,1), xyz_cl(:,2) - refmat(3,2)] / refmat(1:2,:);
rc_cl = round(rc_cl);
rc_cl(rc_cl == 0) = 1;
rc_cl(rc_cl(:,1) > nrows,1) = nrows;
rc_cl(rc_cl(:,2) > ncols,2) = nrows;
crh_cl = [rc_cl(:,[2 1]), xyz_cl(:,3)];
crh = crh_cl;

if arg.Results.verbose
    
    fprintf('done!\n');
    
end


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