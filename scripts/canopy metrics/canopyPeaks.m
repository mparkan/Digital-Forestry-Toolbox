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
%    allometry (optional, default: @(h) 0.5 + 0.25*log(max(h,1))) - anonymous function handle of the form @(h) = ...,
%    specifying the allometric relation between tree height and crown
%    radius in map units when method is set to 'allometric radius'.
%     examples:
%     @(h) (3.09632 + 0.00895 * h.^2)/2; % deciduous forest (Popescu et al, 2004)
%     @(h) (3.75105 - 0.17919 * h + 0.01241 * h.^2)/2; % coniferous forests (Popescu et al, 2004)
%     @(h) (2.51503 + 0.00901 * h.^2)/2; % mixed forests (Popescu et al, 2004)
%     @(h) (1.7425 * h.^0.5566)/2; % mixed forests (Chen et al., 2006)
%     @(h) (1.2 + 0.16 * h)/2; % mixed forests (Pitkänen et al., 2004)
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
%        'allometry', @(h) 0.5 + 0.25*log(max(h,1)), ...
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
% Last revision: March 8, 2018
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
addParameter(arg, 'minHeightDifference', 0.1, @isnumeric);
addParameter(arg, 'fig', true, @(x) islogical(x) && (numel(x) == 1));
addParameter(arg, 'verbose', false, @(x) islogical(x) && (numel(x) == 1));

parse(arg, chm, refmat, varargin{:});

OCTAVE_FLAG = (exist('OCTAVE_VERSION', 'builtin') ~= 0); % determine if system is Matlab or GNU Octave


%% find tree tops

if arg.Results.verbose
   
    fprintf('detecting peaks...');
    
end

gridResolution = abs(refmat(1,2));
chm(chm < 0) = 0;
chm(isnan(chm)) = 0;
chm = double(chm);

switch arg.Results.method
    
    case 'fixedRadius'
        
        r = round(arg.Results.windowRadius / gridResolution);
        SE = bwdist(padarray(true, [r,r])) <= r;
        SE(ceil(size(SE,1)/2), ceil(size(SE,2)/2)) = 0; % set central convolution window value to zero
        idxl_lm = chm > imdilate(chm, SE);
        val_lm = chm(idxl_lm);
        [row_lm, col_lm] = ind2sub(size(idxl_lm), find(idxl_lm));
        
    case 'allometricRadius'
        
        % determine convolution window size for each pixel
        crown_radius = arg.Results.allometry(chm);
        window_radius = max(round(crown_radius ./ gridResolution),1);
        unique_window_radius = unique(window_radius);
        n_radius = length(unique_window_radius);
        
        idxn_lm = cell(n_radius,1);
   
        for j = 1:n_radius
            
            r = unique_window_radius(j);
            SE = bwdist(padarray(true, [r,r])) <= r;
            SE(ceil(size(SE,1)/2), ceil(size(SE,2)/2)) = 0;
            idxn_lm{j,1} = find((chm > imdilate(chm, SE)) & (window_radius == r));
            
        end
        
        [row_lm, col_lm] = ind2sub(size(window_radius), cell2mat(idxn_lm));
        val_lm = chm(sub2ind(size(chm), row_lm, col_lm));
        
    case 'hMaxima'    
        
        % H-maxima transformation (remove local maxima with height < h)
        I_hmax = imreconstruct(chm-arg.Results.minHeightDifference, chm, 8);
        idxl_lm = imregionalmax(I_hmax, 8);
        
        % find centroids of maxima regions
        CC = bwconncomp(idxl_lm, 8);
        stats = regionprops(CC, 'centroid');
        cr_centroid = ceil(cell2mat({stats.Centroid}'));
        row_lm = cr_centroid(:,2);
        col_lm = cr_centroid(:,1);
        val_lm = chm(sub2ind(size(chm), row_lm, col_lm));

end

% transform image to map coordinates
xy = [row_lm, col_lm, ones(size(row_lm))] * refmat;
crh = [col_lm, row_lm, val_lm];
xyh = [xy, val_lm];

if arg.Results.verbose
    
    fprintf('done!\n');
    
end


%% filter tree tops below height threshold

if arg.Results.verbose
   
    fprintf('filtering peaks...');
    
end

idxl_height_filter = (val_lm >= arg.Results.minPeakHeight);

crh = crh(idxl_height_filter,:);
xyh = xyh(idxl_height_filter,:);

% sort peaks by decreasing height
[~, idxn_sort] = sort(xyh(:,3), 1, 'descend');
xyh = xyh(idxn_sort,:);
crh = crh(idxn_sort,:);

if arg.Results.verbose
    
    fprintf('done!\n');
    
end


%% plot individual tree tops

if arg.Results.fig
    
    figure
    imagesc(chm);
    colormap('gray');
    hold on
    plot(crh(:,1), crh(:,2), 'rx', 'MarkerSize', 3);
    axis equal tight
    caxis(quantile(chm(:), [0.01, 0.99]))
    xlabel('col');
    ylabel('row');
    
end