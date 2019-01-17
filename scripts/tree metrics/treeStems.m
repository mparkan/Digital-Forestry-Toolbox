function [label, varargout] = treeStems(xyh, varargin)
% TREESTEMS - detects stems by applying a local maxima filter on a vertical point density raster. 
%
% References:
%
% [1] Rahman, M.Z.A., Gorte, B., 2008. Individual tree detection based on densities of high points of 
% high resolution airborne LiDAR. GEOBIA 350–355.
%
% [2] Ayrey, E., Fraver, S., Jr, J.A.K., Kenefic, L.S., Hayes, D., Weiskittel, A.R., Roth, B.E., 2017. 
% Layer Stacking: A Novel Algorithm for Individual Forest Tree Segmentation from LiDAR Point Clouds. 
% Canadian Journal of Remote Sensing 43, 16–27. https://doi.org/10.1080/07038992.2017.1252907
% 
% Syntax: [label, ...] = treeStems(xyh, ...)
%
% Inputs:
%    xyh - Nx3 numeric matrix, 3D point coordinates (normalized relative to the terrain elevation)
%
%    filter (optional, default: []) - Nx1 boolean vector, logical index indicating the subset of points to be used in the stem detection
%
%    cellSize (optional, default: 0.4) - numeric value, raster cell size
%
%    bandWidth (optional, default: 1.5) - numeric value, raster layer thickness
%
%    verticalStep (optional, default: 0.25) - numeric value, vertical shift between successive raster layers
%
%    searchRadius (optional, default: 2) - numeric value, radius in map units used for the local density maxima search
%
%    minLength (optional, default: 5) - numeric value, minimum layer accumulation length to detect a stem
%
%    verbose (optional, default: true) - boolean value, verbosity switch
%
%    fig (optional, default: false) - boolean value, switch to plot figures
%
% Outputs:
%    label - Nx1 numeric vector, label vector identifying individual stems in the input point cloud
%
%    xyh_stems (optional) - numeric matrix, map coordinates and height (x, y, h) of detected stems
%
% Example:
%
% xyh = alt2height([pc.record.x, pc.record.y, pc.record.z], ...
%     pc.record.classification, ...
%     'classTerrain', [2], ...
%     'method', 'linear');
% 
% idxl_last = pc.record.return_number == pc.record.number_of_returns;
% idxl_veg = ismember(pc.record.classification, [4, 5]);
% idxl_filter = idxl_veg & idxl_last;
% 
% [label_s, xyh_stem] = treeStems(xyh, ...
%     idxl_filter, ...
%     'cellSize', 0.3, ...
%     'bandWidth', 1, ...
%     'verticalStep', 0.5, ...
%     'searchRadius', 2, ...
%     'minLength', 5, ...
%     'fig', true);
%
% Other m-files required: rasterize.m
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
% Last revision: January 17, 2019
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund (WHFF, OFEV), project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'xyh', @(x) (size(x,2) == 3) && isnumeric(x));
addOptional(arg, 'filter', true(size(xyh,1),1), @(x) (size(x,1) == size(xyh,1)) && islogical(x) && any(x));
addParameter(arg, 'cellSize', 0.3, @(x) isnumeric(x) & (numel(x) == 1));
addParameter(arg, 'bandWidth', 1, @(x) isnumeric(x) & (numel(x) == 1) & (x > 0));
addParameter(arg, 'verticalStep', 0.4, @(x) isnumeric(x) & (numel(x) == 1) & (x > 0));
addParameter(arg, 'searchRadius', 2, @(x) isnumeric(x) & (numel(x) == 1) & (x > 0));
addParameter(arg, 'minLength', 5, @(x) isnumeric(x) & (numel(x) == 1) & (x > 0));
addParameter(arg, 'fig', false, @(x) islogical(x) && (numel(x) == 1));
addParameter(arg, 'verbose', true, @(x) islogical(x) && (numel(x) == 1));

parse(arg, xyh, varargin{:});

% check output argument format
nargoutchk(1, 2);


%% rasterize

if arg.Results.verbose
    
    tic
    fprintf('rasterizing point cloud...');
    
end

xyh = double(xyh);
xyh(xyh(:,3) < 0,3) = 0;

[I, refmat, ind] = rasterize(xyh(:,1:2), ...
    xyh(:,3), ...
    'CellSize', arg.Results.cellSize, ...
    'fun', @(x) numel(x), ...
    'fill', nan, ...
    'fig', false);

[nrows, ncols] = size(I);

if arg.Results.verbose
    
    fprintf('done!\n');
    toc
    
end


%% apply filter

xyh_s = xyh(arg.Results.filter,:);
ind_s = ind(arg.Results.filter);


%% compute canopy height

H = accumarray(ind, xyh(:,3), [nrows*ncols,1], @max, 0);
H = reshape(H, nrows, ncols);


%% extract layers

z_min = min(xyh_s(:,3));
z_max = max(xyh_s(:,3));
z_sample = z_min:arg.Results.verticalStep:z_max;
n_layers = length(z_sample);

A = zeros(nrows, ncols, 'single');

reverse_str = '';

for j = 1:n_layers
    
    idxl_layer = (xyh_s(:,3) >= z_sample(j) - arg.Results.bandWidth/2) & (xyh_s(:,3) <= z_sample(j) + arg.Results.bandWidth/2);
    
    % rasterize layer
    BW = false([nrows, ncols]);
    BW(ind_s(idxl_layer)) = true;
    
    % compute distance transform
    [D, ~] = bwdist(BW);
    M = (D <= (0.5/arg.Results.cellSize));
    
    A = A + single(M);
    
    if arg.Results.verbose
        
        msg = sprintf('scanning layer %u/%u...', j, n_layers);
        fprintf([reverse_str, msg]);
        reverse_str = repmat(sprintf('\b'), 1, length(msg));
        
    end
    
end

if arg.Results.verbose
    
    fprintf('done!\n');
    toc
    
end


%% detect local density maxima

r = max(1, ceil(arg.Results.searchRadius / arg.Results.cellSize));
SE = bwdist(padarray(true, [r,r])) <= r;
SE(ceil(size(SE,1)/2), ceil(size(SE,2)/2)) = 0;
idxl_lm = (A >= imdilate(A, SE));
idxl_lm(A <= (arg.Results.minLength/arg.Results.verticalStep)) = false;

% label connected components
L = bwlabel(idxl_lm, 8);


%% compute stem positions

if nargout == 2
    
    if arg.Results.verbose
        
        tic
        fprintf('computing stem positions...');
        
    end
    
    label = L(ind);
    idxl_lab = (label ~= 0);
    
    label(idxl_lab) = grp2idx(label(idxl_lab));
    
    N = length(unique(label(idxl_lab)));

    xyh_proxy = zeros(N,3);
    xyh_proxy(:,1) = accumarray(label(idxl_lab), xyh(idxl_lab,1), [N, 1], @mean, nan);
    xyh_proxy(:,2) = accumarray(label(idxl_lab), xyh(idxl_lab,2), [N, 1], @mean, nan);
    
    rc_proxy = round([xyh_proxy(:,1) - refmat(3,1), xyh_proxy(:,2) - refmat(3,2)] / refmat(1:2,:));
    ind_proxy = sub2ind([nrows, ncols], rc_proxy(:,1), rc_proxy(:,2));
    
    xyh_proxy(:,3) = H(ind_proxy);
    
    varargout{1} = xyh_proxy;
    
    if arg.Results.verbose
        
        fprintf('done!\n');
        toc
        
    end
    
end


%% plot

if arg.Results.fig
     
    figure
    imagesc(H);
    colormap('gray');
    hold on
    plot(rc_proxy(:,2), rc_proxy(:,1), 'rx', 'MarkerSize', 3);
    axis equal tight
    % caxis(quantile(chm(:), [0.01, 0.99]))
    colorbar
    xlabel('col');
    ylabel('row');

end