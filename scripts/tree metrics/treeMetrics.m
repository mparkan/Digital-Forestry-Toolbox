function [metrics, varargout] = treeMetrics(label, xyz, intensity, returnNumber, returnTotal, rgb, dtm, refmat, varargin)
% TREEMETRICS - computes segment metrics
%
% [LABEL, ...] = TREEMETRICS(LABEL, XYZ, INTENSITY, RETURNNUMBER, RETURNTOTAL, RGB, DTM, REFMAT ...) computes various
% segment metrics based on geometry, intensity, opacity (pulse return) and color characteristics
%
% Syntax:  metrics = treeMetrics(label, xyz, intensity, returnNumber, returnTotal, rgb, dtm, refmat ...)
%
% Inputs:
%    label - Nx1 integer vector, point label (individual tree label)
%
%    xyz - Nx3 numeric matrix, point cloud coordinates [x y z]
%
%    intensity - Nx1 numeric vector, intensity
%
%    returnNumber - Nx1 numeric vector, return number
%
%    returnTotal - Nx1 numeric vector, number of returns per pulse
%
%    rgb - Nx3 numeric vector, [red, green, blue] triplets. Set to
%    zeros(N,3) if the point cloud has no color.
%
%    dtm - numeric matrix, terrain elevation model
%
%    refmat - 3x2 numeric matrix, spatial referencing matrix of the terrain model, such that xy_map = [row, col, ones(nrows,1)] * refmat
%
%    metrics (optional, default: {'all'}) - cell array of strings, list of
%    metrics to compute for each segment. The list can contain a
%    combination of individual features (e.g. 'UUID', 'ConvexVolume',
%    'TotalHeight', 'IntensityMedian') or feature categories. Supported
%    feature categories are: 'Identifier', 'Basic', 'PointPatternMetrics',
%    'IntensityMetrics', 'OpacityMetrics', 'ColorMetrics',
%    'ExternalShapeMetrics'. If you want to compute all metrics specify
%    {'all'}.
%
%    treePos (optional, default: []) - Mx3 numeric matrix, xyz coordinates of
%    the tree positions. M is the number of segments (trees). If not specified,
%    the tree positions are automatically computed using the tree top
%    (highest point in segment) as a proxy for tree position
%
%    intensityScaling (optional, default: true) - boolean value, rescale
%    intensity between 0 and 1 using the [0.001, 0.99] quantiles as limits
%
%    scalarOnly (optional, default: true) - boolean value, flag indicating
%    if only scalar values should be returned
%
%    verbose (optional, default: true) - boolean value, verbosiy switch
%
% Outputs:
%    metrics - Mx1 structure, list metrics for each segment.
%
%    scalar - 1xF boolean vector, flag to indicate scalar fields. F is the number of output metrics.
%
%    format - 1xF cell array, print format of scalar fields (non-scalar fields have an empty print format).
%    F is the number of output metrics.
%
% Example:
%
%    metrics = treeMetrics(label_3d, ...
%         [pc.record.x, pc.record.y, pc.record.z], ...
%         pc.record.intensity, ...
%         pc.record.return_number, ...
%         pc.record.number_of_returns, ...
%         [pc.record.red, pc.record.green, pc.record.blue], ...
%         models.terrain.values, ...
%         refmat, ...
%         'metrics', {'Identifier', 'IntensityMetrics', 'OpacityMetrics', 'ConvexVolume', 'ConvexArea'}, ...
%         'intensityScaling', true, ...
%         'verbose', true);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2019b, GNU Octave 5.2.0 (configured for "x86_64-w64-mingw32")
%
% See also:
%
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory (LASIG)
% Website: http://mparkan.github.io/Digital-Forestry-Toolbox/
% Last revision: February 23, 2020
% Acknowledgments: This work was supported by the Swiss Forestry and Wood
% Research Fund, WHFF (OFEV) - project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'label', @(x) (size(x,2) == 1) && isnumeric(x) && any(x ~= 0));
addRequired(arg, 'xyz', @(x) (size(x,2) == 3) && isnumeric(x));
addRequired(arg, 'intensity', @(x) (size(x,2) == 1) && isnumeric(x));
addRequired(arg, 'returnNumber', @(x) (size(x,2) == 1) && isnumeric(x));
addRequired(arg, 'returnTotal', @(x) (size(x,2) == 1) && isnumeric(x));
addRequired(arg, 'rgb', @(x) (size(x,2) == 3) && isnumeric(x));
addRequired(arg, 'dtm', @(x) isnumeric(x));
addRequired(arg, 'refmat', @(x) all(size(x) == [3,2]));
addParameter(arg, 'metrics', {'all'}, @(x) iscell(x) && ~isempty(x));
addParameter(arg, 'treePos', [], @(x) (size(x,2) == 3) && isnumeric(x));
addParameter(arg, 'intensityScaling', true, @(x) islogical(x) && (numel(x) == 1));
addParameter(arg, 'verbose', true, @(x) islogical(x) && (numel(x) == 1));

parse(arg, label, xyz, intensity, returnNumber, returnTotal, rgb, dtm, refmat, varargin{:});

OCTAVE_FLAG = (exist('OCTAVE_VERSION', 'builtin') ~= 0); % determine if system is Matlab or GNU Octave

nargoutchk(1,3);

void_intensity = all(isnan(intensity));
void_return_number = all(isnan(returnNumber));
void_return_total = all(isnan(returnTotal));
void_rgb = all(isnan(rgb(:)));


%% filter points

% remove unlabeled points
idxl_filter = false(size(xyz,1),1);
idxl_filter(:,1) = (label ~= 0);
idxl_filter = all(idxl_filter,2);

[label, ~] = grp2idx(label(idxl_filter)); % reassign labels
n_obs = length(unique(label));

%varargout{1} = zeros(length(label),1);
%varargout{1}(idxl_filter) = label;

xyz = xyz(idxl_filter,:);
intensity = double(intensity(idxl_filter));
returnNumber = double(returnNumber(idxl_filter));
returnTotal = double(returnTotal(idxl_filter));
rgb = double(rgb(idxl_filter,:));


%% sort points by label

[label, idxn_sort] = sort(label);
xyz = xyz(idxn_sort,:);
intensity = intensity(idxn_sort);
returnNumber = returnNumber(idxn_sort);
returnTotal = returnTotal(idxn_sort);
rgb = rgb(idxn_sort,:);


%% compute auxiliary variables

idxl_first = (returnNumber == 1);
idxl_last = (returnNumber == returnTotal);
idxl_single = (returnTotal == 1);

% opacity
opacity = returnNumber ./ returnTotal;

% chromaticity
rg_chromaticity = bsxfun(@rdivide, rgb, sum(rgb,2));

% intensity
intensity_first = intensity;
intensity_first(~idxl_first) = nan;

intensity_last = intensity;
intensity_last(~idxl_last) = nan;

intensity_single = intensity;
intensity_single(~idxl_single) = nan;


%% normalize segment coordinates

if arg.Results.verbose
    
    tic
    fprintf('normalizing coordinates...');
    
end



% position proxy (root, centroid, apex)
if isempty(arg.Results.treePos)
    
    proxy = 'apex'; % change to 'root' or 'centroid', if required
    xyz_proxy = zeros(n_obs,3);
    
    switch proxy
        
        case 'root'
            
            z_min = accumarray(label, xyz(:,3), [n_obs,1], @min, nan);
            [~, idxn_min] = ismember([(1:n_obs)', z_min], [label, xyz(:,3)], 'rows');
            xyz_proxy(:,1:2) = xyz(idxn_min,1:2);
            
        case 'centroid'
            
            xyz_proxy(:,1) = accumarray(label, xyz(:,1), [], @mean, nan);
            xyz_proxy(:,2) = accumarray(label, xyz(:,2), [], @mean, nan);
            
        case 'apex'
            
            z_max = accumarray(label, xyz(:,3), [n_obs,1], @max, nan);
            [~, idxn_max] = ismember([(1:n_obs)', z_max], [label, xyz(:,3)], 'rows');
            xyz_proxy(:,1:2) = xyz(idxn_max,1:2);
            
    end
    
    % convert map coordinates (x,y) to image coordinates (column, row)
    RC = round([xyz_proxy(:,1) - refmat(3,1), xyz_proxy(:,2) - refmat(3,2)] / refmat(1:2,:));
    ind = sub2ind(size(dtm), RC(:,1), RC(:,2));
    
    % find terrain elevation at tree proxy coordinates
    xyz_proxy(:,3) = dtm(ind);
    
else
    
    xyz_proxy = arg.Results.treePos;
    
    % convert map coordinates (x,y) to image coordinates (column, row)
    RC = round([xyz_proxy(:,1) - refmat(3,1), xyz_proxy(:,2) - refmat(3,2)] / refmat(1:2,:));
    ind = sub2ind(size(dtm), RC(:,1), RC(:,2));
    
    % find terrain elevation at tree proxy coordinates
    xyz_proxy(:,3) = dtm(ind);
    
end

xyh = xyz(:,1:3) - xyz_proxy(label,:);
h_tot = accumarray(label, xyh(:,3), [], @max, nan);

% normalized height
h_n = xyh(:,3) ./ h_tot(label);

if arg.Results.verbose
    
    fprintf('done!\n');
    toc
    
end


%% normalize intensity

if arg.Results.intensityScaling
    
    if arg.Results.verbose
        
        tic
        fprintf('rescaling intensity...');
        
    end
    
    qlim = quantile(intensity, [0.01, 0.99]); % [0.005, 0.995]
    intensity(intensity < qlim(1)) = qlim(1);
    intensity(intensity > qlim(2)) = qlim(2);
    intensity = (intensity - min(intensity)) ./ (max(intensity) - min(intensity));
    
    if arg.Results.verbose
        
        fprintf('done!\n');
        toc
        
    end
    
end


%% list feature dependencies
warning off
M = struct;

M.Category = {};
M.Name = {};
M.Description = {};
M.Dependencies = {};
M.Scalar = [];
M.Octave = [];
M.ScaleDependant = [];
M.PrintFormat = {};

k = 1;
M.Category{k} = 'Identifier';
M.Name{k} = 'UUID';
M.Description{k} = 'Universally Unique Identifier';
M.Dependencies{k} = [];
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%s';
k = k + 1;

M.Category{k} = 'Identifier';
M.Name{k} = 'LUID';
M.Description{k} = 'Locally Unique Identifier';
M.Dependencies{k} = [];
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%u';
k = k + 1;

M.Category{k} = 'Basic';
M.Name{k} = 'Random';
M.Description{k} = 'Random Number (between 0 and 1)';
M.Dependencies{k} = [];
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.3f';
k = k + 1;

M.Category{k} = 'Basic';
M.Name{k} = 'NPoints';
M.Description{k} = 'Number of points';
M.Dependencies{k} = [];
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = true;
M.PrintFormat{k} = '%u';
k = k + 1;

M.Category{k} = 'Basic';
M.Name{k} = 'XYZ';
M.Description{k} = 'Raw 3D coordinates (x,y,z)';
M.Dependencies{k} = {'NPoints'};
M.Scalar(k) = false;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = [];
k = k + 1;

M.Category{k} = 'Basic';
M.Name{k} = 'XYH';
M.Description{k} = 'Normalized 3D coordinates (x,y,h)';
M.Dependencies{k} = {'NPoints'};
M.Scalar(k) = false;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = [];
k = k + 1;

M.Category{k} = 'Basic';
M.Name{k} = 'UVW';
M.Description{k} = 'Principal components 3D coordinates (u,v,w)';
M.Dependencies{k} = {'XYH'};
M.Scalar(k) = false;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = [];
k = k + 1;

M.Category{k} = 'Basic';
M.Name{k} = 'XPos';
M.Description{k} = 'X coordinate of the position proxy';
M.Dependencies{k} = {'NPoints'};
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.3f';
k = k + 1;

M.Category{k} = 'Basic';
M.Name{k} = 'YPos';
M.Description{k} = 'Y coordinate of the position proxy';
M.Dependencies{k} = {'NPoints'};
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.3f';
k = k + 1;

M.Category{k} = 'Basic';
M.Name{k} = 'ZPos';
M.Description{k} = 'Z coordinate of the position proxy';
M.Dependencies{k} = {'NPoints'};
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.3f';
k = k + 1;

M.Category{k} = 'Basic';
M.Name{k} = 'BBOX2D';
M.Description{k} = '2D bounding box';
M.Dependencies{k} = {'XYZ'};
M.Scalar(k) = false;
M.Octave(k) = true;
M.ScaleDependant(k) = true;
M.PrintFormat{k} = [];
k = k + 1;

M.Category{k} = 'Basic';
M.Name{k} = 'BBOX3D';
M.Description{k} = '3D bounding box';
M.Dependencies{k} = {'XYZ'};
M.Scalar(k) = false;
M.Octave(k) = true;
M.ScaleDependant(k) = true;
M.PrintFormat{k} = [];
k = k + 1;

M.Category{k} = 'Basic';
M.Name{k} = 'CVH2D';
M.Description{k} = '2D convex hull';
M.Dependencies{k} = {'XYZ'};
M.Scalar(k) = false;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = [];
k = k + 1;

M.Category{k} = 'Basic';
M.Name{k} = 'XCVH2D';
M.Description{k} = 'X coordinates of 2D convex hull';
M.Dependencies{k} = {'CVH2D'};
M.Scalar(k) = false;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
k = k + 1;

M.Category{k} = 'Basic';
M.Name{k} = 'YCVH2D';
M.Description{k} = 'Y coordinates of 2D convex hull';
M.Dependencies{k} = {'CVH2D'};
M.Scalar(k) = false;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
k = k + 1;

M.Category{k} = 'Basic';
M.Name{k} = 'CVH3D';
M.Description{k} = '3D convex hull';
M.Dependencies{k} = {'XYZ'};
M.Scalar(k) = false;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = [];
k = k + 1;

M.Category{k} = 'Basic';
M.Name{k} = 'CCH2D';
M.Description{k} = '2D concave hull';
M.Dependencies{k} = {'XYZ'};
M.Scalar(k) = false;
M.Octave(k) = false;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = [];
k = k + 1;

M.Category{k} = 'Basic';
M.Name{k} = 'CCH3D';
M.Description{k} = '3D concave hull';
M.Dependencies{k} = {'XYZ'};
M.Scalar(k) = false;
M.Octave(k) = false;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = [];
k = k + 1;

M.Category{k} = 'PointPatternMetrics';
M.Name{k} = 'HMean';
M.Description{k} = 'Mean of the point heights';
M.Dependencies{k} = [];
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = true;
M.PrintFormat{k} = '%.3f';
k = k + 1;

M.Category{k} = 'PointPatternMetrics';
M.Name{k} = 'HStdDev';
M.Description{k} = 'Standard deviation of the point heights';
M.Dependencies{k} = [];
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.3f';
k = k + 1;

M.Category{k} = 'PointPatternMetrics';
M.Name{k} = 'HCoefVar';
M.Description{k} = 'Coefficient of variation of the point heights';
M.Dependencies{k} = [];
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.3f';
k = k + 1;

M.Category{k} = 'PointPatternMetrics';
M.Name{k} = 'HKurt';
M.Description{k} = 'Kurtosis of the point heights';
M.Dependencies{k} = [];
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.3f';
k = k + 1;

M.Category{k} = 'PointPatternMetrics';
M.Name{k} = 'HSkew';
M.Description{k} = 'Skewness of the point heights';
M.Dependencies{k} = [];
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.3f';
k = k + 1;

M.Category{k} = 'PointPatternMetrics';
M.Name{k} = 'HQ25';
M.Description{k} = '25th percentile of the point heights';
M.Dependencies{k} = [];
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.3f';
k = k + 1;

M.Category{k} = 'PointPatternMetrics';
M.Name{k} = 'HQ50';
M.Description{k} = '50th percentile (median) of the point heights';
M.Dependencies{k} = [];
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.3f';
k = k + 1;

M.Category{k} = 'PointPatternMetrics';
M.Name{k} = 'HQ75';
M.Description{k} = '75th percentile of the point heights';
M.Dependencies{k} = [];
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.3f';
k = k + 1;

M.Category{k} = 'PointPatternMetrics';
M.Name{k} = 'HQ90';
M.Description{k} = '90th percentile of the point heights';
M.Dependencies{k} = [];
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.3f';
k = k + 1;

M.Category{k} = 'PointPatternMetrics';
M.Name{k} = 'PCVar1';
M.Description{k} = 'Variance of the first principal component';
M.Dependencies{k} = {'UVW'};
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = true;
M.PrintFormat{k} = '%.3f';
k = k + 1;

M.Category{k} = 'PointPatternMetrics';
M.Name{k} = 'PCVar2';
M.Description{k} = 'Variance of the second principal component';
M.Dependencies{k} = {'UVW'};
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = true;
M.PrintFormat{k} = '%.3f';
k = k + 1;

M.Category{k} = 'PointPatternMetrics';
M.Name{k} = 'PCVar3';
M.Description{k} = 'Variance of the third principal component';
M.Dependencies{k} = {'UVW'};
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = true;
M.PrintFormat{k} = '%.3f';
k = k + 1;

M.Category{k} = 'PointPatternMetrics';
M.Name{k} = 'CCH3DRa';
M.Description{k} = 'Ratio of points located on the concave hull';
M.Dependencies{k} = {'CCH3D'};
M.Scalar(k) = true;
M.Octave(k) = false;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.2f';
k = k + 1;

M.Category{k} = 'PointPatternMetrics';
M.Name{k} = 'CVH3DRa';
M.Description{k} = 'Fraction of points located on the convex hull';
M.Dependencies{k} = {'CVH3D'};
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.2f';
k = k + 1;

M.Category{k} = 'PointPatternMetrics';
M.Name{k} = 'CCH3DPD';
M.Description{k} = 'Number of points divided by the 3D concave hull volume';
M.Dependencies{k} = {'NPoints', 'CCH3DVol'};
M.Scalar(k) = true;
M.Octave(k) = false;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.2f';
k = k + 1;

M.Category{k} = 'PointPatternMetrics';
M.Name{k} = 'CVH3DPD';
M.Description{k} = 'Number of points divided by the 3D convex hull volume';
M.Dependencies{k} = {'NPoints', 'CVH3DVol'};
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.2f';
k = k + 1;

M.Category{k} = 'ExternalShapeMetrics';
M.Name{k} = 'H';
M.Description{k} = 'Total height';
M.Dependencies{k} = {'XYZ', 'XPos', 'YPos', 'ZPos'};
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = true;
M.PrintFormat{k} = '%.1f';
k = k + 1;

M.Category{k} = 'ExternalShapeMetrics';
M.Name{k} = 'CCH2DArea';
M.Description{k} = 'Area of the 2D concave hull';
M.Dependencies{k} = {'CCH2D'};
M.Scalar(k) = true;
M.Octave(k) = false;
M.ScaleDependant(k) = true;
M.PrintFormat{k} = '%.1f';
k = k + 1;

M.Category{k} = 'ExternalShapeMetrics';
M.Name{k} = 'CCH3DArea';
M.Description{k} = 'Surface area of the 3D concave hull';
M.Dependencies{k} = {'CCH3D'};
M.Scalar(k) = true;
M.Octave(k) = false;
M.ScaleDependant(k) = true;
M.PrintFormat{k} = '%.1f';
k = k + 1;

M.Category{k} = 'ExternalShapeMetrics';
M.Name{k} = 'CCH3DVol';
M.Description{k} = 'Volume of the 3D concave hull';
M.Dependencies{k} = {'CCH3D'};
M.Scalar(k) = true;
M.Octave(k) = false;
M.ScaleDependant(k) = true;
M.PrintFormat{k} = '%.1f';
k = k + 1;

M.Category{k} = 'ExternalShapeMetrics';
M.Name{k} = 'CCH3DSS';
M.Description{k} = 'Specific surface of the 3D concave hull';
M.Dependencies{k} = {'CCH3DVol', 'CCH3DArea'};
M.Scalar(k) = true;
M.Octave(k) = false;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.1f';
k = k + 1;

M.Category{k} = 'ExternalShapeMetrics';
M.Name{k} = 'CVH2DArea';
M.Description{k} = 'Area of the 2D convex hull';
M.Dependencies{k} = {'CVH2D'};
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = true;
M.PrintFormat{k} = '%.1f';
k = k + 1;

M.Category{k} = 'ExternalShapeMetrics';
M.Name{k} = 'CD';
M.Description{k} = 'Diameter of the equivalent area circle';
M.Dependencies{k} = {'CVH2DArea'};
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = true;
M.PrintFormat{k} = '%.1f';
k = k + 1;

M.Category{k} = 'ExternalShapeMetrics';
M.Name{k} = 'CVH3DArea';
M.Description{k} = 'Surface area of the 3D convex hull';
M.Dependencies{k} = {'CVH3D'};
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = true;
M.PrintFormat{k} = '%.1f';
k = k + 1;

M.Category{k} = 'ExternalShapeMetrics';
M.Name{k} = 'CVH3DVol';
M.Description{k} = 'Volume of the 3D convex hull';
M.Dependencies{k} = {'CVH3D'};
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = true;
M.PrintFormat{k} = '%.1f';
k = k + 1;

M.Category{k} = 'ExternalShapeMetrics';
M.Name{k} = 'Cvty3D';
M.Description{k} = '3D convexity (concave surface area divided by convex surface area)';
M.Dependencies{k} = {'CVH3DArea', 'CCH3DArea'};
M.Scalar(k) = true;
M.Octave(k) = false;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.2f';
k = k + 1;

M.Category{k} = 'ExternalShapeMetrics';
M.Name{k} = 'CVH3DLac';
M.Description{k} = 'Convex hull lacunarity (concave volume divided by convex volume)';
M.Dependencies{k} = {'CVH3DVol', 'CCH3DVol'};
M.Scalar(k) = true;
M.Octave(k) = false;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.2f';
k = k + 1;

M.Category{k} = 'ExternalShapeMetrics';
M.Name{k} = 'CVH3DSS';
M.Description{k} = 'Convex hull specific surface ()';
M.Dependencies{k} = {'CVH3DVol', 'CVH3DArea'};
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.2f';
k = k + 1;

M.Category{k} = 'ExternalShapeMetrics';
M.Name{k} = 'AspRatio';
M.Description{k} = 'Aspect ratio (height divided by convex area)';
M.Dependencies{k} = {'CVH2DArea', 'H'};
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.2f';
k = k + 1;

M.Category{k} = 'ExternalShapeMetrics';
M.Name{k} = 'Sphericity';
M.Description{k} = 'Sphericity';
M.Dependencies{k} = {'CVH3DVol', 'CVH3DArea'};
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.2f';
k = k + 1;

M.Category{k} = 'OpacityMetrics';
M.Name{k} = 'Opacity';
M.Description{k} = 'Opacity ()';
M.Dependencies{k} = {'NPoints'};
M.Scalar(k) = false;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.2f';
k = k + 1;

M.Category{k} = 'OpacityMetrics';
M.Name{k} = 'OpacityQ50';
M.Description{k} = '50th percentile (median) of the opacity';
M.Dependencies{k} = [];
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.2f';
k = k + 1;

M.Category{k} = 'OpacityMetrics';
M.Name{k} = 'SRFrac';
M.Description{k} = 'Fraction of single returns';
M.Dependencies{k} = [];
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.2f';
k = k + 1;

M.Category{k} = 'OpacityMetrics';
M.Name{k} = 'FRFrac';
M.Description{k} = 'Fraction of first returns';
M.Dependencies{k} = [];
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.2f';
k = k + 1;

M.Category{k} = 'OpacityMetrics';
M.Name{k} = 'LRFrac';
M.Description{k} = 'Fraction of last returns';
M.Dependencies{k} = [];
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.2f';
k = k + 1;

M.Category{k} = 'IntensityMetrics';
M.Name{k} = 'I';
M.Description{k} = 'Return intensity values';
M.Dependencies{k} = {'NPoints'};
M.Scalar(k) = false;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = [];
k = k + 1;

M.Category{k} = 'IntensityMetrics';
M.Name{k} = 'IQ25';
M.Description{k} = '25th percentile of the return intensity';
M.Dependencies{k} = [];
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.2f';
k = k + 1;

M.Category{k} = 'IntensityMetrics';
M.Name{k} = 'IQ50';
M.Description{k} = '50th percentile of the return intensity';
M.Dependencies{k} = [];
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.2f';
k = k + 1;

M.Category{k} = 'IntensityMetrics';
M.Name{k} = 'IQ75';
M.Description{k} = '75th percentile of the return intensity';
M.Dependencies{k} = [];
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.2f';
k = k + 1;

M.Category{k} = 'IntensityMetrics';
M.Name{k} = 'IQ90';
M.Description{k} = '90th percentile of the return intensity';
M.Dependencies{k} = [];
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.2f';
k = k + 1;

M.Category{k} = 'IntensityMetrics';
M.Name{k} = 'IMean';
M.Description{k} = 'Mean of the return intensity';
M.Dependencies{k} = [];
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.2f';
k = k + 1;

M.Category{k} = 'IntensityMetrics';
M.Name{k} = 'IMax';
M.Description{k} = 'Max of the return intensity';
M.Dependencies{k} = [];
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.2f';
k = k + 1;

M.Category{k} = 'IntensityMetrics';
M.Name{k} = 'IStdDev';
M.Description{k} = 'Standard deviation of the return intensity';
M.Dependencies{k} = [];
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.2f';
k = k + 1;

M.Category{k} = 'IntensityMetrics';
M.Name{k} = 'ICoefVar';
M.Description{k} = 'Coefficient of variation of the return intensity';
M.Dependencies{k} = [];
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.2f';
k = k + 1;

M.Category{k} = 'IntensityMetrics';
M.Name{k} = 'IKurt';
M.Description{k} = 'Kurtosis of the return intensity';
M.Dependencies{k} = [];
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.2f';
k = k + 1;

M.Category{k} = 'IntensityMetrics';
M.Name{k} = 'ISkew';
M.Description{k} = 'Skewness of the return intensity';
M.Dependencies{k} = [];
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.2f';
k = k + 1;

% M.Category{k} = 'IntensityMetrics';
% M.Name{k} = 'IFRQ50';
% M.Description{k} = '50th percentile (median) of the return intensity of first returns';
% M.Dependencies{k} = [];
% M.Scalar(k) = true;
% M.Octave(k) = true;
% M.ScaleDependant(k) = false;
% M.PrintFormat{k} = '%.2f';
% k = k + 1;
% 
% M.Category{k} = 'IntensityMetrics';
% M.Name{k} = 'ILRQ50';
% M.Description{k} = '50th percentile (median) of the return intensity of last returns';
% M.Dependencies{k} = [];
% M.Scalar(k) = true;
% M.Octave(k) = true;
% M.ScaleDependant(k) = false;
% M.PrintFormat{k} = '%.2f';
% k = k + 1;
% 
% M.Category{k} = 'IntensityMetrics';
% M.Name{k} = 'ISRQ50';
% M.Description{k} = '50th percentile of the return intensity of single returns';
% M.Dependencies{k} = [];
% M.Scalar(k) = true;
% M.Octave(k) = true;
% M.ScaleDependant(k) = false;
% M.PrintFormat{k} = '%.2f';
% k = k + 1;

M.Category{k} = 'ColorMetrics';
M.Name{k} = 'RGB';
M.Description{k} = 'Colors (Red, Green, Blue)';
M.Dependencies{k} = {'NPoints'};
M.Scalar(k) = false;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = [];
k = k + 1;

M.Category{k} = 'ColorMetrics';
M.Name{k} = 'Chroma';
M.Description{k} = 'Chromaticity';
M.Dependencies{k} = [];
M.Scalar(k) = false;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = [];
k = k + 1;

M.Category{k} = 'ColorMetrics';
M.Name{k} = 'RChroQ50';
M.Description{k} = 'Median of red chromaticity';
M.Dependencies{k} = {'Chroma'};
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.2f';
k = k + 1;

M.Category{k} = 'ColorMetrics';
M.Name{k} = 'GChroQ50';
M.Description{k} = 'Median of green chromaticity';
M.Dependencies{k} = {'Chroma'};
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.2f';
k = k + 1;

M.Category{k} = 'ColorMetrics';
M.Name{k} = 'RChroHQ50';
M.Description{k} = 'Median of red chromaticity of convex hull points';
M.Dependencies{k} = {'Chroma', 'CVH3D'};
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.2f';
k = k + 1;

M.Category{k} = 'ColorMetrics';
M.Name{k} = 'GChroHQ50';
M.Description{k} = 'Median of green chromaticity of convex hull points';
M.Dependencies{k} = {'Chroma', 'CVH3D'};
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.2f';
k = k + 1;

M.Category{k} = 'ColorMetrics';
M.Name{k} = 'BChroHQ50';
M.Description{k} = 'Median of blue chromaticity of convex hull points';
M.Dependencies{k} = {'Chroma', 'CVH3D'};
M.Scalar(k) = true;
M.Octave(k) = true;
M.ScaleDependant(k) = false;
M.PrintFormat{k} = '%.2f';
k = k + 1;

% Add additional metrics here


%% list available features

idxl_filter = true(1,length(M.Name));

if OCTAVE_FLAG
    
    idxl_filter = idxl_filter & M.Octave;
    
end

if void_intensity || void_return_number || void_return_total
    
    idxl_filter = idxl_filter & ~ismember([M.Category], 'IntensityMetrics');
    
end

if void_return_number || void_return_total
    
    idxl_filter = idxl_filter & ~ismember([M.Category], 'OpacityMetrics');
    
end

if void_rgb
    
    idxl_filter = idxl_filter & ~ismember([M.Category], 'ColorMetrics');
    
end

if ~any(idxl_filter)
    
    error('none of the specified metrics can be computed');
    
end


%% build dependency graph

if arg.Results.verbose
    
    tic
    fprintf('building dependency graph...');
    
end

% create directed adjacency matrix
n = length(M.Name);
A = false(n);
for j = 1:n
    
    if ~isempty(M.Dependencies{j})
        
        A(:, j) = ismember(M.Name, M.Dependencies{j});
        
    end
    
end

% graph structure
G = struct;

% add node names
G.Nodes.Name = M.Name';

% compute indegree
G.Nodes.Degree = sum(A,1)';
G.Nodes.Indegree = sum(A,1)';

% add edges
if any(A(:))
    
    [row, col] = find(A);
    G.Edges.EndNodes = sortrows([row, col], 1);
    
else
    
    G.Edges.EndNodes = [];
    
end

% mark selected nodes
if ismember('all', lower(arg.Results.metrics))
    
    G.Nodes.Selected = true(n,1);
    
else
    
    G.Nodes.Selected = ismember(lower(M.Category)', lower(arg.Results.metrics)) | ismember(lower(M.Name)', lower(arg.Results.metrics));
    
end

%% determine dependency paths

G.Nodes.DependencyPath = cell(n,1);
G.Nodes.Available = idxl_filter';

idxn_start = find(G.Nodes.Selected & G.Nodes.Available)';

for j = idxn_start
    
    idxn_child = j;
    G.Nodes.DependencyPath{j,1} = idxn_child;
    
    while ~isempty(idxn_child)
        
        idxl_adj = ismember(G.Edges.EndNodes(:,2), idxn_child); % find parent node(s) of current node
        idxn_parent = unique(G.Edges.EndNodes(idxl_adj,1));
        G.Nodes.DependencyPath{j,1} = [idxn_parent; G.Nodes.DependencyPath{j,1}];
        idxn_child = idxn_parent;
        
    end

end

G.Nodes.Include = false(n,1);
G.Nodes.Include(cell2mat(G.Nodes.DependencyPath(G.Nodes.Available & G.Nodes.Selected))) = true;

if arg.Results.verbose

    fprintf('done!\n');
    toc

end

%% schedule tasks with topological sorting

Indegree = G.Nodes.Indegree;

L = []; % Empty list that will contain the sorted elements
S = find(G.Nodes.Degree == 0); % Set of all nodes with no incoming edge
e = G.Edges.EndNodes;

while ~isempty(S)
    
    n = S(1);
    S(1) = []; % curent (queue) - remove a node n from S (queue)
    L = [L; n]; % sorted list - add n to tail of L (sorted list)
    
    idxl_adj = e(:,1) == n; % find nodes that have n as a parent
    Indegree(e(idxl_adj,2)) = Indegree(e(idxl_adj,2)) - 1;
    S = [S; intersect(e(idxl_adj,2), find(Indegree == 0))];
    e(idxl_adj,:) = nan; % remove edge e from the graph
    
end

L = L(ismember(L, find(G.Nodes.Include)));

% figure % Matlab only
% plot(digraph(A), 'NodeLabel', G.Nodes.Name)
%
% figure % Matlab only
% plot(digraph(A))


%% compute metrics
metrics = struct();

for j = 1:length(L)
    
    if arg.Results.verbose
        
        tic
        fprintf('computing "%s"...', G.Nodes.Name{L(j)});
        
    end
    
    switch G.Nodes.Name{L(j)}
        
        case 'Random'
            
            % random control number
            metrics.Random = rand(n_obs,1);
            
        case 'UUID'
            
            % universally unique identifier
            metrics.UUID = repmat({nan}, n_obs, 1); % cell(n_obs,1);
            for k = 1:n_obs
                
                metrics.UUID{k} = lower(strcat(dec2hex(randi(16,32,1)-1)'));
                
            end
            
        case 'LUID'
            
            % locally unique identifier (label)
            metrics.LUID = accumarray(label, label, [], @(x) x(1), nan);
            
        case 'NPoints'
            
            % number of points
            metrics.NPoints = accumarray(label, label, [], @numel, {});
            
        case 'XYZ'
            
            % x, y, z coordinates
            metrics.XYZ = mat2cell(xyz, metrics.NPoints, 3);
            
        case 'RGB'
            
            % red, green, blue triplets
            metrics.RGB = mat2cell(rgb, metrics.NPoints, 3);
            
        case 'XPos'
            
            % x root coordinate
            metrics.XPos = xyz_proxy(:,1);
            
        case 'YPos'
            
            % y root coordinate
            metrics.YPos = xyz_proxy(:,2);
            
        case 'ZPos'
            
            % z root coordinate
            metrics.ZPos = xyz_proxy(:,3);
            
        case 'XYH'
            
            % normalized point coordinates
            metrics.XYH = mat2cell(xyh, metrics.NPoints, 3); % xyh
            
        case 'UVW'
            
            % principal components
            metrics.UVW = repmat({nan}, n_obs, 1); % cell(n_obs,1);
            pc_variance = zeros(n_obs, 3);
            
            if OCTAVE_FLAG
                for k = find(metrics.NPoints > 3)'
                    
                    [~, metrics.UVW{k,1}, pc_variance(k,:)] = princomp(metrics.XYH{k,1});
                    
                end
            else
                for k = find(metrics.NPoints > 3)'
                    
                    [~, metrics.UVW{k,1}, pc_variance(k,:)] = pca(metrics.XYH{k,1});
                    
                end
            end
            
        case 'PCVar1'
            
            % spatial coordinates first principal component variance
            metrics.PCVar1 = pc_variance(:,1);
            
        case 'PCVar2'
            
            % spatial coordinates second principal component variance
            metrics.PCVar2 = pc_variance(:,2);
            
        case 'PCVar3'
            
            % spatial coordinates third principal component variance
            metrics.PCVar3 = pc_variance(:,3);
            
        case 'HMean'
            
            % mean of point heights
            metrics.HMean = accumarray(label, h_n, [], @mean, nan);
            
        case 'HMedian'
            
            % median of point heights
            metrics.HMedian = accumarray(label, h_n, [], @median, nan);
            
        case 'HStdDev'
            
            % standard deviation of point heights
            metrics.HStdDev = accumarray(label, h_n, [], @std, nan);
            
        case 'HCoefVar'
            
            % coefficient of variation of point heights
            metrics.HCoefVar = accumarray(label, h_n, [], @(x) std(x)/mean(x), nan);
            
        case 'HKurt'
            
            % kurtosis of point heights
            metrics.HKurt = accumarray(label, h_n, [], @kurtosis, nan);
            
        case 'HSkew'
            
            % skewness of point heights
            metrics.HSkew = accumarray(label, h_n, [], @skewness, nan);
            
        case 'HQ25'
            
            % 25% quantile of point heights
            metrics.HQ25 = accumarray(label, h_n, [], @(x) quantile(x, 0.25), nan);
            
        case 'HQ50'
            
            % 50% quantile of point heights
            metrics.HQ50 = accumarray(label, h_n, [], @(x) quantile(x, 0.5), nan);
            
        case 'HQ75'
            
            % 75% quantile of point heights
            metrics.HQ75 = accumarray(label, h_n, [], @(x) quantile(x, 0.75), nan);
            
        case 'HQ90'
            
            % 90% quantile of point heights
            metrics.HQ90 = accumarray(label, h_n, [], @(x) quantile(x, 0.9), nan);
            
        case 'H'
            
            % total height
            metrics.H = h_tot;
            
        case 'BBOX2D'
            
            metrics.BBOX2D = cell(n_obs,1);
            
            for k = 1:n_obs
                
                metrics.BBOX2D{k,1} = [min(metrics.XYZ{k}(:,1:2), [], 1); max(metrics.XYZ{k}(:,1:2), [], 1)]; % [min(x), min(y); max(x), max(y)]
                
            end
            
        case 'BBOX3D'
            
            metrics.BBOX3D = cell(n_obs,1);
            
            for k = 1:n_obs
                
                metrics.BBOX3D{k,1} = [min(metrics.XYZ{k}(:,1:3), [], 1); max(metrics.XYZ{k}(:,1:3), [], 1)]; % [min(x), min(y), min(z); max(x), max(y), max(z)]
                
            end
            
        case 'CVH2D'
            
            % 2D convex alpha shape
            metrics.CVH2D = repmat({nan}, n_obs, 1); % cell(n_obs,1);
            a_convex_2d = nan(n_obs,1);
            idxl_convex_hull_2d = false(n_obs,1);
            
            for k = find(metrics.NPoints >= 3)' % at least 3 points are required to compute the 2D alpha shape
                
                try
                    
                    [metrics.CVH2D{k,1}, a_convex_2d(k,1)] = convhull(metrics.XYZ{k}(:,1:2));
                    idxl_convex_hull_2d(k,1) = true;
                    
                catch
                   
                    idxl_convex_hull_2d(k,1) = false;
                    
                end
                
            end
            
        case 'XCVH2D'
            
            metrics.XCVH2D = repmat({nan}, n_obs, 1);
            for k = find(idxl_convex_hull_2d)'
                
                metrics.XCVH2D{k,1} = metrics.XYZ{k}(metrics.CVH2D{k},1)';
                
            end
            % metrics.XCVH2D = cellfun(@(x,k) x(k,1)', metrics.XYZ, metrics.CVH2D, 'UniformOutput', false);
            
        case 'YCVH2D'
            
            metrics.YCVH2D = repmat({nan}, n_obs, 1);
            for k = find(idxl_convex_hull_2d)'
                
                metrics.YCVH2D{k,1} = metrics.XYZ{k}(metrics.CVH2D{k},2)';
                
            end
            % metrics.YCVH2D = cellfun(@(x,k) x(k,2)', metrics.XYZ, metrics.CVH2D, 'UniformOutput', false);
            
        case 'CVH3D'
            
            %  3D convex alpha shape
            metrics.CVH3D = repmat({nan}, n_obs, 1);
            v_convex_3d = nan(n_obs,1);
            
            for k = find(metrics.NPoints >= 4)' % at least 4 points are required to compute the 3D alpha shape
                
                try
                    
                    [metrics.CVH3D{k,1}, v_convex_3d(k,1)] = convhulln(metrics.XYZ{k});
                    
                catch
                    
                    metrics.CVH3D{k,1} = nan;
                    
                end
                
            end
            
        case 'CCH2D'
            
            % single region 2D alpha shape
            metrics.CCH2D = repmat({nan}, n_obs, 1);
            
            for k = find(metrics.NPoints >= 3)' % at least 3 points are required to compute the 2D alpha shape
                
                try
                    
                    shp = alphaShape(metrics.XYZ{k}(:,1:2), inf);
                    alpha = criticalAlpha(shp, 'one-region');
                    metrics.CCH2D{k,1} = alphaShape(metrics.XYZ{k}(:,1:2), ...
                        alpha, ...
                        'HoleThreshold', 10^9);
                    
                catch
                    
                    metrics.CCH2D{k,1} = nan;
                    
                end
                
                
            end
            
        case 'CCH3D'
            
            % single region 3D alpha shape
            metrics.CCH3D = repmat({nan}, n_obs, 1); % cell(n_obs,1);
            
            for k = find(metrics.NPoints >= 4)' % at least 4 points are required to compute the 3D alpha shape
                
                try
                    
                    shp = alphaShape(metrics.XYZ{k}, inf);
                    alpha = criticalAlpha(shp, 'one-region');
                    metrics.CCH3D{k,1} = alphaShape(metrics.XYZ{k}, alpha);
                    
                catch
                    
                    metrics.CCH3D{k,1} = nan;
                    
                end
                
            end
            
        case 'CCH2DArea'
            
            % concave area
            metrics.CCH2DArea = nan(n_obs,1);
            idxl_valid = cellfun(@(x) isa(x, 'alphaShape'), metrics.CCH2D);
            for k = find(idxl_valid)'
                
                metrics.CCH2DArea(k) = area(metrics.CCH2D{k});
                
            end
            
        case 'CD'
            
            % crown diameter
            metrics.CD = 2 * sqrt(metrics.CVH2DArea ./ pi);
            
        case 'CCH3DArea'
            
            % concave surface area
            metrics.CCH3DArea = nan(n_obs,1);
            idxl_valid = cellfun(@(x) isa(x, 'alphaShape'), metrics.CCH3D);
            for k = find(idxl_valid)'
                
                metrics.CCH3DArea(k) = surfaceArea(metrics.CCH3D{k});
                
            end
            
        case 'CCH3DVol'
            
            % concave volume
            metrics.CCH3DVol = nan(n_obs,1);
            idxl_valid = cellfun(@(x) isa(x, 'alphaShape'), metrics.CCH3D);
            for k = find(idxl_valid)'
                
                metrics.CCH3DVol(k) = volume(metrics.CCH3D{k});
                
            end
            
        case 'CCH3DSS'
            
            % concave specific surface
            metrics.CCH3DSS = metrics.CCH3DArea ./ metrics.CCH3DVol;
            
        case 'CVH3DArea'
            
            % convex surface area
            metrics.CVH3DArea = nan(n_obs,1);
            idxl_valid = cellfun(@(x) ~isnan(x(1)), metrics.CVH3D);
            
            for k = find(idxl_valid)'
                
                P1 = metrics.XYZ{k,1}(metrics.CVH3D{k,1}(:,1),:);
                P2 = metrics.XYZ{k,1}(metrics.CVH3D{k,1}(:,2),:);
                P3 = metrics.XYZ{k,1}(metrics.CVH3D{k,1}(:,3),:);
                
                % compute triangle normals
                v_n = cross(P1-P2, P1-P3, 2);
                
                % compute triangle areas (= half the length of normal vectors)
                metrics.CVH3DArea(k,1) = sum(0.5 * sqrt(sum(v_n.^2,2)));
                
            end
            
        case 'CCH3DRa'
            
            % concave boundary fraction
            metrics.CCH3DRa = nan(n_obs,1);
            idxl_valid = cellfun(@(x) isa(x, 'alphaShape'), metrics.CCH3D);
            for k = find(idxl_valid)'
                
                [P, ~] = boundaryFacets(metrics.CCH3D{k});
                metrics.CCH3DRa(k,1) = size(P,1);
                
            end
            
            metrics.CCH3DRa = metrics.CCH3DRa ./ metrics.NPoints;
            
        case 'CVH2DArea'
            
            % convex area
            metrics.CVH2DArea = a_convex_2d;
            
        case 'CVH3DVol'
            
            % convex volume
            metrics.CVH3DVol = v_convex_3d;
            
        case 'CVH3DSS'
            
            % convex specific surface
            metrics.CVH3DSS = metrics.CVH3DArea ./ metrics.CVH3DVol;
            
        case 'Cvty3D'
            
            % convexity
            metrics.Cvty3D = metrics.CCH3DArea ./ metrics.CVH3DArea;
            
        case 'Sphericity'
            
            % convex sphericity
            metrics.Sphericity = 3 * sqrt(4*pi) * metrics.CVH3DVol ./ (metrics.CVH3DArea .^(3/2));
            
        case 'AspRatio'
            
            % aspect ratio
            metrics.AspRatio = 2*sqrt(metrics.CVH2DArea/pi) ./ metrics.H;
            
        case 'CVH3DLac'
            
            % convex hull lacunarity
            metrics.CVH3DLac = (metrics.CVH3DVol - metrics.CCH3DVol) ./ metrics.CVH3DVol;
            
        case 'CVH3DRa'
            
            % convex boundary fraction
            metrics.CVH3DRa = cellfun(@(x) size(x,1), metrics.CVH3D) ./ metrics.NPoints;
            
        case 'Opacity'
            
            % opacity
            metrics.Opacity = accumarray(label, opacity, [], @(x) {x}, {nan});
            
        case 'OpacityQ50'
            
            % opacity (50% quantile)
            metrics.OpacityQ50 = accumarray(label, opacity, [], @(x) quantile(x, 0.5), nan);
            
        case 'SRFrac'
            
            % fraction of single returns
            metrics.SRFrac = accumarray(label, idxl_single, [], @(x) nnz(x)/length(x), nan);
            
        case 'FRFrac'
            
            % fraction of first returns
            metrics.FRFrac = accumarray(label, idxl_first, [], @(x) nnz(x)/length(x), nan);
            
        case 'LRFrac'
            
            % fraction of last returns
            metrics.LRFrac = accumarray(label, idxl_last, [], @(x) nnz(x)/length(x), nan);
            
        case 'I'
            
            % intensity
            metrics.I = accumarray(label, intensity, [], @(x) {x}, {nan});
            
        case 'IQ25'
            
            % intensity 25% quantile
            metrics.IQ25 = accumarray(label, intensity, [], @(x) quantile(x, 0.25), nan);
            
        case 'IQ50'
            
            % intensity 50% quantile
            metrics.IQ50 = accumarray(label, intensity, [], @(x) quantile(x, 0.5), nan);
            
        case 'IQ75'
            
            % intensity 75% quantile
            metrics.IQ75 = accumarray(label, intensity, [], @(x) quantile(x, 0.75), nan);
            
        case 'IQ90'
            
            % intensity 90% quantile
            metrics.IQ90 = accumarray(label, intensity, [], @(x) quantile(x, 0.9), nan);
            
        case 'IMean'
            
            % intensity mean
            metrics.IMean = accumarray(label, intensity, [], @mean, nan);
            
        case 'IMax'
            
            % intensity maximum
            metrics.IMax = accumarray(label, intensity, [], @max, nan);
            
        case 'IStdDev'
            
            % intensity standard deviation
            metrics.IStdDev = accumarray(label, intensity, [], @std, nan);
            
        case 'ICoefVar'
            
            % intensity coefficient of variation
            metrics.ICoefVar = accumarray(label, intensity, [], @(x) std(x)/mean(x), nan);
            
        case 'IKurt'
            
            % intensity kurtosis
            metrics.IKurt = accumarray(label, intensity, [], @kurtosis, nan);
            
        case 'ISkew'
            
            % intensity skewness
            metrics.ISkew = accumarray(label, intensity, [], @skewness, nan);
            
        case 'IFRQ50'
            
            % intensity of first returns 50% quantile
            metrics.IFRQ50 = accumarray(label, intensity_first, [], @(x) nanmedian([x; nan]), nan); % @(x) median(x(~isnan(x)))
            
        case 'ILRQ50'
            
            % intensity of last returns 50% quantile
            metrics.ILRQ50 = accumarray(label, intensity_last, [], @(x) nanmedian([x; nan]), nan);
            
        case 'ISRQ50'
            
            % intensity of single returns 50% quantile
            metrics.ISRQ50 = accumarray(label, intensity_single, [], @(x) nanmedian([x; nan]), nan); % @(x) median([x(~isnan(x)), nan])
            
        case 'CCH3DPD'
            
            % point density in concave hull
            metrics.CCH3DPD = metrics.NPoints ./ metrics.CCH3DVol;
            
        case 'CVH3DPD'
            
            % point density in convex hull
            metrics.CVH3DPD = metrics.NPoints ./ metrics.CVH3DVol;
            
        case 'Chroma'
            
            % chromaticity
            metrics.Chroma = mat2cell(rg_chromaticity, metrics.NPoints, 3);
            
        case 'RChroQ50'
            
            % red chromaticity median
            metrics.RChroQ50 = accumarray(label, rg_chromaticity(:,1), [], @median, nan);
            
        case 'GChroQ50'
            
            % green chromaticity median
            metrics.GChroQ50 = accumarray(label, rg_chromaticity(:,2), [], @median, nan);
            
            
    end
    
    if arg.Results.verbose
        
        fprintf('done!\n');
        toc
        
    end
    
end

% remove dependency fields
metrics = rmfield(metrics, setdiff(fieldnames(metrics), M.Name(logical(G.Nodes.Available & G.Nodes.Selected))));

% reorder field names
metrics = orderfields(metrics, M.Name(logical(G.Nodes.Available & G.Nodes.Selected)));

% set print format and scalar flag
varargout{1} = M.PrintFormat(G.Nodes.Available & G.Nodes.Selected); % print_format
varargout{2} = logical(M.Scalar(G.Nodes.Available & G.Nodes.Selected)); % scalar_flag

% convert to non-scalar structure
sfields = fieldnames(metrics);
m = length(sfields);
n = length(metrics.(sfields{1}));

% convert structure to cell array
C = cell(m, n);
for j = 1:m
    
    if isnumeric(metrics.(sfields{j}))
        
        C(j,:) = num2cell(metrics.(sfields{j}));
        
    else
        
        C(j,:) = metrics.(sfields{j});
        
    end
    
end


% delete records containing NaN
fclass = cellfun(@class, C(:,1), 'UniformOutput', false);
idxl_alphashape = ismember(fclass, 'alphaShape');

idxl_nan = any(cellfun(@(x) any(isnan(x(:))), C(~idxl_alphashape,:), 'UniformOutput', true), 1);
idxl_void = any(cellfun(@(x) isempty(x), C, 'UniformOutput', true), 1);

C(:,idxl_nan | idxl_void) = [];

clear metrics
metrics = cell2struct(C, sfields, 1);

