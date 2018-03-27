function [metrics, varargout] = treeMetrics(label, xyz, classification, intensity, returnNumber, returnTotal, rgb, varargin)
% TREEMETRICS - computes segment metrics
%
% [LABEL, ...] = TREEMETRICS(LABEL, XYZ, CLASSIFICATION, INTENSITY, RETURNNUMBER, RETURNTOTAL, RGB, ...) computes various
% segment metrics based on geometry, intensity, opacity (pulse return) and color characteristics
%
% Syntax:  metrics = treeMetrics(label, xyz, classification, intensity, returnNumber, returnTotal, rgb, ...)
%
% Inputs:
%    label - Nx1 integer vector, point label (individual tree label)
%
%    xyz - Nx3 numeric matrix, point cloud coordinates [x y z]
%
%    classification - Nx1 integer vector, classification
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
%    classTerrain (optional, default: 2) - integer vector, terrain class
%    number(s)
%
%    metrics (optional, default: {'all'}) - cell array of strings, list of
%    metrics to compute for each segment. The list can contain a
%    combination of individual features (e.g. 'UUID', 'ConvexVolume',
%    'TotalHeight', 'IntensityMedian'), feature categories (i.e. 'identifier', 'Basic', 
%    'Spatial point pattern', 'External shape', 'Intensity statistics', 'Opacity statistics', 'Color'), 
%    or 'all' (compute all features)
%
%    treePos (optional, default: []) - Mx3 numeric matrix, xyz coordinates of
%    the tree positions. M is the number of segments (trees). If not specified, 
%    the tree positions are automatically computed using the tree top 
%    (highest point in segment) as a proxy for tree position
%
%    intensityScaling (optional, default: true) - boolean value, rescale
%    intensity between 0 and 1 using the [0.001, 0.99] quantiles as limits
%
%    dependencies (optional, default: false) - boolean value, flag
%    indicating if feature dependencies are included in the output
%    (e.g. the 3D convex hull is a dependency of the volume)
%
%    scalarOnly (optional, default: true) - boolean value, flag indicating
%    if only scalar values should be returned
%
%    verbose (optional, default: true) - boolean value, verbosiy switch
%
% Outputs:
%    metrics - structure, metrics for each segment
%
%    label - Nx1 integer vector, updated point label (individual tree labels)
%
% Example:
%
%    metrics = treeMetrics(label_3d, ...
%         [pc.record.x, pc.record.y, pc.record.z], ...
%         pc.record.classification, ...
%         pc.record.intensity, ...
%         pc.record.return_number, ...
%         pc.record.number_of_returns, ...
%         [pc.record.red, pc.record.green, pc.record.blue], ...
%         'metrics', {'Identifiers', 'Intensity statistics', 'Opacity statistics', 'ConvexVolume', 'ConvexArea'}, ...
%         'intensityScaling', true, ...
%         'scalarOnly', true, ...
%         'dependencies', true, ...
%         'verbose', true);
%
% Other m-files required:
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2017b, GNU Octave 4.2.1 (configured for "x86_64-w64-mingw32")
%
% See also:
%
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory (LASIG)
% Website: http://mparkan.github.io/Digital-Forestry-Toolbox/
% Last revision: March 27, 2018
% Acknowledgments: This work was supported by the Swiss Forestry and Wood
% Research Fund, WHFF (OFEV) - project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'label', @(x) (size(x,2) == 1) && isnumeric(x) && any(x ~= 0));
addRequired(arg, 'xyz', @(x) (size(x,2) == 3) && isnumeric(x));
addRequired(arg, 'classification', @(x) (size(x,2) == 1) && isnumeric(x));
addRequired(arg, 'intensity', @(x) (size(x,2) == 1) && isnumeric(x));
addRequired(arg, 'returnNumber', @(x) (size(x,2) == 1) && isnumeric(x));
addRequired(arg, 'returnTotal', @(x) (size(x,2) == 1) && isnumeric(x));
addRequired(arg, 'rgb', @(x) (size(x,2) == 3) && isnumeric(x));
addParameter(arg, 'classTerrain', 2, @(x) isnumeric(x));
addParameter(arg, 'metrics', {'all'}, @(x) iscell(x) && ~isempty(x));
addParameter(arg, 'treePos', [], @(x) (size(x,2) == 3) && isnumeric(x));
addParameter(arg, 'intensityScaling', true, @(x) islogical(x) && (numel(x) == 1));
addParameter(arg, 'dependencies', false, @(x) islogical(x) && (numel(x) == 1));
addParameter(arg, 'scalarOnly', false, @(x) islogical(x) && (numel(x) == 1));
addParameter(arg, 'verbose', true, @(x) islogical(x) && (numel(x) == 1));

parse(arg, label, xyz, classification, intensity, returnNumber, returnTotal, rgb, varargin{:});

OCTAVE_FLAG = (exist('OCTAVE_VERSION', 'builtin') ~= 0); % determine if system is Matlab or GNU Octave

nargoutchk(1,2);


%% extract terrain points

idxl_ter = ismember(classification, arg.Results.classTerrain);
xyz_ter = unique(xyz(idxl_ter,:), 'rows');


%% filter points

% remove unlabeled points
idxl_filter = false(size(xyz,1),1);
idxl_filter(:,1) = (label ~= 0);
idxl_filter = all(idxl_filter,2);

[label, ~] = grp2idx(label(idxl_filter)); % reassign labels
n_obs = length(unique(label));

varargout{1} = zeros(length(label),1);
varargout{1}(idxl_filter) = label;

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

opacity = returnNumber ./ returnTotal;
rg_chromaticity = bsxfun(@rdivide, rgb, sum(rgb,2));

idxl_first = (returnNumber == 1);
idxl_last = (returnNumber == returnTotal);
idxl_single = (returnTotal == 1);

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
    
    proxy = 'apex';
    % margin = 0.5;
    
    xyz_proxy = zeros(n_obs,3);
    
    switch proxy
        
        case 'root'
            
            [~, idxn_sort] = sort(xyz(:,3), 1, 'ascend'); % sort points by elevation
            [~, idxn_first, ~] = unique(label(idxn_sort));
            xyz_proxy = xyz(idxn_sort(idxn_first),:);
            
        case 'centroid'
            
            xyz_proxy(:,1) = accumarray(label, xyz(:,1), [], @mean, nan);
            xyz_proxy(:,2) = accumarray(label, xyz(:,2), [], @mean, nan);
            
        case 'apex'
            
            [~, idxn_sort] = sort(xyz(:,3), 1, 'descend'); % sort points by elevation
            [~, idxn_first, ~] = unique(label(idxn_sort));
            xyz_proxy = xyz(idxn_sort(idxn_first),:);

    end
    
    xyz_proxy(:,3) = griddata(xyz_ter(:,1), ...
        xyz_ter(:,2), ...
        xyz_ter(:,3), ...
        xyz_proxy(:,1), ...
        xyz_proxy(:,2), ...
        'linear');
    idxl_nan = isnan(xyz_proxy(:,3));
    xyz_proxy(idxl_nan,3) = griddata(xyz_ter(:,1), ...
        xyz_ter(:,2), ...
        xyz_ter(:,3), ...
        xyz_proxy(idxl_nan,1), ...
        xyz_proxy(idxl_nan,2), ...
        'nearest');
    clear xyz_ter
    
else
    
    xyz_proxy = arg.Results.treePos;
    
end

xyh = xyz - xyz_proxy(label,:);
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
meta = struct;

meta.Category = {};
meta.Name = {};
meta.Dependencies = {};
meta.Scalar = [];
meta.Octave = [];

k = 1;
meta.Category{k} = 'Identifier';
meta.Name{k} = 'UUID';
meta.Dependencies{k} = [];
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Identifier';
meta.Name{k} = 'LUID';
meta.Dependencies{k} = [];
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Basic';
meta.Name{k} = 'RandomControl';
meta.Dependencies{k} = [];
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Basic';
meta.Name{k} = 'NPoints';
meta.Dependencies{k} = [];
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Basic';
meta.Name{k} = 'XYZ';
meta.Dependencies{k} = {'NPoints'};
meta.Scalar(k) = false;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Basic';
meta.Name{k} = 'RGB';
meta.Dependencies{k} = {'NPoints'};
meta.Scalar(k) = false;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Basic';
meta.Name{k} = 'X';
meta.Dependencies{k} = {'NPoints'};
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Basic';
meta.Name{k} = 'Y';
meta.Dependencies{k} = {'NPoints'};
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Basic';
meta.Name{k} = 'Z';
meta.Dependencies{k} = {'NPoints'};
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Basic';
meta.Name{k} = 'XYH';
meta.Dependencies{k} = {'NPoints'};
meta.Scalar(k) = false;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Basic';
meta.Name{k} = 'UVW';
meta.Dependencies{k} = {'XYH'};
meta.Scalar(k) = false;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Basic';
meta.Name{k} = 'Opacity';
meta.Dependencies{k} = {'NPoints'};
meta.Scalar(k) = false;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Basic';
meta.Name{k} = 'Intensity';
meta.Dependencies{k} = {'NPoints'};
meta.Scalar(k) = false;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Basic';
meta.Name{k} = 'ConvexHull2D';
meta.Dependencies{k} = {'XYH'};
meta.Scalar(k) = false;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Basic';
meta.Name{k} = 'ConvexHull3D';
meta.Dependencies{k} = {'XYH'};
meta.Scalar(k) = false;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Basic';
meta.Name{k} = 'ConcaveHull2D';
meta.Dependencies{k} = {'XYH'};
meta.Scalar(k) = false;
meta.Octave(k) = false;
k = k + 1;

meta.Category{k} = 'Basic';
meta.Name{k} = 'ConcaveHull3D';
meta.Dependencies{k} = {'XYH'};
meta.Scalar(k) = false;
meta.Octave(k) = false;
k = k + 1;

meta.Category{k} = 'Spatial point pattern';
meta.Name{k} = 'TotalHeight';
meta.Dependencies{k} = {'XYZ', 'X', 'Y', 'Z'};
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Spatial point pattern';
meta.Name{k} = 'HeightMean';
meta.Dependencies{k} = [];
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Spatial point pattern';
meta.Name{k} = 'HeightMedian';
meta.Dependencies{k} = [];
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Spatial point pattern';
meta.Name{k} = 'HeightSD';
meta.Dependencies{k} = [];
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Spatial point pattern';
meta.Name{k} = 'HeightCV';
meta.Dependencies{k} = [];
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Spatial point pattern';
meta.Name{k} = 'HeightKurtosis';
meta.Dependencies{k} = [];
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Spatial point pattern';
meta.Name{k} = 'HeightSkewness';
meta.Dependencies{k} = [];
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Spatial point pattern';
meta.Name{k} = 'HeightQ25';
meta.Dependencies{k} = [];
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Spatial point pattern';
meta.Name{k} = 'HeightQ50';
meta.Dependencies{k} = [];
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Spatial point pattern';
meta.Name{k} = 'HeightQ75';
meta.Dependencies{k} = [];
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Spatial point pattern';
meta.Name{k} = 'HeightQ90';
meta.Dependencies{k} = [];
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Spatial point pattern';
meta.Name{k} = 'PrinCompVar1';
meta.Dependencies{k} = {'UVW'};
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Spatial point pattern';
meta.Name{k} = 'PrinCompVar2';
meta.Dependencies{k} = {'UVW'};
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Spatial point pattern';
meta.Name{k} = 'PrinCompVar3';
meta.Dependencies{k} = {'UVW'};
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'External shape';
meta.Name{k} = 'ConcaveArea';
meta.Dependencies{k} = {'ConcaveHull2D'};
meta.Scalar(k) = true;
meta.Octave(k) = false;
k = k + 1;

meta.Category{k} = 'External shape';
meta.Name{k} = 'ConcaveSurfaceArea';
meta.Dependencies{k} = {'ConcaveHull3D'};
meta.Scalar(k) = true;
meta.Octave(k) = false;
k = k + 1;

meta.Category{k} = 'External shape';
meta.Name{k} = 'ConcaveVolume';
meta.Dependencies{k} = {'ConcaveHull3D'};
meta.Scalar(k) = true;
meta.Octave(k) = false;
k = k + 1;

meta.Category{k} = 'External shape';
meta.Name{k} = 'ConcaveSpecificSurface';
meta.Dependencies{k} = {'ConcaveVolume', 'ConcaveSurfaceArea'};
meta.Scalar(k) = true;
meta.Octave(k) = false;
k = k + 1;

meta.Category{k} = 'External shape';
meta.Name{k} = 'ConcaveBoundaryFraction';
meta.Dependencies{k} = {'ConcaveHull3D'};
meta.Scalar(k) = true;
meta.Octave(k) = false;
k = k + 1;

meta.Category{k} = 'External shape';
meta.Name{k} = 'ConvexArea';
meta.Dependencies{k} = {'ConvexHull2D'};
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'External shape';
meta.Name{k} = 'CrownDiameter';
meta.Dependencies{k} = {'ConvexArea'};
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'External shape';
meta.Name{k} = 'ConvexSurfaceArea';
meta.Dependencies{k} = {'ConvexHull3D'};
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'External shape';
meta.Name{k} = 'ConvexVolume';
meta.Dependencies{k} = {'ConvexHull3D'};
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'External shape';
meta.Name{k} = 'Convexity';
meta.Dependencies{k} = {'ConvexSurfaceArea', 'ConcaveSurfaceArea'};
meta.Scalar(k) = true;
meta.Octave(k) = false;
k = k + 1;

meta.Category{k} = 'External shape';
meta.Name{k} = 'ConvexHullLacunarity';
meta.Dependencies{k} = {'ConvexVolume', 'ConcaveVolume'};
meta.Scalar(k) = true;
meta.Octave(k) = false;
k = k + 1;

meta.Category{k} = 'External shape';
meta.Name{k} = 'ConvexSpecificSurface';
meta.Dependencies{k} = {'ConvexVolume', 'ConvexSurfaceArea'};
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'External shape';
meta.Name{k} = 'ConvexBoundaryFraction';
meta.Dependencies{k} = {'ConvexHull3D'};
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'External shape';
meta.Name{k} = 'AspectRatio';
meta.Dependencies{k} = {'ConvexArea', 'TotalHeight'};
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'External shape';
meta.Name{k} = 'ConcavePointDensity';
meta.Dependencies{k} = {'NPoints', 'ConcaveVolume'};
meta.Scalar(k) = true;
meta.Octave(k) = false;
k = k + 1;

meta.Category{k} = 'External shape';
meta.Name{k} = 'ConvexPointDensity';
meta.Dependencies{k} = {'NPoints', 'ConvexVolume'};
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Opacity statistics';
meta.Name{k} = 'OpacityQ50';
meta.Dependencies{k} = [];
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Opacity statistics';
meta.Name{k} = 'SingleReturnFraction';
meta.Dependencies{k} = [];
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Opacity statistics';
meta.Name{k} = 'FirstReturnFraction';
meta.Dependencies{k} = [];
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Opacity statistics';
meta.Name{k} = 'LastReturnFraction';
meta.Dependencies{k} = [];
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Intensity statistics';
meta.Name{k} = 'IntensityQ25';
meta.Dependencies{k} = [];
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Intensity statistics';
meta.Name{k} = 'IntensityQ50';
meta.Dependencies{k} = [];
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Intensity statistics';
meta.Name{k} = 'IntensityQ75';
meta.Dependencies{k} = [];
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Intensity statistics';
meta.Name{k} = 'IntensityQ90';
meta.Dependencies{k} = [];
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Intensity statistics';
meta.Name{k} = 'IntensityMean';
meta.Dependencies{k} = [];
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Intensity statistics';
meta.Name{k} = 'IntensityMax';
meta.Dependencies{k} = [];
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Intensity statistics';
meta.Name{k} = 'IntensitySD';
meta.Dependencies{k} = [];
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Intensity statistics';
meta.Name{k} = 'IntensityCV';
meta.Dependencies{k} = [];
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Intensity statistics';
meta.Name{k} = 'IntensityKurtosis';
meta.Dependencies{k} = [];
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Intensity statistics';
meta.Name{k} = 'IntensitySkewness';
meta.Dependencies{k} = [];
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Intensity statistics';
meta.Name{k} = 'IntensityFirstReturnQ50';
meta.Dependencies{k} = [];
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Intensity statistics';
meta.Name{k} = 'IntensityLastReturnQ50';
meta.Dependencies{k} = [];
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Intensity statistics';
meta.Name{k} = 'IntensitySingleReturnQ50';
meta.Dependencies{k} = [];
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Color';
meta.Name{k} = 'Chromaticity';
meta.Dependencies{k} = [];
meta.Scalar(k) = false;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Color';
meta.Name{k} = 'ChromaticityRedMedian';
meta.Dependencies{k} = {'Chromaticity'};
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

meta.Category{k} = 'Color';
meta.Name{k} = 'ChromaticityGreenMedian';
meta.Dependencies{k} = {'Chromaticity'};
meta.Scalar(k) = true;
meta.Octave(k) = true;
k = k + 1;

% Add additional metrics here


%% build dependency graph

if arg.Results.verbose
    
    tic
    fprintf('building dependency graph...');
    
end

if ismember('all', arg.Results.metrics)
    
    idxl_filter = true(length(meta.Name),1);
    
else
    
    idxl_filter = false(length(meta.Name),1);
    
    for j = 1:length(meta.Name)
        
        if any(ismember([meta.Category(j), meta.Name(j)], arg.Results.metrics))
            
            idxl_filter(j) = true;
            if ~isempty(meta.Dependencies{j})
                
                idxl_filter(ismember(meta.Name, meta.Dependencies{j})) = true;
                
            end
            
        end
        
    end
    
    idxl_filter = false(length(meta.Name),1);
    idxl_filter(ismember(meta.Category, arg.Results.metrics)) = true;
    idxl_filter(ismember(meta.Name, arg.Results.metrics)) = true;
    
end

if OCTAVE_FLAG
    
    idxl_filter(~meta.Octave) = false;
    
end

if ~any(idxl_filter)
    
    error('none of the specified metrics can be computed');
    
end

% create directed adjacency matrix
n = length(meta.Name);
A = false(n);
for j = 1:n
    
    if ~isempty(meta.Dependencies{j})
        
        A(:, j) = ismember(meta.Name, meta.Dependencies{j});
        
    end
    
end

% graph structure
G = struct;

% add node names
G.Nodes.Name = meta.Name;

% compute indegree
G.Nodes.Degree = sum(A,1);

% add edges
if any(A(:))
    
    [row, col] = find(A);
    G.Edges.EndNodes = sortrows([row, col], 1);
    
else
    
    G.Edges.EndNodes = [];
    
end

% topological sorting with depth first search
idxn_start = find(idxl_filter);
order = idxn_start(1);
for j = 1:length(idxn_start)
    
    discovered = false(n,1);
    S = idxn_start(j);
    order = [order; idxn_start(j)];
    
    while ~isempty(S) % while S is not empty
        
        v = S(end);
        S(end) = [];
        
        if ~discovered(v)
            
            discovered(v) = true;
            idxl_adj = G.Edges.EndNodes(:,2) == v;
            S = [S; G.Edges.EndNodes(idxl_adj,1)];
            order = [order; G.Edges.EndNodes(idxl_adj,1)];
            
        end
        
    end
    
end

L = unique(flipud(order), 'first');

if arg.Results.verbose
    
    fprintf('done!\n');
    toc
    
end


%% compute metrics

metrics = struct();

for j = 1:length(L)
    
    if arg.Results.verbose
        
        tic
        fprintf('computing "%s"...', G.Nodes.Name{L(j)});
        
    end
    
    switch G.Nodes.Name{L(j)}
        
        case 'RandomControl'
            
            % random control number
            metrics.RandomControl = rand(n_obs,1);
            
        case 'UUID'
            
            % universally unique identifier
            metrics.UUID = cell(n_obs,1);
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
            
        case 'X'
            
            % x root coordinate
            metrics.X = xyz_proxy(:,1);
            
        case 'Y'
            
            % y root coordinate
            metrics.Y = xyz_proxy(:,2);
            
        case 'Z'
            
            % z root coordinate
            metrics.Z = xyz_proxy(:,3);
            
        case 'XYH'
            
            % normalized point coordinates
            metrics.XYH = mat2cell(xyh, metrics.NPoints, 3); % xyh
            
        case 'UVW'
            
            % principal components
            metrics.UVW = cell(n_obs,1);
            pc_variance = zeros(n_obs, 3);
            
            for k = 1:n_obs
                
                % [~, metrics.UVW{k,1}, ~] = pca(metrics.XYH{k,1}); Matlab only
                [~, metrics.UVW{k,1}, pc_variance(k,:)] = princomp(metrics.XYH{k,1});
                
            end
            
        case 'PrinCompVar1'
            
            % spatial coordinates first principal component variance
            metrics.PrinCompVar1 = pc_variance(:,1);
            
        case 'PrinCompVar2'
            
            % spatial coordinates second principal component variance
            metrics.PrinCompVar2 = pc_variance(:,2);
            
        case 'PrinCompVar3'
            
            % spatial coordinates third principal component variance  
            metrics.PrinCompVar3 = pc_variance(:,3); 
            
        case 'HeightMean'
            
            % mean of point heights
            metrics.HeightMean = accumarray(label, h_n, [], @mean, nan);
            
        case 'HeightMedian'
            
            % median of point heights
            metrics.HeightMedian = accumarray(label, h_n, [], @median, nan);
            
        case 'HeightSD'
            
            % standard deviation of point heights
            metrics.HeightSD = accumarray(label, h_n, [], @std, nan);
            
        case 'HeightCV'
            
            % coefficient of variation of point heights
            metrics.HeightCV = accumarray(label, h_n, [], @(x) std(x)/mean(x), nan);
            
        case 'HeightKurtosis'
            
            % kurtosis of point heights
            metrics.HeightKurtosis = accumarray(label, h_n, [], @kurtosis, nan);
            
        case 'HeightSkewness'
            
            % skewness of point heights
            metrics.HeightSkewness = accumarray(label, h_n, [], @skewness, nan);
            
        case 'HeightQ25'
            
            % 25% quantile of point heights
            metrics.HeightQ25 = accumarray(label, h_n, [], @(x) quantile(x, 0.25), nan);
            
        case 'HeightQ50'
            
            % 50% quantile of point heights
            metrics.HeightQ50 = accumarray(label, h_n, [], @(x) quantile(x, 0.5), nan);
            
        case 'HeightQ75'
            
            % 75% quantile of point heights
            metrics.HeightQ75 = accumarray(label, h_n, [], @(x) quantile(x, 0.75), nan);
            
        case 'HeightQ90'
            
            % 90% quantile of point heights
            metrics.HeightQ90 = accumarray(label, h_n, [], @(x) quantile(x, 0.9), nan);
            
        case 'TotalHeight'
            
            % total height
            metrics.TotalHeight = h_tot;
            
        case 'ConvexHull2D'
            
            % 2D convex alpha shape
            metrics.ConvexHull2D = cell(n_obs,1);
            a_convex_2d = nan(n_obs,1);
            
            for k = 1:n_obs
                
                try
                    
                    [metrics.ConvexHull2D{k,1}, a_convex_2d(k,1)] = convhulln(metrics.XYH{k}(:,1:2));
                    
                catch
                    
                end
                
            end
            
        case 'ConvexHull3D'
            
            %  3D convex alpha shape
            metrics.ConvexHull3D = cell(n_obs,1);
            v_convex_3d = nan(n_obs,1);
            
            for k = 1:n_obs
                
                try
                    
                    [metrics.ConvexHull3D{k,1}, v_convex_3d(k,1)] = convhulln(metrics.XYH{k});
                    
                catch
                    
                    metrics.ConvexHull3D{k,1} = nan;
                    
                end
                
            end
            
        case 'ConcaveHull2D'
            
            % single region 2D alpha shape
            metrics.ConcaveHull2D = cell(n_obs,1);
            
            for k = 1:n_obs
                
                try
                    
                    shp = alphaShape(metrics.XYH{k}(:,1:2), 5);
                    alpha = criticalAlpha(shp, 'one-region');
                    metrics.ConcaveHull2D{k,1} = alphaShape(metrics.XYH{k}(:,1:2), ...
                        alpha, ...
                        'HoleThreshold', 10^9);
                    
                catch
                    
                    metrics.ConcaveHull2D{k,1} = nan;
                    
                end
                
                
            end
            
        case 'ConcaveHull3D'
            
            % single region 3D alpha shape
            metrics.ConcaveHull3D = cell(n_obs,1);
            
            for k = 1:n_obs
                
                try
                    
                    shp = alphaShape(metrics.XYH{k}, 5);
                    alpha = criticalAlpha(shp, 'one-region');
                    metrics.ConcaveHull3D{k,1} = alphaShape(metrics.XYH{k}, alpha);
                    
                catch
                    
                    metrics.ConcaveHull3D{k,1} = nan;
                    
                end
                
            end
            
        case 'ConcaveArea'
            
            % concave area
            metrics.ConcaveArea = nan(n_obs,1);
            idxl_valid = cellfun(@(x) isa(x, 'alphaShape'), metrics.ConcaveHull2D);
            for k = find(idxl_valid)'
                
                metrics.ConcaveArea(k) = area(metrics.ConcaveHull2D{k});
                
            end

        case 'CrownDiameter'
            
            % crown diameter
            metrics.CrownDiameter = round(2 * sqrt(metrics.ConvexArea ./ pi),1);
            
        case 'ConcaveSurfaceArea'
     
            % concave surface area
            metrics.ConcaveSurfaceArea = nan(n_obs,1);
            idxl_valid = cellfun(@(x) isa(x, 'alphaShape'), metrics.ConcaveHull3D);
            for k = find(idxl_valid)'
                
                metrics.ConcaveSurfaceArea(k) = surfaceArea(metrics.ConcaveHull3D{k});
                
            end
            
        case 'ConcaveVolume'
            
            % concave volume
            metrics.ConcaveVolume = nan(n_obs,1);
            idxl_valid = cellfun(@(x) isa(x, 'alphaShape'), metrics.ConcaveHull3D);
            for k = find(idxl_valid)'
                
                metrics.ConcaveVolume(k) = volume(metrics.ConcaveHull3D{k});
                
            end

        case 'ConcaveSpecificSurface'
            
            % concave specific surface
            metrics.ConcaveSpecificSurface = metrics.ConcaveSurfaceArea ./ metrics.ConcaveVolume;
            
        case 'ConvexSurfaceArea'
            
            % convex surface area
            metrics.ConvexSurfaceArea = nan(n_obs,1);
            idxl_valid = cellfun(@(x) ~isnan(x(1)), metrics.ConvexHull3D);
             
            for k = find(idxl_valid)'
               
                P1 = metrics.XYH{k,1}(metrics.ConvexHull3D{k,1}(:,1),:);
                P2 = metrics.XYH{k,1}(metrics.ConvexHull3D{k,1}(:,2),:);
                P3 = metrics.XYH{k,1}(metrics.ConvexHull3D{k,1}(:,3),:);
                
                % compute triangle normals
                v_n = cross(P1-P2, P1-P3, 2);
                
                % compute triangle areas (= half the length of normal vectors)
                metrics.ConvexSurfaceArea(k,1) = sum(0.5 * sqrt(sum(v_n.^2,2)));
                
            end
            
        case 'ConcaveBoundaryFraction'
            
            % concave boundary fraction
            metrics.ConcaveBoundaryFraction = nan(n_obs,1);
            idxl_valid = cellfun(@(x) isa(x, 'alphaShape'), metrics.ConcaveHull3D);
            for k = find(idxl_valid)'
                
                [P, ~] = boundaryFacets(metrics.ConcaveHull3D{k});
                metrics.ConcaveBoundaryFraction(k,1) = size(P,1);
                
            end
            
            metrics.ConcaveBoundaryFraction = metrics.ConcaveBoundaryFraction ./ metrics.NPoints;
            
        case 'ConvexArea'
            
            % convex area
            metrics.ConvexArea = a_convex_2d;

        case 'ConvexVolume'
            
            % convex volume
            metrics.ConvexVolume = v_convex_3d;
            
        case 'ConvexSpecificSurface'
            
            % convex specific surface
            metrics.ConvexSpecificSurface = metrics.ConvexSurfaceArea ./ metrics.ConvexVolume;
            
        case 'Convexity'
            
            % convexity
            metrics.Convexity = metrics.ConcaveSurfaceArea ./ metrics.ConvexSurfaceArea;
            
        case 'Sphericity'
            
            % sphericity
            metrics.Sphericity = 3 * sqrt(4*pi) * metrics.Volume ./ (metrics.SurfaceArea .^(3/2));
            
        case 'AspectRatio'
            
            % aspect ratio
            metrics.AspectRatio = 2*sqrt(metrics.ConvexArea/pi) ./ metrics.TotalHeight;
            
        case 'ConvexHullLacunarity'
            
            % convex hull lacunarity
            metrics.ConvexHullLacunarity = (metrics.ConvexVolume - metrics.ConcaveVolume) ./ metrics.ConvexVolume;
            
        case 'ConvexBoundaryFraction'
            
            % convex boundary fraction
            metrics.ConvexBoundaryFraction = cellfun(@(x) size(x,1), metrics.ConvexHull3D) ./ metrics.NPoints;

        case 'Opacity'
            
            % opacity
            metrics.Opacity = accumarray(label, opacity, [], @(x) {x}, {nan});
            
        case 'OpacityQ50'
            
            % opacity (50% quantile)
            metrics.OpacityQ50 = accumarray(label, opacity, [], @(x) quantile(x, 0.5), nan);
            
        case 'SingleReturnFraction'
            
            % fraction of single returns
            metrics.SingleReturnFraction = accumarray(label, idxl_single, [], @(x) nnz(x)/length(x), nan);
            
        case 'FirstReturnFraction'
            
            % fraction of first returns
            metrics.FirstReturnFraction = accumarray(label, idxl_first, [], @(x) nnz(x)/length(x), nan);
            
        case 'LastReturnFraction'
            
            % fraction of last returns
            metrics.LastReturnFraction = accumarray(label, idxl_last, [], @(x) nnz(x)/length(x), nan);
            
        case 'Intensity'
            
            % intensity
            metrics.Intensity = accumarray(label, intensity, [], @(x) {x}, {nan});
            
        case 'IntensityQ25'
            
            % intensity 25% quantile
            metrics.IntensityQ25 = accumarray(label, intensity, [], @(x) quantile(x, 0.25), nan);
            
        case 'IntensityQ50'
            
            % intensity 50% quantile
            metrics.IntensityQ50 = accumarray(label, intensity, [], @(x) quantile(x, 0.5), nan);
            
        case 'IntensityQ75'
            
            % intensity 75% quantile
            metrics.IntensityQ75 = accumarray(label, intensity, [], @(x) quantile(x, 0.75), nan);
            
        case 'IntensityQ90'
            
            % intensity 90% quantile
            metrics.IntensityQ90 = accumarray(label, intensity, [], @(x) quantile(x, 0.9), nan);
            
        case 'IntensityMean'
            
            % intensity mean
            metrics.IntensityMean = accumarray(label, intensity, [], @mean, nan);
            
        case 'IntensityMax'
            
            % intensity maximum
            metrics.IntensityMax = accumarray(label, intensity, [], @max, nan); 
            
        case 'IntensitySD'
            
            % intensity standard deviation
            metrics.IntensitySD = accumarray(label, intensity, [], @std, nan); 
            
        case 'IntensityCV'
            
            % intensity coefficient of variation
            metrics.IntensityCV = accumarray(label, intensity, [], @(x) std(x)/mean(x), nan);
            
        case 'IntensityKurtosis'
            
            % intensity kurtosis
            metrics.IntensityKurtosis = accumarray(label, intensity, [], @kurtosis, nan);
            
        case 'IntensitySkewness'
            
            % intensity skewness
            metrics.IntensitySkewness = accumarray(label, intensity, [], @skewness, nan);
            
        case 'IntensityFirstReturnQ50'
            
            % intensity of first returns 50% quantile
            metrics.IntensityFirstReturnQ50 = accumarray(label, intensity_first, [], @(x) nanmedian([x; nan]), nan); % @(x) median(x(~isnan(x)))
            
        case 'IntensityLastReturnQ50'
            
            % intensity of last returns 50% quantile
            metrics.IntensityLastReturnQ50 = accumarray(label, intensity_last, [], @(x) nanmedian([x; nan]), nan);
            
        case 'IntensitySingleReturnQ50'
            
            % intensity of single returns 50% quantile
            metrics.IntensitySingleReturnQ50 = accumarray(label, intensity_single, [], @(x) nanmedian([x; nan]), nan); % @(x) median([x(~isnan(x)), nan])
            
        case 'ConcavePointDensity'
            
            % point density in concave hull
            metrics.ConcavePointDensity = metrics.NPoints ./ metrics.ConcaveVolume;
            
        case 'ConvexPointDensity'
            
            % point density in convex hull
            metrics.ConvexPointDensity = metrics.NPoints ./ metrics.ConvexVolume;
            
        case 'Chromaticity'
            
            % chromaticity
            metrics.Chromaticity = mat2cell(rg_chromaticity, metrics.NPoints, 3);
            
        case 'ChromaticityRedMedian'
            
            % red chromaticity median
            metrics.ChromaticityRedMedian = accumarray(label, rg_chromaticity(:,1), [], @median, nan);
            
        case 'ChromaticityGreenMedian'
            
            % green chromaticity median
            metrics.ChromaticityGreenMedian = accumarray(label, rg_chromaticity(:,2), [], @median, nan);
            
    end
    
    if arg.Results.verbose
        
        fprintf('done!\n');
        toc
        
    end
    
end

% remove dependencies from metrics
if  ~arg.Results.dependencies
   
   metrics = rmfield(metrics, intersect(G.Nodes.Name(L), G.Nodes.Name(~idxl_filter)));
   
end

% remove non-scalar metrics
if arg.Results.scalarOnly
    
    metrics = rmfield(metrics, setdiff(fieldnames(metrics), meta.Name(logical(meta.Scalar))));
    
end
