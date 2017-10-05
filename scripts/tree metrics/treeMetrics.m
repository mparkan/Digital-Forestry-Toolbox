function [metrics, varargout] = treeMetrics(label, xyz, classification, intensity, returnNumber, returnTotal, varargin)
% TREEMETRICS - computes segment features
%
% [LABEL, ...] = TREEMETRICS(LABEL, XYZ, CLASSIFICATION, INTENSITY, RETURNNUMBER, RETURNTOTAL ...) computes various 
% segment metrics based on geometry, intensity and pulse return characteristics
%
% Syntax:  metrics = treeMetrics(label, xyz, classification, intensity, returnNumber, ...)
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
%    classTerrain (optional, default: 2) - integer vector, ground class
%    number(s)
%
%    metrics (optional, default: {'all'}) - cell array of strings, list of
%    metrics to compute for each segment
%    
%    heightCorrection (optional, default: true) - boolean value, height
%    correction flag (adjusts tree height in slopes)
%
%    intensityScaling (optional, default: true) - boolean value, rescale
%    intensity between 0 and 1 using the [0.005, 0.995] quantiles as limits. 
%
%    verbose (optional, default: true) - boolean value, verbosiy switch
%
%    fig (optional, default: true) - boolean value, switch to plot figures
%
% Outputs:
%    metrics - table, metrics for each segment
%
%    label - Nx1 integer vector, updated point label (individual tree labels)
%
% Example:
%
%    metrics = treeMetrics(label, ...
%        xyz, ...
%        classification, ...
%        intensity, ...
%        zeros(size(classification)), ...
%        zeros(size(classification)), ...
%        'metrics', {'all'}, ... 
%        'heightCorrection', false, ...
%        'intensityScaling', true, ...
%        'fig', false, ...
%        'verbose', true);
%
%
% Other m-files required:
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2016b
%
% See also:
%
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory (LASIG)
% Website: http://mparkan.github.io/Digital-Forestry-Toolbox/
% Last revision: October 10, 2017
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
addParameter(arg, 'classTerrain', 2, @(x) isnumeric(x));
addParameter(arg, 'metrics', {'all'}, @(x) iscell(x) && ~isempty(x));
addParameter(arg, 'heightCorrection', true, @(x) islogical(x) && (numel(x) == 1));
addParameter(arg, 'intensityScaling', true, @(x) islogical(x) && (numel(x) == 1));
addParameter(arg, 'fig', false, @(x) islogical(x) && (numel(x) == 1));
addParameter(arg, 'verbose', true, @(x) islogical(x) && (numel(x) == 1));

parse(arg, label, xyz, classification, intensity, returnNumber, returnTotal, varargin{:});


%% list feature dependencies

warning off
meta = table;

k = 1;

meta.Name{k,1} = 'UUID';
meta.Dependencies{k,1} = [];
k = k + 1;

meta.Name{k,1} = 'Label';
meta.Dependencies{k,1} = [];
k = k + 1;

meta.Name{k,1} = 'IndexTree';
meta.Dependencies{k,1} = [];
k = k + 1;

meta.Name{k,1} = 'NPoints';
meta.Dependencies{k,1} = [];
k = k + 1;

meta.Name{k,1} = 'XYZ';
meta.Dependencies{k,1} = [];
k = k + 1;

meta.Name{k,1} = 'Pos';
meta.Dependencies{k,1} = {'XYZ'};
k = k + 1;

meta.Name{k,1} = 'XYZ0';
meta.Dependencies{k,1} = {'XYZ', 'Pos'};
k = k + 1;

meta.Name{k,1} = 'TotalHeight';
meta.Dependencies{k,1} = {'XYZ', 'Pos'};
k = k + 1;

meta.Name{k,1} = 'ConvexHull3D';
meta.Dependencies{k,1} = {'XYZ0'};
k = k + 1;

meta.Name{k,1} = 'ConcaveHull2D';
meta.Dependencies{k,1} = {'XYZ0', 'TotalHeight'};
k = k + 1;

meta.Name{k,1} = 'ConcaveHull3D';
meta.Dependencies{k,1} = {'XYZ0', 'TotalHeight'};
k = k + 1;

meta.Name{k,1} = 'Area';
meta.Dependencies{k,1} = {'ConcaveHull2D'};
k = k + 1;

meta.Name{k,1} = 'CrownDiameter';
meta.Dependencies{k,1} = {'Area'};
k = k + 1;

meta.Name{k,1} = 'SurfaceArea';
meta.Dependencies{k,1} = {'ConcaveHull3D'};
k = k + 1;

meta.Name{k,1} = 'Volume';
meta.Dependencies{k,1} = {'ConcaveHull3D'};
k = k + 1;

meta.Name{k,1} = 'SpecificSurface';
meta.Dependencies{k,1} = {'Volume', 'SurfaceArea'};
k = k + 1;

meta.Name{k,1} = 'ConvexSurfaceArea';
meta.Dependencies{k,1} = {'ConvexHull3D'};
k = k + 1;

meta.Name{k,1} = 'ConvexVolume';
meta.Dependencies{k,1} = {'ConvexHull3D'};
k = k + 1;

meta.Name{k,1} = 'Convexity';
meta.Dependencies{k,1} = {'ConvexSurfaceArea', 'SurfaceArea'};
k = k + 1;

meta.Name{k,1} = 'Sphericity';
meta.Dependencies{k,1} = {'SurfaceArea', 'Volume'};
k = k + 1;

meta.Name{k,1} = 'Opacity';
meta.Dependencies{k,1} = {'Opacity'};
k = k + 1;

meta.Name{k,1} = 'OpacityQ50';
meta.Dependencies{k,1} = {'Opacity'};
k = k + 1;

meta.Name{k,1} = 'Intensity';
meta.Dependencies{k,1} = {'Intensity'};
k = k + 1;

meta.Name{k,1} = 'IntensityQ25';
meta.Dependencies{k,1} = {'Intensity'};
k = k + 1;

meta.Name{k,1} = 'IntensityQ50';
meta.Dependencies{k,1} = {'Intensity'};
k = k + 1;

meta.Name{k,1} = 'IntensityQ75';
meta.Dependencies{k,1} = {'Intensity'};
k = k + 1;

% Add you own custom metrics here


%% build dependency graph

if arg.Results.verbose
    
    tic
    fprintf('building dependency graph...');
    
end

G = digraph;

% add nodes
G = addnode(G, meta.Name);

% add edges
for k = 1:G.numnodes
    
    if ~isempty(meta.Dependencies{k,1})
        
        for j = 1:length(meta.Dependencies{k,1})
            
            G = addedge(G, ...
                meta.Name(k,1), ...
                meta.Dependencies{k,1}(j));
            
        end
        
    end
    
end

H = transreduction(G);

if arg.Results.verbose
    
    fprintf('done!\n');
    toc
    
end


%% determine evaluation order

if ismember('all', arg.Results.metrics)
    
    metrics_list = meta.Name;
    
else
    
    metrics_list = arg.Results.metrics;
    
end

dependencies = cell(length(metrics_list),1);

for j = 1:length(metrics_list)
    
    dependencies{j,1} = bfsearch(H, metrics_list{j});
    
end

id_nodes = unique(vertcat(dependencies{:}));

HS = subgraph(H, id_nodes);

idxn_exec = flip(toposort(HS));

if arg.Results.fig
    
    figure
    plot(HS)
    title('metric dependencies')
    
end


%% reassign labels

idxl_null = (label == 0);
newlabel = zeros(size(label));
[newlabel(~idxl_null), ~] = grp2idx(label(~idxl_null));
label = newlabel;


%% compute terrain model

idxl_ter = ismember(classification, arg.Results.classTerrain);

F = scatteredInterpolant(xyz(idxl_ter,1), ...
    xyz(idxl_ter,2), ...
    xyz(idxl_ter,3), ...
    'linear', ...
    'nearest');


%% filter out abnormal segments

if arg.Results.heightCorrection
    
    % terrain altitude
    z_ter = F(xyz(:,1), xyz(:,2));
    
    htol = 1.5;
    idxn_error = splitapply(@(j,z,zt) {j(z < (max(zt)-htol))}, ...
        (1:length(label))', ...
        xyz(:,3), ...
        z_ter, ...
        label+1);
    
    idxn_error = cell2mat(idxn_error(2:end));
    label(idxn_error) = 0;
    
    varargout{1} = label;
    
end


%% normalize intensity

if ~isa(intensity, 'double')
    
    intensity = double(intensity);
    fprintf('Waring: intensity was cast to double storage format\n');
    
end

if arg.Results.intensityScaling
    
    qlim = quantile(intensity, [0.005, 0.995]);
    intensity(intensity < qlim(1)) = qlim(1);
    intensity(intensity > qlim(2)) = qlim(2);
    intensity = (intensity - min(intensity)) ./ (max(intensity) - min(intensity));
    
end


%% compute opacity

opacity = returnNumber ./ returnTotal;


%% compute metrics

metrics = table;

for j = 1:length(idxn_exec)
    
    if arg.Results.verbose
        
        tic
        fprintf('computing "%s"...', HS.Nodes.Name{idxn_exec(j)});
        
    end
    
    switch HS.Nodes.Name{idxn_exec(j)}
        
        case 'UUID'
            
            metrics.UUID = cell(height(metrics),1);
            
            nseg = length(unique(label(label ~= 0)));
            for k = 1:nseg
                
                metrics.UUID{k} = lower(strrep(char(java.util.UUID.randomUUID), '-', ''));
                
            end
            
        case 'Label'
            
            LABEL = splitapply(@(x) {x}, label, label+1);
            LABEL = LABEL(2:end);
            
            metrics.Label = cellfun(@(x) x(1), ...
                LABEL, ...
                'UniformOutput', true);
            
        case 'IndexTree'
            
            IndexTree = splitapply(@(x) {x}, (1:length(label))', label+1);
            metrics.IndexTree = IndexTree(2:end);
            
        case 'NPoints'
            
            NPoints = splitapply(@(x) numel(x), label, label+1);
            metrics.NPoints = NPoints(2:end);
            
        case 'XYZ'
            
            XYZ = splitapply(@(x) {x}, xyz, label+1);
            metrics.XYZ = XYZ(2:end);
            
        case 'Pos'
            
            [z_max, idxn_apex] = cellfun(@(x) max(x(:,3)), ...
                metrics.XYZ, ...
                'UniformOutput', false);
            
            xyz_apex = cellfun(@(x,k) x(k,:), ...
                metrics.XYZ, ...
                idxn_apex, ...
                'UniformOutput', false);
            
            xyz_root = cellfun(@(x) [x(1:2), F(x(1), x(2))], ...
                xyz_apex, ...
                'UniformOutput', false);
            
            if arg.Results.heightCorrection
                
                z_terr_max = splitapply(@(x) max(x), ...
                    z_ter, ...
                    label+1);
                z_terr_max = z_terr_max(2:end);

                xyz_root_mat = cell2mat(xyz_root);
                z_diff = xyz_root_mat(:,3) - z_terr_max;
                
                z_diff_outlier = median(abs(z_diff)) + 1.5 * iqr(abs(z_diff));
                idxl_outlier = abs(z_diff) > z_diff_outlier;
                
                xyz_root_mat(idxl_outlier,3) = z_terr_max(idxl_outlier);
                xyz_root = mat2cell(xyz_root_mat, ones(1,size(xyz_root_mat,1)), 3);
                
                % height correction flag
                metrics.HeightCorr = idxl_outlier;
                
            end
            
            % root coordinates
            xyz_root_mat = cell2mat(xyz_root);
            metrics.X = xyz_root_mat(:,1);
            metrics.Y = xyz_root_mat(:,2);
            metrics.Z = xyz_root_mat(:,3);
            
        case 'XYZ0'
            
            % normalized point coordinates
            metrics.XYZ0 = cellfun(@(x1,x2) bsxfun(@minus, x1, x2), ...
                metrics.XYZ, ...
                xyz_root, ...
                'UniformOutput', false);
            
        case 'TotalHeight'
            
            % total height (oblique height)
            metrics.TotalHeight = cellfun(@(x1,x2) norm(x1 - x2), ...
                xyz_root, ...
                xyz_apex, ...
                'UniformOutput', true);

        case 'ConvexHull3D'
            
            % compute 3D convex alpha shape
            metrics.ConvexHull3D = cell(height(metrics),1);
            
            for k = 1:height(metrics)
                
                metrics.ConvexHull3D{k,1} = alphaShape([metrics.XYZ0{k}; 0,0,0], ...
                    inf);
                
            end
            
        case 'ConcaveHull2D'
            
            % compute complete 2D alpha shape
            metrics.ConcaveHull2D = cell(height(metrics),1);
            
            for k = 1:height(metrics)
                
                metrics.ConcaveHull2D{k,1} = alphaShape(metrics.XYZ0{k}(:,1:2), ...
                    1.5);
                
            end
            
        case 'ConcaveHull3D'
            
            % compute 3D concave alpha shape
            metrics.ConcaveHull3D = cell(height(metrics),1);
            
            for k = 1:height(metrics)
                
                metrics.ConcaveHull3D{k,1} = alphaShape(metrics.XYZ0{k}, ...
                    min(max([2, 0.1*metrics.TotalHeight(k,1)]), 3));
                
            end
            
        case 'Area'
            
            % total 2D projected area
            metrics.Area = zeros(height(metrics),1);
            
            for k = 1:height(metrics)
                
                if size(metrics.ConcaveHull2D{k,1}.Points,1) > 6
                    
                    metrics.Area(k,1) = round(cellfun(@(x) area(x,1), metrics.ConcaveHull2D(k,1)), 1);
                    
                end
                
            end
            
        case 'CrownDiameter'
            
            % total 2D projected area
            metrics.CrownDiameter = round(2 * sqrt(metrics.Area ./ pi),1);
            
        case 'SurfaceArea'
            
            metrics.SurfaceArea = zeros(height(metrics),1);
            
            for k = 1:height(metrics)
                
                if size(metrics.ConcaveHull3D{k,1}.Points,1) > 6
                    
                    metrics.SurfaceArea(k,1) = round(cellfun(@(x) surfaceArea(x,1), metrics.ConcaveHull3D(k,1)), 1);
                    
                end
                
            end
            
        case 'Volume'
            
            metrics.Volume = zeros(height(metrics),1);
            
            for k = 1:height(metrics)
                
                if size(metrics.ConcaveHull3D{k,1}.Points,1) > 6
                    
                    metrics.Volume(k,1) = round(cellfun(@(x) volume(x), metrics.ConcaveHull3D(k,1)), 1);
                    
                end
                
            end

        case 'SpecificSurface'
            
            metrics.SpecificSurface = metrics.SurfaceArea ./ metrics.Volume;

        case 'ConvexSurfaceArea'
            
            metrics.ConvexSurfaceArea = zeros(height(metrics),1);
            
            for k = 1:height(metrics)
                
                if size(metrics.ConvexHull3D{k,1}.Points,1) > 6
                    
                    metrics.ConvexSurfaceArea(k,1) = round(cellfun(@(x) surfaceArea(x,1), metrics.ConvexHull3D(k,1)), 1);
                    
                end
                
            end

        case 'ConvexVolume'
            
            metrics.ConvexVolume = zeros(height(metrics),1);
            
            for k = 1:height(metrics)
                
                if size(metrics.ConvexHull3D{k,1}.Points,1) > 6
                    
                    metrics.ConvexVolume(k,1) = round(cellfun(@(x) volume(x,1), metrics.ConvexHull3D(k,1)), 1);
                    
                end
                
            end
            
        case 'ConvexSpecificSurface'
            
            metrics.ConvexSpecificSurface = metrics.ConvexSurfaceArea ./ metrics.ConvexVolume;
            
        case 'Convexity'
            
            metrics.Convexity = metrics.ConvexSurfaceArea ./ metrics.SurfaceArea;
            
        case 'Sphericity'
            
            metrics.Sphericity = 3 * sqrt(4*pi) * metrics.Volume ./ (metrics.SurfaceArea .^(3/2));
            
        case 'Opacity'
            
            Opacity = splitapply(@(x) {x}, opacity, label+1);
            metrics.Opacity = Opacity(2:end);
            
        case 'OpacityQ50'
            
            metrics.OpacityQ50 = cellfun(@(x) quantile(x,0.5), metrics.Opacity);
            
        case 'Intensity'
            
            Intensity = splitapply(@(x) {x}, intensity, label+1);
            metrics.Intensity = Intensity(2:end);
            
        case 'IntensityQ25'
            
            metrics.IntensityQ25 = cellfun(@(x) quantile(x,0.25), metrics.Intensity);
            
        case 'IntensityQ50'
            
            metrics.IntensityQ50 = cellfun(@(x) quantile(x,0.5), metrics.Intensity);
            
        case 'IntensityQ75'
            
            metrics.IntensityQ75 = cellfun(@(x) quantile(x,0.75), metrics.Intensity);
            
        % Add you own custom metrics here
            
    end
    
    if arg.Results.verbose
        
        fprintf('done!\n');
        toc
        
    end
    
end