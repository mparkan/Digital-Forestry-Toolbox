function [labels, tree] = treeRegionGrowing(xyz, classification, varargin)
%TREEREGIONGROWING - Attempts to extract individual tree crowns from a 3D point cloud
% using a modified version of the top down region growing method described in Li et al. (2012) [1]
% [LABELS, TREE] = TREEREGIONGROWING(XYZ, CLASSIFICATION, ...) segments individual tree
% crowns from the 3D point cloud XYZ with classification labels CLASSIFICATION
% and returns the label vector LABELS and a table TREE containing tree data and metrics.
%
% [1] Li, Wenkai, Qinghua Guo, Marek K. Jakubowski, and Maggi Kelly,
% "A New Method for Segmenting Individual Trees from the Lidar Point Cloud",
% Photogrammetric Engineering and Remote Sensing 78, no. 1 (2012): 75–84.
%
% Syntax:  [labels, trees] = treeRegionGrowing(xyz, classification, ...)
%
% Inputs:
%    xyz - Nx3 numeric matrix, 3D point cloud x,y,z coordinates
%
%    classification - Nx1 numeric matrix, 3D point cloud classification
%
%    classTerrain (optional, default: 2) - numeric vector, terrain class
%    number(s)
%
%    classHighVegetation (optional, default: [4, 5, 12]) - numeric vector,
%    high vegetation class number(s)
%
%    coordinateResolution (optional, default: 0.5) - numeric value, spatial resolution used to rasterize 
%    the point cloud. If this parameter is set to zero, the point cloud is not rasterized.
%
%    normalizeElevation (optional, default: true) - boolean value, normalize point
%    elevation by terrain elevation
%
%    peakSearchRadius (optional, default: 1.5) - numeric value, search radius used to
%    define local maxima
%
%    minPeakSpacing (optional, default: [0 1.5; 15 3]) - Kx2 numeric array, minimimum horizontal distance
%    between local maxima as a function of height (e.g. at 15 m height -> 3 m horizontal separation)
%
%    minSampleSize (optional, default: 20) - numeric value, minimum number of points in the sampling
%    buffer to proceed with segmentation
%
%    minSamplingRadius (optional, default: 6) - numeric value, minimum circular sampling
%    radius around each seed point, the sampling radius for each seed is
%    sampling_radius = min(minSamplingRadius + 0.8 * h_seed, maxSamplingRadius);
%
%    maxSamplingRadius (optional, default: 16) - numeric value, maximum circular sampling
%    radius around each seed point, the sampling radius for each seed is
%    sampling_radius = min(minSamplingRadius + 0.8 * h_seed, maxSamplingRadius);
%
%    verbose (optional, default: true) - boolean value, verbosiy switch
%
%    fig (optional, default: true) - boolean value, switch to plot figures
%
% Outputs:
%    labels - Nx1 numeric vector, individual tree crown label for each point
%
%    tree - Table, data and metrics for each tree crown
%
% Example:
%
% pc = LASread('..\data\measurements\vector\als\zh_6995_2710_coniferous.las');
% xyz = [pc.record.x, pc.record.y, pc.record.z];
% classification = pc.record.classification;
%
% [labels, tree] = treeRegionGrowing([x y z], ...
% classification, ...
% 'classTerrain', [2], ...
% 'classHighVegetation', [4,5,12], ...
% 'coordinateResolution', 0.5, ...
% 'normalizeElevation', true, ...
% 'peakSearchRadius', 1.5, ...
% 'minPeakSpacing', [0 1.5; 15 3], ...
% 'minSampleSize', 20, ...
% 'minSamplingRadius', 5, ...
% 'maxSamplingRadius', 18, ...
% 'verbose', true, ...
% 'fig', true);
%
% Other m-files required: pcsel.m
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2016a
%
% See also: TREEGEODESICVOTE
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory
% Website: http://lasig.epfl.ch/
% Last revision: May 23, 2016
% Acknowledgments: This work was supported by the Swiss Forestry and Wood
% Research Fund, WHFF (OFEV) - project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% setup constants

OCTAVE_FLAG = (exist('OCTAVE_VERSION', 'builtin') ~= 0); % determine if system is Matlab or GNU Octave
N_COLORS = 7;


%% check argument validity

arg = inputParser;

addRequired(arg, 'xyz', @(x) (size(x,2) == 3) && isnumeric(x));
addRequired(arg, 'classification', @(x) (size(x,2) == 1) && (size(x,1) == size(xyz,1)) && isnumeric(x));
addParamValue(arg, 'classTerrain', 2, @(x) isnumeric(x));
addParamValue(arg, 'classHighVegetation', [4 5], @(x) isnumeric(x));
addParamValue(arg, 'normalizeElevation', true, @(x) islogical(x) && (numel(x) == 1));
addParamValue(arg, 'coordinateResolution', 0.5, @(x) isnumeric(x) && (numel(x) == 1));
addParamValue(arg, 'peakSearchRadius', 1.5, @(x) isnumeric(x) && (numel(x) == 1));
addParamValue(arg, 'minPeakSpacing', [0 1.5; 15 3], @(x) isnumeric(x) && (size(x,2) == 2));
addParamValue(arg, 'minSampleSize', 20, @(x) isnumeric(x) && (numel(x) == 1));
addParamValue(arg, 'minSamplingRadius', 6, @(x) isnumeric(x) && (numel(x) == 1));
addParamValue(arg, 'maxSamplingRadius', 16, @(x) isnumeric(x) && (numel(x) == 1));
addParamValue(arg, 'fig', true, @(x) islogical(x) && (numel(x) == 1));
addParamValue(arg, 'verbose', true, @(x) islogical(x) && (numel(x) == 1));

parse(arg, xyz, classification, varargin{:});

minPeakSpacing = arg.Results.minPeakSpacing;
minSampleSize = arg.Results.minSampleSize;

x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);


%% compute spatial extents

x_min = min(xyz(:,1));
x_max = max(xyz(:,1));
y_min = min(xyz(:,2));
y_max = max(xyz(:,2));
z_min = min(xyz(:,3));
z_max = max(xyz(:,3));


%% extract relvant point classes

idxl_classTerrain = ismember(classification, arg.Results.classTerrain);
idxl_classHighVegetation = ismember(classification, arg.Results.classHighVegetation);

xyz_veg = xyz(idxl_classHighVegetation,:);
xyz_terr = xyz(idxl_classTerrain,:);


%% adjust coordinate resolution

if arg.Results.coordinateResolution > 0
    
    if arg.Results.verbose
        
        tic
        fprintf('setting coordinate resolution to %2.2f...', arg.Results.coordinateResolution);
        
    end
    
    scale = arg.Results.coordinateResolution;
    
    % rasterize
    xv = x_min-scale:scale:x_max+scale;
    yv = y_min-scale:scale:y_max+scale;
    zv = z_min-scale:scale:z_max+scale;
    
    [sub_crl] = rasterize(xyz_veg, xv, yv, zv);
    nrows = length(yv);
    ncols = length(xv);
    
    % 2D raster
    ind_cr = sub2ind([nrows, ncols], sub_crl(:,2), sub_crl(:,1)); % linear index for 2D raster cells
    [idxn_cells, ~, idxn_cr_unique] = unique(ind_cr);
    [row, col] = ind2sub([nrows, ncols], idxn_cells);
    
    x_voxel = x_min - arg.Results.coordinateResolution/2 + col * arg.Results.coordinateResolution;
    y_voxel = y_min - arg.Results.coordinateResolution/2 + row * arg.Results.coordinateResolution;
    
    z_voxel_sparse = accumarray(ind_cr, xyz_veg(:,3), [nrows*ncols, 1], @(x) max(x), 0, true); % find highest point in each 2D raster cell
    z_voxel = full(z_voxel_sparse(idxn_cells));
    
    x = x_voxel;
    y = y_voxel;
    z = z_voxel;

    if arg.Results.verbose
        
        fprintf('done!\n');
        toc
        
    end
    
end


%% compute terrain model

if arg.Results.verbose
    
    tic
    fprintf('computing terrain model...');
    
end

interpolant = scatteredInterpolant(xyz_terr(:,1), xyz_terr(:,2), xyz_terr(:,3), 'linear', 'nearest');


if arg.Results.verbose
    
    fprintf('done!\n');
    toc
    
end


%% normalize point cloud

if arg.Results.normalizeElevation
    
    if arg.Results.verbose
        
        tic
        fprintf('normalizing point cloud elevation...');
        
    end
    
    h = z - interpolant(x, y);
    
    if arg.Results.verbose
        
        fprintf('done!\n');
        toc
        
    end
    
else
    
    h = z;
    
end

h(h < 0) = 0;


%% segment individual tree crowns

tic

% sort points by height
[h, ix_sort] = sort(h,'descend');
ix_unsort(ix_sort,1) = 1:length(h);

x = x(ix_sort);
y = y(ix_sort);
z = z(ix_sort);
xyzh = [x, y, z, h];
n_points = length(h);
id_point = (1:n_points)';

minPeakSpacing(:,2) = minPeakSpacing(:,2).^2;
idxl_labelled = false(n_points,1);
id_tree = -ones(n_points,1);
id_tree_current = 1;

reverse_str = '';

while ~all(idxl_labelled)
    
    n_remaining = sum(~idxl_labelled);
    
    % define seed point
    idx_seed = find(~idxl_labelled, 1, 'first');
    x_seed = x(idx_seed);
    y_seed = y(idx_seed);
    h_seed = h(idx_seed);
    
    % extract sample around seed
    sampling_radius = min(arg.Results.minSamplingRadius + 0.8 * h_seed, arg.Results.maxSamplingRadius);
    xc = x_seed + sampling_radius * sin(linspace(0,2*pi,10));
    yc = y_seed + sampling_radius * cos(linspace(0,2*pi,10));
    
    % insert dummy point
    xd = xc'; %x_seed + search_radius;
    yd = yc'; %y_seed + search_radius;
    % hd = h_seed;
    
    % extract point cloud sample around seed point
    idxl_sample = inpolygon(x, y, xc, yc) & ~idxl_labelled;
    %sample = [xyh(idxl_sample,:) idx_lm(idxl_sample) id_point(idxl_sample)];
    sample = [xyzh(idxl_sample,:) id_point(idxl_sample)];
    
    sample_size = size(sample,1);
    
    % find all local maxima in sample
    [idx_lm_range,~] = rangesearch(sample(:,1:2), sample(:,1:2), arg.Results.peakSearchRadius);
    sample(:,6) = arrayfun(@(k) all(sample(idx_lm_range{k,1}(2:end),4) < sample(k,4)), 1:sample_size);
    
    if sample_size >= minSampleSize
        
        % initialize P and N
        idxl_P = false(sample_size,1); % set of points belonging to tree
        idxl_labelled(sample(1,5)) = true;
        id_tree(sample(1,5)) = id_tree_current;
        idxl_P(1) = true;
        
        idxl_N = false(sample_size,1); % set of points not belonging to tree
        
        for j = 2:length(sample)
            
            P = sample(idxl_P,1:2);
            N = [[xd yd]; sample(idxl_N,1:2)];
            
            dmin1 = min(sum((P - repmat(sample(j,1:2), size(P,1), 1)).^2, 2)); % minimum distance from u to any point in P_i
            dmin2 = min(sum((N - repmat(sample(j,1:2), size(N,1), 1)).^2, 2)); % minimum distance from u to any point in N_i
            
            if sample(j,6) % if the point is a local maxima

                dt = minPeakSpacing(find(sample(j,4) >= minPeakSpacing(:,1), 1, 'first'), 2);
                
                if dmin1 > dt % if separation distance is larger than peak spacing threshold
                    
                    % point does not belong to tree (N cluster)
                    idxl_N(j) = true;
                    
                elseif (dmin1 <= dt) && (dmin1 <= dmin2) % if point is closer to P cluster and separation distance is smaller than peak spacing threshold
                    
                    % point belongs to tree (P cluster)
                    idxl_labelled(sample(j,5)) = true;
                    id_tree(sample(j,5)) = id_tree_current;
                    idxl_P(j) = true;
                    
                elseif (dmin1 <= dt) && (dmin1 > dmin2) % if point is closer to N cluster and separation distance is smaller than peak spacing threshold
                    
                    % point does not belong to tree (N cluster)
                    idxl_N(j) = true;
                    
                end
                
            else
                
                if dmin1 <= dmin2 % if point is closer to P cluster
                    
                    % point belongs to tree (P cluster)
                    idxl_labelled(sample(j,5)) = true;
                    id_tree(sample(j,5)) = id_tree_current;
                    idxl_P(j) = true;
                    
                else % if point is closer to N cluster
                    
                    % point does not belong to tree (N cluster)
                    idxl_N(j) = true;
                    
                end
                
            end
            
        end
        
        % increment current tree id
        id_tree_current = id_tree_current + 1;
        
    else
        
        idxl_labelled(idxl_sample) = true;
        
    end
    
    msg = sprintf('segmenting %u/%u points...', n_points-n_remaining, n_points);
    fprintf([reverse_str, msg]);
    reverse_str = repmat(sprintf('\b'), 1, length(msg));
    
end

fprintf('done!\n');
toc


%% create label vector

labels = zeros(size(xyz,1), 1, 'uint32');

if arg.Results.coordinateResolution > 0
    
    label_voxels = id_tree(ix_unsort);
    labels(idxl_classHighVegetation) = uint32(label_voxels(idxn_cr_unique));
    
else

    labels = uint32(id_tree(ix_unsort));
    
end


%% build tree objects

if arg.Results.verbose
    
    tic
    fprintf('computing tree metrics...');
    
end

n_trees = max(labels);

data = table();
metrics = table();

warning off

for j = 1:n_trees
    
    idxl_tree = (labels == j);
    
    % point data
    data.nodes(j,1) = {xyz(idxl_tree, :)};
    [z_max, idxn_top] = max(data.nodes{j}(:,3));
    z_min = interpolant(data.nodes{j}(idxn_top, 1), data.nodes{j}(idxn_top, 2));
    data.root(j,1) = {[data.nodes{j}(idxn_top, 1:2), z_min]};
        
    % point metrics
    metrics.centroid(j,1) = {mean(data.nodes{j}(:,1:3),1)};
    metrics.height_vertical(j,1) = z_max - z_min;

end


%% metrics = [tree.metrics];

% assign distinctive colors to adjacent trees
[color_number, ~] = pcsel(reshape(cell2mat(metrics.centroid), 3, [])', N_COLORS); % generate distinctive colors in each neighbourhood
color_number = color_number';
colors = hsv(N_COLORS);

% merge data, metrics and color tables into tree table
tree = table(data, metrics, colors(color_number,:), 'VariableNames',{'data' 'metrics', 'color'});

if arg.Results.verbose
    
    fprintf('done!\n');
    toc
    
end


%% plot segments

if arg.Results.fig
    
    xyz_trees = cell2mat(tree.data.nodes);
    rgb_colors = cell2mat(cellfun(@(color, nodes) repmat(color, size(nodes,1), 1), mat2cell(tree(:,:).color, ones(height(tree),1), 3), tree.data.nodes, 'UniformOutput', false));
    
    figure
    scatter3(xyz_trees(:,1), xyz_trees(:,2), xyz_trees(:,3), 12, ...
        rgb_colors, ...
        'Marker', '.');
    
    axis equal tight vis3d
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('Individual tree crowns')
    
end

end