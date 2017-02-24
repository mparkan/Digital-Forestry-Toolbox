function [color, rgb] = clusterColor(X, label, varargin)
% CLUSTERCOLOR - assigns distinct colors to adjacent 3D point clusters (topological graph colouring)
% using the largest degree first heuristic described in [1] and [2].
% [COLOR, RGB] = CLUSTERCOLOR(X, LABEL, ...) assigns a distinct color index (COLOR) and RGB triplet (RGB) to each
% 3D labelled point cluster specified with X and LABEL. Clusters are considered adjacent if their 2D or 3D extents 
% (plus an optional buffer) intersect. The colormaps are based on the values provided in [3] and [4].
%
% [1] Welsh Dominic JA, and B. Powell Martin. "An upper bound for the chromatic number of a graph and its application to timetabling problems." 
% The Computer Journal 10.1 (1967): 85-86.
% [2] Kosowski Adrian and Manuszewski, Krzysztof. "Classical coloring of graphs." Contemporary Mathematics 352 (2004): 1-20.
% [3] Brewer Cynthia and Harrower Marc, "ColorBrewer 2.0 - Color advice for
% cartography", Pennsylvania State University, http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=12
% [4] Jacomy Mathieu, "I want hue - Colors for data scientists", Sciences-Po Medialab, http://tools.medialab.sciences-po.fr/iwanthue/
%
% Syntax:  [color, rgb] = clusterColor(X, label, ...)
%
% Inputs:
%    X - Nx2 or Nx3 numeric matrix, 2D or 3D point coordinates [x y] or [x y z]
%
%    label - Nx1 integer vector, individual cluster label
%
%    adjacency (optional, default: '3d') - string ('2d' or '3d') -
%    adjacency mode. In 2D mode only the horizontal distance between clusters is used.
%
%    buffer (optional, default: 2) - float, width of the buffer around each cluster. 
%    Points are voxelized to the buffer resolution.
%
%    colormap (optional, default: 'cmap12') - string, 'hsv', 'cmap12', 'cmap25'
%
%    unlabelledColor (optional, default: [0,0,0]) - 3x1 numeric matrix,
%    RGB triplet associated with label 0 (i.e. unlabelled)
%
%    verbose (optional, default: true) - boolean value, verbosiy switch
%
%    fig (optional, default: false) - boolean value, switch to plot figures
%
% Outputs:
%    color - Nx1 integer vector, distinct color index
%
%    rgb - Mx3 numeric matrix, RGB triplets associated with the 3D points,
%    unlabelled points are in black by default [0,0,0]
%
% Example:
%
%    [color, rgb] = clusterColor([x y z], ...
%                label, ...
%                'adjacency', '3d', ...   
%                'buffer', 3, ...
%                'colormap', 'cmap12', ...
%                'unlabelledColor', [0.1, 0.1, 0.1], ...
%                'fig', true, ...
%                'verbose', true);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2016b
%
% See also:
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory (LASIG)
% Website: http://mparkan.github.io/Digital-Forestry-Toolbox/
% Last revision: February 24, 2017
% Acknowledgments: This work was supported by the Swiss Forestry and Wood
% Research Fund, WHFF (OFEV) - project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'X', @(x) (size(x,2) >= 2) && (size(x,2) <= 3) && isnumeric(x));
addRequired(arg, 'label', @(x) (size(x,2) == 1) && isnumeric(x));
addParameter(arg, 'adjacency', '3d', @(x) ismember(x, {'2d', '3d'}));
addParameter(arg, 'buffer', 2, @(x) isnumeric(x) && (numel(x) == 1) && x >= 0);
addParameter(arg, 'colormap', 'cmap12', @(x) ismember(x, {'auto', 'hsv', 'cmap12', 'cmap25'}));
addParameter(arg, 'unlabelledColor', [0,0,0], @(x) all(size(x) == [1 3]) && all(x >= 0) && all(x <= 1) && isnumeric(x)); 
addParameter(arg, 'fig', false, @(x) islogical(x) && (numel(x) == 1));
addParameter(arg, 'verbose', true, @(x) islogical(x) && (numel(x) == 1));

parse(arg, X, label, varargin{:});

% check point and label size consistency

if size(X,1) ~= size(label,1)
    
    error('X and label arrays do not have consistent dimensions.');
    
end


%% compute adjacency

% select labelled points only
idxl_labelled = (label ~= 0);
X = X(idxl_labelled,:);
label = label(idxl_labelled);


if ~any(idxl_labelled)
    
    fprintf('Warning: no labelled points.\n');
    color = zeros(size(X,1), 1, 'uint8');
    rgb = arg.Results.unlabelledColor(color+1,:);
    return
    
end

% voxelize


switch arg.Results.adjacency
    
    case '2d'
        
        [X_r, ~, idxn_cell] = unique(round(X(:,1:2) / arg.Results.buffer) * arg.Results.buffer, 'rows');
        
    case '3d' % & size(X,2) == 3
        
        [X_r, ~, idxn_cell] = unique(round(X(:,1:3) / arg.Results.buffer) * arg.Results.buffer, 'rows');
        
end

% find unique labels in each voxel
% [G,ID] = findgroups(A)
Y = splitapply(@(x1){unique(x1)}, label, idxn_cell);

% number of adjacent nodes in each voxel
n_adj = cellfun(@numel, Y);

% define graph edges
idxl_adj = n_adj > 1;
cliques = Y(idxl_adj);

% find all possible node pairs (combinations)
pairs = cell2mat(cellfun(@(x1) nchoosek(x1,2), cliques, 'UniformOutput', false));

if isempty(pairs)
    
    fprintf('Warning: no adjacent clusters, assigning same color to all\n')
    color = uint8(idxl_labelled);
    cmap = [arg.Results.unlabelledColor; 0.8235, 0.4196, 0.4039];
    rgb = cmap(color+1,:);
    return
    
else
    
    % define adjacency matrix
    n_clusters = length(unique(label));
    
    A = false(n_clusters);
    linearInd = sub2ind(size(A), [pairs(:,1); pairs(:,2)], [pairs(:,2); pairs(:,1)]);
    A(linearInd) = true;
    
    % create graph
    G = graph(A);
    G.Nodes.Label = unique(label);
    
end

%% reorder nodes by degree

if arg.Results.verbose
    
    fprintf('reordering nodes by degree...');
    tic
    
end

% G.Nodes.Label = Cluster.Label;
G.Nodes.Degree = degree(G,1:height(G.Nodes))';

[~, idxn_sort] = sort(G.Nodes.Degree);
G = reordernodes(G, idxn_sort(end:-1:1));

if arg.Results.verbose
    
    fprintf('done!\n');
    toc
    
end


%% assign color index to graph nodes

if arg.Results.verbose
    
    fprintf('assigning color index to graph nodes...');
    tic
    
end

n = G.numnodes;
G.Nodes.Color = zeros(n,1, 'uint8');
n_colors = min(G.Nodes.Degree(1), 256);

% traverse nodes
for k = 1:n
    
    idxn_adj = neighbors(G, k);
    idxn_color_adj = unique(G.Nodes.Color(idxn_adj));  % color index in node neighbourhood
    
    if isempty(idxn_color_adj)
        
        idxn_color_max = 0;
        
    else
        
        idxn_color_max = idxn_color_adj(length(idxn_color_adj)); % max color index in node neighbourhood
        
    end
    
    idxn_color_pool = setdiff(1:idxn_color_max, idxn_color_adj);
    
    if isempty(idxn_color_pool)
        
        if idxn_color_max > n_colors % error, color pool is empty
            
            error('\nToo many adjacent clusters');
            
        else
            
            G.Nodes.Color(k) = idxn_color_max + 1; % increment color index
            
        end
        
    else
        
        G.Nodes.Color(k) = min(idxn_color_pool); % use min color index avaiable in pool
        
    end
    
end

m = length(unique(G.Nodes.Color));

if arg.Results.verbose
    
    fprintf('done!\n');
    toc
    
end


%% assign color index to 3D points

if arg.Results.verbose
    
    fprintf('assigning color index to 3D points...');
    tic
    
end

[idxl_sample_lab, locb_sample] = ismember(label, G.Nodes.Label);

color_cluster = zeros(nnz(idxl_labelled), 1, 'uint8');
color_cluster(idxl_sample_lab) = G.Nodes.Color(locb_sample(idxl_sample_lab));

color = zeros(length(idxl_labelled), 1, 'uint8');
color(idxl_labelled) = color_cluster;

if arg.Results.verbose
    
    fprintf('done!\n');
    toc
    
end


%% compute RGB colors

if arg.Results.verbose
    
    fprintf('assigning RGB colors to 3D points...');
    tic
    
end

switch arg.Results.colormap
   
    case 'hsv' % M distinct colors (default Matlab HSV colormap)
        
        cmap = hsv(m);
    
    case 'cmap12' % 12 distinct colors - source: http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=12
        
        cmap = [166,206,227;
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
            177,89,40; 0,0,0] ./ 255;
        
        color_cluster(color_cluster > 12) = mod(color_cluster(color_cluster > 12), 12);
        color(color > 12) = mod(color(color > 12), 12);
        
    case 'cmap25' % 25 distinct colors - source: http://tools.medialab.sciences-po.fr/iwanthue/
        
        cmap = [210,107,103;
            93,192,62;
            202,72,221;
            179,180,58;
            105,71,213;
            87,165,89;
            217,62,166;
            91,169,149;
            216,75,41;
            74,48,146;
            214,143,60;
            152,66,174;
            149,147,78;
            106,123,212;
            137,75,35;
            97,163,203;
            202,52,82;
            65,87,41;
            201,129,214;
            206,149,124;
            72,74,119;
            199,80,133;
            111,47,51;
            192,144,181;
            117,49,107] ./ 255;
        
        color_cluster(color_cluster > 25) = mod(color_cluster(color_cluster > 25), 25);
        color(color > 25) = mod(color(color > 25), 25);
        
end

if m > size(cmap,1)
    
   warning('Number of available colors in chosen colormap (%u) is smaller than required (%u)', size(cmap,1), m)
   
end

rgb_cluster = cmap(color_cluster, :);

rgb = repmat(arg.Results.unlabelledColor, length(idxl_labelled), 1);
rgb(idxl_labelled,:) = cmap(color(idxl_labelled), :);

if arg.Results.verbose
    
    fprintf('done!\n');
    toc
    
end


%% plot figures

if arg.Results.fig
    
    % plot clusters
    switch size(X,2)
        
        case 2
            
            figure
            scatter(X(:,1), ...
                X(:,2), ...
                6, ...
                rgb_cluster, ...
                'Marker', '.')
            axis equal tight vis3d
            title('2D point coloring')
            xlabel('x')
            ylabel('y')
            zlabel('z')
            
        case 3
            
            figure
            scatter3(X(:,1), ...
                X(:,2), ...
                X(:,3), ...
                6, ...
                rgb_cluster, ...
                'Marker', '.')
            axis equal tight vis3d
            title('3D point coloring')
            xlabel('x')
            ylabel('y')
            zlabel('z')
            
    end
    
end
