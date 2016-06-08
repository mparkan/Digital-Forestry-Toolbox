function cover = canopyCover(chm, refmat, varargin)
% CANOPYCOVER - computes a canopy crown cover map
%
% COVER = CANOPYCOVER(CHM, REFMAT, ...) computes a canopy crown cover map using 2D convolution or a variant of the Delaunay 
% triangulation based method described in [1]. The cover map is computed from the georeferenced raster canopy height 
% model CHM and tree locations STEMS (for the triangulation method).
%
% [1] L. Eysn, M. Hollaus, K. Schadauer, et N. Pfeifer, « Forest Delineation Based on Airborne LIDAR Data », 
% Remote Sensing, vol. 4, p. 762-783, march 2012.
%
% Syntax:  cover = canopyCover(chm, refmat, ...)
%
% Inputs:
%    chm - RxC numeric matrix, raster Canopy Height Model (CHM) with R rows
%    and C columns
%
%    refmat - 3x2 numeric matrix, spatial referencing matrix such that [map_x map_y] = [row col 1] * refmat
%
%    method (optional, default: 'convolution') - A string indicating which method to use for the
%    crown cover computation, either 'convolution' (default) or 'triangulation'
%
%    stems - Mx3 numeric matrix, column, row and height values for
%    each stem, [col, row, height]
%
%    crowns - RxC numeric matrix, individual tree crown label matrix with R rows
%    and C columns
%
%    minHeight (optional, default: 0.6) - numeric value, minimum canopy
%    height, all values below this threshold are set to zero.
%
%    coverThreshold (optional, default: 0.9) - numeric value between 0 and 1, threshold
%    of Delaunay triangle filling ratio above which the cover is set to 1  
%
%    filter (optional, default: ones(20,20)) - numeric matrix, convolution window (for convolution method only)
%
%    verbose (optional, default: true) - boolean value, verbosiy switch
%
%    fig (optional, default: true) - boolean value, switch to plot figures
%
% Outputs:
%    cover - RxC numeric matrix, canopy cover ratio
%
% Example:
%    info = geotiffinfo('..\data\measurements\raster\chm\so_2014_woodland_pasture.tif');
%    [chm, ~] = geotiffread('..\data\measurements\raster\chm\so_2014_woodland_pasture.tif');
%    
%    [crh, xyh] = treePeaks(chm, info.RefMatrix, ...
%        'method', 'fixedRadius', ...
%        'windowRadius', 3, ...
%        'fig', true);
%
%    labels = treeWatershed(chm, crh, ...
%        'minHeight', 1, ...
%        'fig', true);
%
%    cover = canopyCover(chm, info.RefMatrix, ...
%        'method', 'triangulation', ...
%        'stems', crh, ...
%        'crowns', labels, ...
%        'verbose', true, ...
%        'fig', true);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2016a
%
% See also: treePeaks.m, treeWatershed.m
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory
% Website: http://lasig.epfl.ch/
% Last revision: June 8, 2016
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund (WHFF, OFEV), project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'chm', @isnumeric);
addRequired(arg, 'refmat', @isnumeric);
addParamValue(arg, 'method', 'convolution', @(x) ismember(x, {'convolution', 'triangulation'}));
addParamValue(arg, 'stems', [], @(x) isnumeric(x) && size(x,1) >= 3);
addParamValue(arg, 'crowns', [], @(x) isnumeric(x) && isequal(size(x), size(chm)));
addParamValue(arg, 'minHeight', 1, @(x) isnumeric(x) && (numel(x) == 1));
addParamValue(arg, 'coverThreshold', 0.9, @(x) isnumeric(x) && (numel(x) == 1) && (x > 0) && (x <= 1));
addParamValue(arg, 'filter', ones(20,20), @(x) isnumeric(x));
addParamValue(arg, 'verbose', true, @(x) islogical(x) && (numel(x) == 1));
addParamValue(arg, 'fig', true, @(x) islogical(x) && (numel(x) == 1));

parse(arg, chm, refmat, varargin{:});


%% compute crown cover

chm(chm <= arg.Results.minHeight) = 0;
mask = single(chm > 0);

switch arg.Results.method
    
    case 'convolution'
        
        if arg.Results.verbose
            
            fprintf('computing convolution...');
            
        end
        
        % compute 2D convolution
        cover = imfilter(mask, arg.Results.filter, 'replicate') ./ nnz(arg.Results.filter);
        
        if arg.Results.verbose
            
            fprintf('done!\n');
            
        end
        
        if arg.Results.fig
            
            figure
            imagesc(cover)
            colormap(gray)
            title('Canopy cover ratio')
            xlabel('col')
            ylabel('row')
            axis equal tight
            colorbar
            
        end
        
    case 'triangulation'
        
        if isempty(arg.Results.stems)
            
            error('stem locations were not provided');
            
        end
        
        [nrows, ncols] = size(mask);
        col_stems = double(arg.Results.stems(:,1));
        row_stems = double(arg.Results.stems(:,2));
        n_trees = length(col_stems);
        
        % compute Delaunay triangulation
        if arg.Results.verbose
            
            fprintf('computing Delaunay triangulation...');
            
        end
        
        dt = delaunayTriangulation(col_stems, row_stems);
        
        if arg.Results.verbose
            
            fprintf('done!\n');
            
        end
        
        % determine pixel col/row coordinates
        [col_grid, row_grid] = meshgrid(1:ncols, 1:nrows);
        row_grid = reshape(row_grid, [], 1);
        col_grid = reshape(col_grid, [], 1);
        ind_grid = sub2ind([nrows, ncols], row_grid, col_grid);
        
        if arg.Results.verbose
            
            fprintf('computing triangulation label matrix...');
            
        end
        
        % compute triangulation label matrix
        id_triangulation = uint32(pointLocation(dt, [col_grid row_grid])); %
        
        if arg.Results.verbose
            
            fprintf('done!\n');
            
        end
        
        % compute triangle fill ratio
        if arg.Results.verbose
            
            fprintf('computing canopy / triangle area ratio...');
            
        end
        
        area_triangle_canopy = accumarray(id_triangulation + 1, reshape(chm > 0, 1, []));
        area_triangle_total = accumarray(id_triangulation + 1, true(size(id_triangulation)));
        canopy_fill_ratio = area_triangle_canopy ./ area_triangle_total;
        
        % set canopy / triangle area ratio threshold
        idxl_high_cover = (canopy_fill_ratio >= arg.Results.coverThreshold);
        
        if arg.Results.verbose
            
            fprintf('done!\n');
            
        end
        
        % compute crown triplet areas
        if ~isempty(arg.Results.crowns)
            
            if arg.Results.verbose
                
                fprintf('computing crown areas...');
                
            end
            
            labels_crowns = reshape(arg.Results.crowns,[],1);
            labels_crowns(sub2ind([nrows, ncols], row_stems, col_stems)) = 1:n_trees; % add marker pixel to label matrix
            
            area_crowns_chm = accumarray(labels_crowns + 1, reshape(chm > 0, 1, []));
            area_crowns = accumarray(labels_crowns + 1, true(size(labels_crowns)));
            
            % extract crown pixel indices
            [~, ind_sort] = sort(labels_crowns);
            col_grid_sort = col_grid(ind_sort);
            row_grid_sort = row_grid(ind_sort);
            col_pixel_crowns = mat2cell(col_grid_sort, area_crowns, 1);
            col_pixel_crowns = col_pixel_crowns(2:end,:);
            row_pixel_crowns = mat2cell(row_grid_sort, area_crowns, 1);
            row_pixel_crowns = row_pixel_crowns(2:end,:);
            
            % compute crown triplets areas
            area_crowns_chm = area_crowns_chm(2:end);
            area_crown_triplets = area_crowns_chm(dt.ConnectivityList);
            
            idxl_low_cover = ~idxl_high_cover(2:end);
            row_crown_triplets = cellfun(@(x1,x2,x3) [x1; x2; x3], row_pixel_crowns(dt.ConnectivityList(idxl_low_cover,1)), row_pixel_crowns(dt.ConnectivityList(idxl_low_cover,2)), row_pixel_crowns(dt.ConnectivityList(idxl_low_cover,3)), 'UniformOutput', 0);
            col_crown_triplets = cellfun(@(x1,x2,x3) [x1; x2; x3], col_pixel_crowns(dt.ConnectivityList(idxl_low_cover,1)), col_pixel_crowns(dt.ConnectivityList(idxl_low_cover,2)), col_pixel_crowns(dt.ConnectivityList(idxl_low_cover,3)), 'UniformOutput', 0);
            
            if arg.Results.verbose
                
                fprintf('done!\n');
                
            end
            
            % compute crown triplets convex hull
            if arg.Results.verbose
                
                fprintf('computing crown triplets convex hulls...');
                
            end
            
            [~, V] = cellfun(@convhull, col_crown_triplets, row_crown_triplets, 'UniformOutput', false); % subset
            
            if arg.Results.verbose
                
                fprintf('done!\n');
                
            end
            
            if arg.Results.verbose
                
                fprintf('computing crown cover...');
                
            end
            
            cover_crown_triplets = zeros(size(area_crowns));
            cover_crown_triplets(idxl_high_cover(2:end),1) = 1;
            cover_crown_triplets(idxl_low_cover,1) = sum(area_crown_triplets(idxl_low_cover,:),2) ./ [V{:}]'; % subset
            
            cover = zeros(size(chm), 'single');
            cover(ind_grid(id_triangulation > 0)) = cover_crown_triplets(id_triangulation(id_triangulation > 0));
            cover(cover > 1) = 1;
            
            if arg.Results.fig
                
                figure
                imagesc(cover)
                colormap(gray)
                title('Canopy cover ratio')
                xlabel('col')
                ylabel('row')
                axis equal tight
                colorbar
                
            end
            
            if arg.Results.verbose
                
                fprintf('done!\n');
                
            end
            
        end
        
end