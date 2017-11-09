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
%    [crh, xyh] = canopyPeaks(chm, ...
%        info.RefMatrix, ...
%        'method', 'allometricRadius', ...
%        'allometry', @(h) 1 + 0.5*log(max(h,1)), ...
%        'adjacency', @(h) min(0.5 + 0.5*log(max(h,1)),4), ...
%        'fig', true, ...
%        'verbose', true);
%
%    label = treeWatershed(chm, ...
%        'markers', crh, ...
%        'minHeight', 1, ...
%        'fig', true, ...
%        'verbose', true);
%
%    cover = canopyCover(chm, ...
%        info.RefMatrix, ...
%        'method', 'triangulation', ...
%        'stems', crh, ...
%        'crowns', label, ...
%        'verbose', true, ...
%        'fig', true);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2016b
%
% See also: canopyPeaks.m, treeWatershed.m
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory (lasig.epfl.ch)
% Website: http://mparkan.github.io/Digital-Forestry-Toolbox/
% Last revision: November 11, 2017
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund (WHFF, OFEV), project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'chm', @isnumeric);
addRequired(arg, 'refmat', @isnumeric);
addParameter(arg, 'method', 'convolution', @(x) ismember(x, {'convolution', 'triangulation'}));
addParameter(arg, 'stems', [], @(x) isnumeric(x) && size(x,1) >= 3);
addParameter(arg, 'crowns', [], @(x) isnumeric(x) && isequal(size(x), size(chm)));
addParameter(arg, 'minHeight', 1, @(x) isnumeric(x) && (numel(x) == 1));
addParameter(arg, 'coverThreshold', 0.9, @(x) isnumeric(x) && (numel(x) == 1) && (x > 0) && (x <= 1));
addParameter(arg, 'filter', ones(20,20), @(x) isnumeric(x));
addParameter(arg, 'verbose', true, @(x) islogical(x) && (numel(x) == 1));
addParameter(arg, 'fig', true, @(x) islogical(x) && (numel(x) == 1));

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
        
        if arg.Results.verbose
            
            fprintf('computing triangulation label matrix...');
            
        end
        
        % compute triangulation label matrix
        id_triangulation = uint32(pointLocation(dt, [col_grid row_grid])); %
        
        if arg.Results.verbose
            
            fprintf('done!\n');
            
        end
        
        % compute crown triplet areas
        if ~isempty(arg.Results.crowns)
            
            if arg.Results.verbose
                
                fprintf('computing crown triplet areas...');
                
            end
            
            dt_ind = sub2ind([nrows, ncols], dt.Points(:,2), dt.Points(:,1));
            dt_labels = arg.Results.crowns(dt_ind);
            dt_labels_triplet = dt_labels(dt.ConnectivityList);
            idxl_void = any(dt_labels_triplet == 0,2);
            dt_labels_triplet = dt_labels_triplet(~idxl_void,:);
            
            crown_area = accumarray(arg.Results.crowns(:) + 1, true(size(arg.Results.crowns(:))));
            crown_area = crown_area(2:end);
            crown_area_triplet = sum(crown_area(dt_labels_triplet), 2);
            
            if arg.Results.verbose
                
                fprintf('done!\n');
                
            end
            
            % compute crown triplets convex hull
            if arg.Results.verbose
                
                fprintf('computing crown triplets convex hull areas...');
                
            end
            
            convex_crown_area = zeros(size(crown_area_triplet,1),1);
            
            for j = 1:size(crown_area_triplet,1)
                
                [row_crown_triplets, col_crown_triplets] = ind2sub([nrows, ncols], find(ismember(arg.Results.crowns, dt_labels_triplet(j,:))));
                [~, convex_crown_area(j,1)] = convhull(col_crown_triplets, row_crown_triplets);
                
            end
            
            cover_crown_triplets = -ones(size(dt,1),1);
            cover_crown_triplets(~idxl_void) = crown_area_triplet ./ convex_crown_area;
            cover_crown_triplets = [-1; cover_crown_triplets];
            
            if arg.Results.verbose
                
                fprintf('done!\n');
                
            end
            
            cover = reshape(cover_crown_triplets(id_triangulation+1), nrows, ncols);
            
            if arg.Results.fig
                
                figure
                imagesc(cover)
                colormap(parula)
                title('Canopy cover ratio')
                xlabel('col')
                ylabel('row')
                axis equal tight
                colorbar
                
            end
            
        end
        
end