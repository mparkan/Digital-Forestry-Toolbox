function cover = canopyCover(chm, refmat, varargin)
% CANOPYCOVER - computes a canopy crown cover map
%
% COVER = CANOPYCOVER(CHM, REFMAT, ...) computes a canopy crown cover map using a 2D convolution or a variant of the Delaunay
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
%    label - RxC numeric matrix, individual tree crown label matrix with R rows
%    and C columns
%
%    minHeight (optional, default: 0.6) - numeric value, minimum canopy
%    height, all values below this threshold are set to zero.
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
%    pc = LASread('..\data\measurements\vector\als\so_2014_woodland_pasture.las');
%
%    [models, refmat] = elevationModels([pc.record.x, pc.record.y, pc.record.z], ...
%        pc.record.classification, ...
%        'classTerrain', 2, ...
%        'classSurface', 18, ...
%        'cellResolution', 1, ...
%        'closing', 10, ...
%        'smoothingFilter', fspecial('gaussian', [3, 3], 1), ...
%        'outputModels', {'terrain', 'surface', 'height'}, ...
%        'fig', true, ...
%        'verbose', true);
%
%    [crh, ~] = canopyPeaks(models.height.values, ...
%        refmat, ...
%        'method', 'allometricRadius', ...
%        'allometry', @(h) 0.5 + 0.25*log(max(h,1)), ...
%        'fig', true, ...
%        'verbose', true);
%
%    label = treeWatershed(models.height.values, ...
%        'markers', crh, ...
%        'minHeight', 1, ...
%        'fig', true, ...
%        'verbose', true);
%
%    cover_1 = canopyCover(models.height.values, ...
%        refmat, ...
%        'method', 'triangulation', ...
%        'stems', crh, ...
%        'label', label, ...
%        'verbose', true, ...
%        'fig', true);
%
%    cover_2 = canopyCover(models.height.values, ...
%        refmat, ...
%        'method', 'convolution', ...
%        'filter', ones(10,10), ...
%        'verbose', true, ...
%        'fig', true);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2017b, GNU Octave 4.2.2 (configured for "x86_64-w64-mingw32")
%
% See also: elevationModels, canopyPeaks.m, treeWatershed.m
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory (LASIG)
% Website: http://mparkan.github.io/Digital-Forestry-Toolbox/
% Last revision: April 12, 2018
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund (WHFF, OFEV), project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'chm', @isnumeric);
addRequired(arg, 'refmat', @isnumeric);
addParameter(arg, 'method', 'convolution', @(x) ismember(x, {'convolution', 'triangulation'}));
addParameter(arg, 'stems', [], @(x) isnumeric(x) && size(x,1) >= 3);
addParameter(arg, 'label', [], @(x) isnumeric(x) && isequal(size(x), size(chm)));
addParameter(arg, 'minHeight', 1, @(x) isnumeric(x) && (numel(x) == 1));
addParameter(arg, 'filter', ones(20,20), @(x) isnumeric(x));
addParameter(arg, 'verbose', true, @(x) islogical(x) && (numel(x) == 1));
addParameter(arg, 'fig', true, @(x) islogical(x) && (numel(x) == 1));

parse(arg, chm, refmat, varargin{:});

%% load Octave packages

OCTAVE_FLAG = (exist('OCTAVE_VERSION', 'builtin') ~= 0); % determine if system is Matlab or GNU Octave

if OCTAVE_FLAG
    
    pkg load statistics
    pkg load image
    more off % disable pager
    
end

%% compute crown cover

switch arg.Results.method
    
    case 'convolution'
        
        if arg.Results.verbose
            
            fprintf('computing convolution...');
            
        end
        
        chm(chm < arg.Results.minHeight) = 0;
        mask = single(chm > 0);

        % compute 2D convolution
        cover = imfilter(mask, arg.Results.filter, 'replicate') ./ nnz(arg.Results.filter);
        
        if arg.Results.verbose
            
            fprintf('done!\n');
            
        end
        
        if arg.Results.fig
            
            figure
            imagesc(cover)
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
        
        % remove crowns with height below threshold
        idxl_filter = arg.Results.stems(:,3) >= arg.Results.minHeight;
        label = arg.Results.label;
        label(~ismember(label, find(idxl_filter))) = 0;

        [nrows, ncols] = size(chm);
        col_stems = double(arg.Results.stems(idxl_filter,1));
        row_stems = double(arg.Results.stems(idxl_filter,2));
        
        % compute alpha triangulation
        if arg.Results.verbose
            
            fprintf('computing delaunay triangulation...');
            
        end
        
        col_corners = [1; ncols; ncols; 1];
        row_corners = [1; 1; nrows; nrows];
        
        col = [col_corners; col_stems];
        row = [row_corners; row_stems];
        
        dt = delaunay(col, row);
        
        if arg.Results.verbose
            
            fprintf('done!\n');
            
        end
        
        if arg.Results.verbose
            
            fprintf('computing triangulation label matrix...');
            
        end
        
        % determine pixel col/row coordinates
        [col_grid, row_grid] = meshgrid(1:ncols, 1:nrows);
        row_grid = row_grid(:);
        col_grid = col_grid(:);
        
        % closest triangulation simplex search
        id_triangulation = tsearchn([col, row], dt, [col_grid, row_grid]);
        
        if arg.Results.verbose
            
            fprintf('done!\n');
            
        end
        
        % compute crown triplet areas
        if ~isempty(label)
            
            if arg.Results.verbose
                
                fprintf('computing crown triplet areas...');
                
            end
            
            %dt_ind = sub2ind([nrows, ncols], row_stems, col_stems);
            dt_ind = sub2ind([nrows, ncols], row, col);
            
            % use distance transform to fill missing ITC labels
            [~, idxn_nnz] = bwdist(label ~= 0);
            itc_labels = label(idxn_nnz);
            dt_labels = itc_labels(dt_ind);
            
            % compute cumulative crown triplet area
            dt_labels_triplet = dt_labels(dt);
            
            crown_area = accumarray(label(:) + 1, true(size(label(:))));
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
                
                % convex_crown_area(j,1) =
                % nnz(bwconvhull(ismember(label, dt_labels_triplet(j,:)))); % Matlab only
                
                [row_crown_triplets, col_crown_triplets] = ind2sub([nrows, ncols], find(ismember(label, dt_labels_triplet(j,:))));
                [~, convex_crown_area(j,1)] = convhulln([col_crown_triplets, row_crown_triplets]);
                
            end
            
            cover_crown_triplets = min(crown_area_triplet ./ convex_crown_area, 1);
            
            if arg.Results.verbose
                
                fprintf('done!\n');
                
            end
            
            idxl_void = isnan(id_triangulation) | (id_triangulation == 0);
            cover = zeros(size(id_triangulation));
            cover(~idxl_void) = cover_crown_triplets(id_triangulation(~idxl_void));
            cover = reshape(cover, nrows, ncols);
            
            [d_nn, idxn_nn] = bwdist(cover ~= 0);
            cover(d_nn == 1) = cover(idxn_nn(d_nn == 1));
            
            if arg.Results.fig
                
                figure
                imagesc(cover)
                title('Canopy cover ratio')
                xlabel('col')
                ylabel('row')
                axis equal tight
                colorbar
                caxis([0,1])
                
            end
            
        end
        
end