function [models, refmat] = elevationModels(xyz, classification, varargin)
% ELEVATIONMODELS - computes a Digital Terrain Model (DTM), a Digital Surface Model
% (DSM), a Digital Height Model (DHM) and point density maps from the 3D classified points defined by XYZ and CLASSIFICATION.
%
% [MODELS, REFMAT] = ELEVATIONMODELS(XYZ, CLASSIFICATION, ...) returns a structure containing
% elevation model values (MODELS.TERRAIN.VALUES, MODELS.SURFACE.VALUES, MODELS.HEIGHT.VALUES) and a boolean mask with
% missing data areas (MODELS.MASK) and point density maps (MODELS.DENSITY.OVERALL, MODELS.DENSITY.TERRAIN, MODELS.DENSITY.SURFACE)
% from the classified 3D point cloud defined by XYZ and CLASSIFICATION.
%
% Syntax:  [models, refmat] = elevationModels(xyz, classification, ...)
%
% Inputs:
%    xyz - Nx3 numeric matrix, containing the 3D point coordinates, i.e. [x, y, z]
%
%    classification - Nx1 numeric vector, containing point classification
%
%    cellSize (optional, default: 1) - numeric value, pixel (grid cell) size (in the same units as XY). If only cellSize is provided,
%    the upper left xy corner point is used as the center of the first
%    pixel. The cellSize argument is ignored if refmat is provided.
%
%    refmat (optional, default: []) - 3x2 numeric matrix, spatial referencing matrix, such that xy_map = [row, col, ones(nrows,1)] * refmat.
%    The cellSize argument is ignored if refmat is provided.
%
%    rasterSize (optional, default: []) - 1x2 numeric vector, target size [nrows, ncols] of the
%    output raster. Points that fall outside the extent of the raster are ignored.
%
%    classTerrain (optional, default: [2, 9]) - numeric vector, point classes used to represent the terrain
%
%    classSurface (optional, default: [3, 4, 5, 6]) - numeric vector, point
%    classes used to represent the surface (note that points with terrain class are also used)
%
%    interpolation (optional, default: 'idw') - interpolation method, 'linear', 'nearest', 'natural', 'idw'
%
%    searchRadius (optional, default: inf) -, inf, numeric value, search radius (in map units) used with the inverse distance 
%    weighting interpolation method
%    
%    weightFunction (optional, default: @(d) d^-3, @(x)) - anonymous function handle of the form @(d) = ..., 
%    weight function (in map units) used with the inverse distance weighting interpolation method
%
%    closing (optional, default: 1) - numeric value, the maximum distance (in map units) from data mask at which missing values (holes) are filled (i.e
%    interpolated/extrapolated). Use inf to fill all missing values.
%
%    smoothingFilter (optional, default: []) - numeric matrix, two-dimensional smoothing
%    filter, e.g. fspecial('gaussian', [4 4], 1)
%
%    outputModels (optional, default: {'terrain', 'surface', 'height', 'density'}) - cell array of strings, output models i.e {'terrain', 'surface', 'height', 'density'}
%
%    verbose (optional, default: true) - boolean value, verbosiy switch
%
%    fig (optional, default: true) - boolean value, switch to plot figures
%
% Outputs:
%    models - structure, boolean mask (data boundary), terrain, surface, height and point density models
%
%    refmat - 3x2 numeric matrix, spatial referencing matrix, such that xy_map = [row, col, ones(nrows,1)] * refmat
%
% Example:
%    pc = LASread('..\data\measurements\vector\als\so_2014_woodland_pasture.las', false, true);
%
%    [models, refmat] = elevationModels([pc.record.x, pc.record.y, pc.record.z], ...
%        pc.record.classification, ...
%        'classTerrain', [2], ...
%        'classSurface', [4,5], ...
%        'cellSize', 1, ...
%        'interpolation', 'idw', ...
%        'closing', 5, ...
%        'smoothingFilter', fspecial('gaussian', [3, 3], 0.5), ...
%        'outputModels', {'terrain', 'surface', 'height', 'density'}, ...
%        'fig', true, ...
%        'verbose', true);
%
%     spatialRef = refmatToMapRasterReference(refmat, size(models.height.values), ...
%       'rasterInterpretation', 'cells');
%
%     geotiffwrite(..\data\measurements\raster\als\so_2014_woodland_pasture.tif', ...
%       single(models.height.values), spatialRef, ...
%       'CoordRefSysCode', 'EPSG:21781');
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
% Last revision: February 21, 2019
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund (WHFF, OFEV), project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'xyz', @(x) (size(x,2) == 3) && isnumeric(x));
addRequired(arg, 'classification', @(x) (size(x,2) == 1) && isnumeric(x));
addParameter(arg, 'classTerrain', [2, 9], @isnumeric);
addParameter(arg, 'classSurface', [3, 4, 5, 6], @isnumeric);
addParameter(arg, 'cellSize', 1, @(x) isnumeric(x) & (numel(x) == 1));
addParameter(arg, 'refmat', [], @(x) (all(size(x) == [3, 2]) && isnumeric(x)) || isempty(x));
addParameter(arg, 'rasterSize', [], @(x) (isnumeric(x) && (numel(x) == 2) && all(x > 0)) || isempty(x));
addParameter(arg, 'closing', 1, @(x) isnumeric(x) & (numel(x) == 1) & (x >= 0));
addParameter(arg, 'smoothingFilter', [], @isnumeric);
addParameter(arg, 'interpolation', 'linear', @(x) ismember(x, {'nearest', 'natural', 'linear', 'idw'}));
addParameter(arg, 'searchRadius', inf, @(x) isnumeric(x) & (numel(x) == 1) & (x >= 1));
addParameter(arg, 'weightExponent', 3, @(x) isnumeric(x) & (numel(x) == 1) & (x >= 1));
addParameter(arg, 'weightFunction', @(d) d^-3, @(x) strfind(func2str(x),'@(d)') == 1);
addParameter(arg, 'outputModels', {'terrain', 'surface', 'height', 'density'}, @(x) iscell(x) & any(ismember(x, {'terrain', 'surface', 'height', 'density'})));
addParameter(arg, 'fig', true, @(x) islogical(x) && (numel(x) == 1));
addParameter(arg, 'verbose', true, @(x) islogical(x) && (numel(x) == 1));

parse(arg, xyz, classification, varargin{:});


%% setup spatial reference matrix

if isempty(arg.Results.refmat)
    
    dx = arg.Results.cellSize;
    dy = -arg.Results.cellSize;
    refmat = [0, dy; dx, 0; min(xyz(:,1))-dx, max(xyz(:,2))-dy];
    
else
    
    refmat = arg.Results.refmat;
    dx = refmat(2,1);
    dy = refmat(1,2);
    
end


%% convert map to image coordinates

P = [xyz(:,1) - refmat(3,1), xyz(:,2) - refmat(3,2)] / refmat(1:2,:); % center of upper left pixel -> [1,1]
row = round(P(:,1));
col = round(P(:,2));


if ~isempty(arg.Results.rasterSize)
    
    nrows = arg.Results.rasterSize(1);
    ncols = arg.Results.rasterSize(2);
    
else
    
    nrows = max(row);
    ncols = max(col);
    
end

%% compute boolean mask

if arg.Results.verbose
    
    fprintf('computing mask...');
    
end

[mask, refmat] = rasterize(xyz(:,1:2), ...
    xyz(:,3), ...
    'refmat', refmat, ...
    'rasterSize', [nrows, ncols], ...
    'fun', @any, ...
    'fill', false, ...
    'fig', false);

% fill holes in mask (with distance transform)
if arg.Results.closing > 0
    
    [D, ~] = bwdist(mask);
    mask = D <= ceil(arg.Results.closing / abs(dx));
    
end

models.mask = mask;

if arg.Results.verbose
    
    fprintf('done!\n');
    
end


%% compute overall point density

if any(ismember(arg.Results.outputModels, {'density'}))
    
    if arg.Results.verbose
        
        fprintf('computing overall (all classes) point density...');
        
    end
    
    [density_overall, refmat] = rasterize(xyz(:,1:2), ...
        xyz(:,3), ...
        'refmat', refmat, ...
        'rasterSize', [nrows, ncols], ...
        'fun', @numel, ...
        'fill', 0, ...
        'fig', false);
    
    models.density.overall = single(density_overall);
    models.density.overall(~models.mask) = nan;
    
    if arg.Results.verbose
        
        fprintf('done!\n');
        
    end
    
    % print overall (all classes) density statistics
    if arg.Results.verbose
        
        % compute point density statistics
        fprintf('median point density (all classes): %u\n', median(models.density.overall(models.mask)));
        fprintf('mean point density (all classes): %.1f\n', mean(models.density.overall(models.mask)));
        fprintf('max point density (all classes): %u\n', max(models.density.overall(models.mask)));
        fprintf('standard deviation point density (all classes): %.1f\n', std(models.density.overall(models.mask)));
        
    end
    
    % plot overall point density
    if arg.Results.fig
        
        figure
        gh01 = imagesc(models.density.overall);
        set(gh01, 'AlphaData', double(models.mask))
        title('Overall (all classes) point density')
        axis equal tight
        colorbar
        
    end
    
end


%% compute Terrain Model (DTM)

if any(ismember(arg.Results.outputModels, {'terrain', 'surface', 'height'}))
    
    idxl_terrain = ismember(classification, arg.Results.classTerrain);
    
    if any(idxl_terrain)
        
        [col_grid, row_grid] = meshgrid(1:ncols, 1:nrows);
        t = [row_grid(:), col_grid(:)] * refmat(1:2,:);
        x_grid = reshape(t(:,1) + refmat(3,1), nrows, ncols);
        y_grid = reshape(t(:,2) + refmat(3,2), nrows, ncols);
        
        OCTAVE_FLAG = (exist('OCTAVE_VERSION', 'builtin') ~= 0); % determine if system is Matlab or GNU Octave
        
        if OCTAVE_FLAG && ~strcmp(arg.Results.interpolation, 'idw')
            
            interpolation = 'idw';
            fprintf('\nWarning: Octave only supports inverse distance weighting interpolation\n');
            % Octave depends on the QHull library to compute the Delaunay triangulation required for
            % linear, nearest and natural neighbour interpolation. In its current version QHull downsamples 
            % the triangulation when the Qz option is specified (see https://github.com/qhull/qhull/issues/25) 
            % and is too slow for most practical purposes when the Qz option is disabled.
            
        else
            
            interpolation = arg.Results.interpolation;
            
        end
        
        if arg.Results.verbose
            
            fprintf('computing terrain model with %s interpolation...', interpolation);
            
        end

        switch interpolation
            
            case {'nearest', 'natural', 'linear'}
                
                % interpolation
                F = scatteredInterpolant(xyz(idxl_terrain,1), ...
                    xyz(idxl_terrain,2), ...
                    xyz(idxl_terrain,3), ...
                    interpolation, ...
                    'linear');
                
                dtm = F(x_grid, y_grid);
                
            case {'idw'}
                
                [dtm, ~] = rasterize(xyz(idxl_terrain,1:2), ...
                    xyz(idxl_terrain,3), ...
                    'refmat', refmat, ...
                    'rasterSize', [nrows, ncols], ...
                    'fun', @(x) sum(x)/numel(x), ...
                    'fill', nan, ...
                    'fig', false);
                
                % compute search radius in pixels 
                if isinf(arg.Results.searchRadius)
                    
                    r = ceil(max(max(bwdist(~isnan(dtm))))); % radius to fill all holes
                    
                else
                    
                    r = ceil(arg.Results.searchRadius ./ abs(dx));
                    
                end
                
                % create structuring element for convolution
                SE = double(bwdist(padarray(true, double([r,r]))));
                SE(r+1,r+1) = 0.001; % central element of convolution window
                SE(SE > r) = 0;
                
                % apply weight function
                SE = arrayfun(arg.Results.weightFunction, SE * abs(dx));
                SE(isinf(SE)) = 0;
                
                idxl_nan = isnan(dtm);
                dtm_i = dtm;
                dtm_i(idxl_nan) = 0;
                
                % sum of values
                V = imfilter(dtm_i, SE);
                
                % sum of weights
                W = imfilter(double(~idxl_nan), SE);
                
                % normalize
                dtm = V ./ W;
                
        end
        
        % extrapolation (nearest neighbour only)
        [~, idxn_nn] = bwdist(~isnan(dtm));
        dtm = dtm(idxn_nn);
    
        % smooth terrain model values
        if ~isempty(arg.Results.smoothingFilter)
            
            dtm = imfilter(dtm, arg.Results.smoothingFilter, 'replicate', 'same');
            
        else
            
            fprintf('\nWarning: no smoothing filter was applied to the terrain model\n');
            
        end
        
        % flip verticaly
        dtm = single(dtm);
        dtm(~models.mask) = nan;
        
        % store values and interpolant in structure
        models.terrain.values = dtm;
        
        % display terrain model plot
        if arg.Results.fig
            
            figure
            gh11 = imagesc(models.terrain.values);
            colormap(gray)
            set(gh11, 'AlphaData', double(models.mask)) % transparency not yet available in Octave;
            title('Terrain Model')
            axis equal tight
            colorbar
            
        end

        % compute terrain model point density
        if any(ismember(arg.Results.outputModels, {'density'}))
            
            [density_terrain, refmat] = rasterize(xyz(idxl_terrain,1:2), ...
                xyz(idxl_terrain,3), ...
                'refmat', refmat, ...
                'rasterSize', [nrows, ncols], ...
                'fun', @numel, ...
                'fill', 0, ...
                'fig', false);
            
            models.density.terrain = single(density_terrain);
            models.density.terrain(~models.mask) = nan;
            
            % display terrain model point density
            if arg.Results.fig
                
                figure
                gh12 = imagesc(models.density.terrain);
                set(gh12, 'AlphaData', double(models.mask))
                title('Terrain Model - Point density')
                axis equal tight
                colorbar
                
            end
            
        end
        
        if arg.Results.verbose
            
            fprintf('done!\n');
            
        end
        
    else
        
        error('There are no terrain points in point data record');
        
    end
    
end


%% compute Surface Model (DSM)

if any(ismember(arg.Results.outputModels, {'surface', 'height'}))
    
    if arg.Results.verbose
        
        fprintf('computing surface model...');
        
    end
    
    idxl_surface = ismember(classification, [arg.Results.classTerrain, arg.Results.classSurface]);
    
    [dsm, refmat] = rasterize(xyz(idxl_surface,1:2), ...
        xyz(idxl_surface,3), ...
        'refmat', refmat, ...
        'rasterSize', [nrows, ncols], ...
        'fun', @max, ...
        'fill', nan, ...
        'fig', false);
    
    % extrapolation (nearest neighbour only)
    [~, idxn_nn] = bwdist(~isnan(dsm));
    dsm = dsm(idxn_nn);
    
    if ~isempty(arg.Results.smoothingFilter)
        
        dsm = imfilter(dsm, arg.Results.smoothingFilter, 'replicate', 'same');
        
    else
        
        fprintf('\nWarning: no smoothing filter was applied to surface model\n');
        
    end
    
    dsm(~models.mask) = nan;
    models.surface.values = single(dsm);
    
    % display surface model
    if arg.Results.fig
        
        figure
        gh21 = imagesc(models.surface.values);
        set(gh21, 'AlphaData', double(models.mask))
        colormap(gray)
        title('Surface Model')
        axis equal tight
        colorbar
        
    end
    
    % compute surface model point density
    if any(ismember(arg.Results.outputModels, {'density'}))
        
        [density_surface, refmat] = rasterize(xyz(idxl_surface,1:2), ...
            xyz(idxl_surface,3), ...
            'refmat', refmat, ...
            'rasterSize', [nrows, ncols], ...
            'fun', @numel, ...
            'fill', 0, ...
            'fig', false);
        
        models.density.surface = single(density_surface);
        models.density.surface(~models.mask) = nan;
        
        % display surface model point density
        if arg.Results.fig
            
            figure
            gh22 = imagesc(models.density.surface);
            set(gh22, 'AlphaData', double(models.mask))
            title('Surface Model - Point density')
            axis equal tight
            colorbar
            
        end
        
    end
    
    if arg.Results.verbose
        
        fprintf('done!\n');
        
    end
    
end


%% compute Height Model (DHM)

if any(ismember(arg.Results.outputModels, 'height'))
    
    if arg.Results.verbose
        
        fprintf('computing height model...');
        
    end
    
    models.height.values = models.surface.values - models.terrain.values;
    models.height.values(models.height.values < 0) = 0;
    
    % display height model
    if arg.Results.fig
        
        figure
        gh31 = imagesc(models.height.values);
        set(gh31, 'AlphaData', double(models.mask))
        colormap(gray)
        title('Height Model')
        axis equal tight
        colorbar
        
    end
    
    if arg.Results.verbose
        
        fprintf('done!\n');
        
    end
    
end
