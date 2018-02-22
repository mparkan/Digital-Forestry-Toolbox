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
%    cellResolution (optional, default: 1) - numeric value, xy pixel
%    spatial resolution (i.e. grid cell size) of the output raster models
%  
%    xv (optional, default: []) - Cx1 numeric matrix, containing the ordered x coordinates of the
%    grid cells column centers
%
%    yv (optional, default: []) - Rx1 matrix, containing the ordered y coordinates of the grid
%    cells row centers
%
%    classTerrain (optional, default: [2]) - numeric vector, point classes used to represent the terrain
%
%    classSurface (optional, default: [4, 5, 6]) - numeric vector, point classes used to represent the surface
%
%    closing (optional, default: 0) - integer value, the number of morphological closing operations applied to 
%    the mask to suppress holes in the raster
%
%    smoothingFilter (optional, default: []) - numeric matrix, two-dimensional smoothing
%    filter, e.g. fspecial('gaussian', [4 4], 1)
%
%    alpha (optional, default: 5) - numeric value, minimum edge length (in map units, not pixel units) used
%    when computing the alpha shape  that defines the data mask
%
%    outputModels (optional, default: {'terrain', 'surface', 'height', 'density'}) - cell array of strings, output models i.e {'terrain', 'surface', 'height', 'density'}
%
%    verbose (optional, default: true) - boolean value, verbosiy switch
%
%    fig (optional, default: true) - boolean value, switch to plot figures
%
% Outputs:
%    models - structure, boolean mask (data boundary), terrain,
%    surface, height and point density models
%
%    refmat - 3x2 numeric matrix, spatial referencing matrix, such that xy_map = [row, col, ones(nrows,1)] * refmat
%
% Example:
%    pc = LASread('..\data\measurements\vector\als\so_2014_woodland_pasture.las', false, true);
%
%    [models, refmat] = elevationModels([pc.record.x, pc.record.y, pc.record.z], ...
%     pc.record.classification, ...
%     'classTerrain', [2], ...
%     'classSurface', [4,5], ...
%     'cellResolution', 1, ...
%     'closing', 3, ...
%     'smoothingFilter', fspecial('gaussian', [3, 3], 0.5), ...
%     'alpha', 5, ...
%     'outputModels', {'terrain', 'surface', 'height', 'density'}, ...
%     'fig', true, ...
%     'verbose', true);
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
% Compatibility: tested on Matlab R2017b, GNU Octave 4.2.1 (configured for "x86_64-w64-mingw32")
%
% See also: rasterize.m
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory (LASIG)
% Website: http://mparkan.github.io/Digital-Forestry-Toolbox/
% Last revision: February 22, 2018
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund (WHFF, OFEV), project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'xyz', @(x) (size(x,2) == 3) && isnumeric(x));
addRequired(arg, 'classification', @(x) (size(x,2) == 1) && isnumeric(x));
addParameter(arg, 'classTerrain', [2, 9], @isnumeric);
addParameter(arg, 'classSurface', [4, 5, 6], @isnumeric);
addParameter(arg, 'cellResolution', 1, @(x) isnumeric(x) & (numel(x) == 1));
addParameter(arg, 'xv', [], @isnumeric);
addParameter(arg, 'yv', [], @isnumeric);
addParameter(arg, 'closing', inf, @(x) isnumeric(x) & (numel(x) == 1) & (x >= 0));
addParameter(arg, 'smoothingFilter', [], @isnumeric);
addParameter(arg, 'outputModels', {'terrain', 'surface', 'height', 'density'}, @(x) iscell(x) & any(ismember(x, {'terrain', 'surface', 'height', 'density'})));
addParameter(arg, 'fig', true, @(x) islogical(x) && (numel(x) == 1));
addParameter(arg, 'verbose', true, @(x) islogical(x) && (numel(x) == 1));

parse(arg, xyz, classification, varargin{:});


%% create regular grid

if isempty(arg.Results.xv) || isempty(arg.Results.yv)
    
    dx = arg.Results.cellResolution;
    dy = arg.Results.cellResolution;

    XY = round(xyz(:,1:2) / arg.Results.cellResolution) * arg.Results.cellResolution;
    xy_min = min(XY, [], 1);
    xy_max = max(XY, [], 1);
    xv = xy_min(1):arg.Results.cellResolution:xy_max(1);
    yv = xy_min(2):arg.Results.cellResolution:xy_max(2);

else
    
    xv = arg.Results.xv;
    yv = arg.Results.yv;
    dx = abs(xv(1) - xv(2));
    dy = abs(yv(1) - yv(2));
    
end

% coordinates of cell (pixel) centers
[grid_x, grid_y] = meshgrid(xv, yv);


%% compute boolean mask

if arg.Results.verbose
    
    fprintf('computing mask...');
    
end

[~, ~, mask] = rasterize(xyz, xv, yv, [], xyz(:,3), @any, false);
n = ceil(arg.Results.closing);
if n > 0
    
    BW = bwmorph(padarray(mask, [n, n], false), 'close', n);
    mask = BW(n+1:end-n,n+1:end-n);
    
end

models.mask = flipud(mask);

if arg.Results.verbose
    
    fprintf('done!\n');
    
end
    

%% compute overall point density

if any(ismember(arg.Results.outputModels, {'density'}))
    
    if arg.Results.verbose
        
        fprintf('computing overall (all classes) point density...');
        
    end
    
    [~, ~, density_overall] = rasterize(xyz, xv, yv, [], xyz(:,3), @numel, 0);
    models.density.overall = single(flipud(density_overall));
    
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
        
        if arg.Results.verbose
            
            fprintf('computing terrain model...');
            
        end
        
        xyz_terrain = xyz(idxl_terrain, :);
        
        % interpolation
        dtm = griddata(xyz_terrain(:,1), ...
            xyz_terrain(:,2), ...
            xyz_terrain(:,3), ...
            grid_x, ...
            grid_y, ...
            'linear');
        
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
        dtm = single(flipud(dtm));
        
        % store values and interpolant in structure
        models.terrain.values = dtm;

        % display terrain model plot
        if arg.Results.fig
            
            figure
            gh11 = imagesc(models.terrain.values);
            colormap(gray)
            set(gh11, 'AlphaData', double(models.mask))
            title('Terrain Model')
            axis equal tight
            colorbar
            
        end
        
        % compute terrain model point density
        if any(ismember(arg.Results.outputModels, {'density'}))
            
            [~, ~, density_terrain] = rasterize(xyz_terrain, xv, yv, [], xyz_terrain(:,3), @numel, 0);
            models.density.terrain = single(flipud(density_terrain));
            
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
    
    xyz_surface = xyz(ismember(classification, [arg.Results.classTerrain, arg.Results.classSurface]), :);
    
    % rasterize
    [idxl_in, ~, dsm] = rasterize(xyz_surface, xv, yv, [], xyz_surface(:,3), @max, NaN);
    dsm = flipud(dsm); % flip verticaly
    idxl_dsm_missing = isnan(dsm);

    % interpolate missing areas with terrain model values
    dsm(idxl_dsm_missing) = dtm(idxl_dsm_missing);
    
    if ~isempty(arg.Results.smoothingFilter)
        
        dsm = imfilter(dsm, arg.Results.smoothingFilter, 'replicate', 'same');
        
    else
        
        fprintf('\nWarning: no smoothing filter was applied to surface model\n');
        
    end
        
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
        
        [~, ~, density_surface] = rasterize(xyz_surface(idxl_in,:), xv, yv, [], xyz_surface(idxl_in,3), @numel, 0);
        models.density.surface = single(flipud(density_surface));
        
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


%% create spatial referencing matrix

refmat = [0, -dy; dx, 0; grid_x(1,1)-dx, grid_y(end,1)+dy];
