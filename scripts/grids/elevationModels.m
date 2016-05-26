function [models, refmat] = elevationModels(xyz, classification, varargin)
% ELEVATIONMODELS - computes a Digital Terrain Model (DTM), a Digital Surface Model
% (DSM) and a Digital Height Model (DHM) from the 3D classified points defined by XYZ and CLASSIFICATION.
%
% [MODELS, REFMAT] = ELEVATIONMODELS(XYZ, CLASSIFICATION, ...) returns a structure containing 
% elevation model values (MODELS.TERRAIN.VALUES, MODELS.SURFACE.VALUES and MODELS.HEIGHT.VALUES), their respective interpolant 
% functions (MODELS.TERRAIN.INTERPOLANT) and optionally point density maps (MODELS.TERRAIN.POINT_DENSITY, MODELS.SURFACE.POINT_DENSITY) 
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
%    classTerrain (optional, default: [2, 3, 9]) - numeric vector, point classes used to represent the terrain
%
%    classSurface (optional, default: [4, 5, 6]) - numeric vector, point classes used to represent the surface
%
%    maxFillArea (optional, default: inf) - numeric value, maximum contiguous missing data area that will
%    be interpolated. Missing data areas larger than this value will be filled with NaN.
%
%    smoothingFilter (optional, default: []) - numeric matrix, two-dimensional smoothing
%    filter, e.g. fspecial('gaussian', [4 4], 1)
%
%    outputModels (optional, default: {'dtm','dsm','dhm'}) - cell array of strings, output models i.e {'dtm','dsm','dhm'}
%
%    outputDensity (optional, default: true) - boolean value, flag indicating if point density maps are computed
%
%    verbose (optional, default: true) - boolean value, verbosiy switch
%
%    fig (optional, default: true) - boolean value, switch to plot figures
%
% Outputs:
%    models - structure, terrain, surface and height models with their respective
%    interpolant functions and optionally point density maps
%
%    refmat - 3x2 numeric matrix, spatial referencing matrix, such that xy_map = [row, col, ones(nrows,1)] * refmat
%
% Example:
%    pc = LASread('..\data\measurements\vector\als\so_2014_woodland_pasture.las', false, true);
%    x = pc.record.x;
%    y = pc.record.y;
%    z = pc.record.z;
%    classification = pc.record.classification;
%
%    [models, refmat] = elevationModels([x, y, z], classification, ...
%     'classTerrain', [2], ...
%     'classSurface', [18], ...
%     'cellResolution', 1, ...
%     'maxFillArea', inf, ...
%     'smoothingFilter', [], ...
%     'outputModels', {'dtm','dsm','dhm'}, ...
%     'outputDensity', true, ...
%     'fig', true, ...
%     'verbose', true);
%
% Other m-files required: rasterize.m
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2016a
%
% See also: rasterize.m
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory
% Website: http://lasig.epfl.ch/
% Last revision: May 26, 2016
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund (WHFF, OFEV), project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'xyz', @(x) (size(x,2) == 3) && isnumeric(x));
addRequired(arg, 'classification', @(x) (size(x,2) == 1) && isnumeric(x));
addParamValue(arg, 'classTerrain', [2, 3, 9], @isnumeric);
addParamValue(arg, 'classSurface', [4, 5, 6], @isnumeric);
addParamValue(arg, 'cellResolution', 1, @(x) isnumeric(x) & (numel(x) == 1));
addParamValue(arg, 'xv', [], @isnumeric);
addParamValue(arg, 'yv', [], @isnumeric);
addParamValue(arg, 'maxFillArea', inf, @(x) isnumeric(x) & (numel(x) == 1));
addParamValue(arg, 'smoothingFilter', [], @isnumeric); % fspecial('gaussian', [6, 6], 1)
addParamValue(arg, 'outputModels', {'dtm', 'dsm', 'dhm'}, @(x) iscell(x) & any(ismember(x, {'dtm', 'dsm', 'dhm'})));
addParamValue(arg, 'outputDensity', false, @(x) islogical(x) & (numel(x) == 1));
addParamValue(arg, 'fig', true, @(x) islogical(x) && (numel(x) == 1));
addParamValue(arg, 'verbose', true, @(x) islogical(x) && (numel(x) == 1));

parse(arg, xyz, classification, varargin{:}); % classTerrain, classSurface


%% create regular grid

if isempty(arg.Results.xv) || isempty(arg.Results.yv)
    
    dx = arg.Results.cellResolution;
    dy = arg.Results.cellResolution;
    
    x_min = min(xyz(:,1));
    x_max = max(xyz(:,1));
    y_min = min(xyz(:,2));
    y_max = max(xyz(:,2));
    
    xv = x_min:dx:x_max;
    yv = y_min:dy:y_max;
    
    xv = [xv, xv(end)+ceil(mod(x_max, dx))*dx];
    yv = [yv, yv(end)+ceil(mod(y_max, dy))*dy];

else
    
    xv = arg.Results.xv;
    yv = arg.Results.yv;
    dx = abs(xv(1) - xv(2));
    dy = abs(yv(1) - yv(2));
    
end

[grid_x, grid_y] = meshgrid(xv, yv);
[nrows, ncols] = size(grid_x);


%% compute Terrain Model (DTM)

if any(ismember(arg.Results.outputModels, {'dtm', 'dsm', 'dhm'}));
    
    idxl_terrain = ismember(classification, arg.Results.classTerrain);
    
    if any(idxl_terrain)
        
        if arg.Results.verbose
            
            fprintf('computing terrain model...');
            
        end
        
        xyz_terrain = xyz(idxl_terrain, :);
        
        % rasterize
        [idxl_in, dtm_cr, dtm] = rasterize(xyz_terrain(:,1:3), xv, yv, [], xyz_terrain(:,3), @min, NaN);
        
        % interpolate missing areas
        idxl_dtm_missing = isnan(dtm);
        idxl_dtm_valid = ~idxl_dtm_missing;
        dtm_connected_components = bwconncomp(idxl_dtm_missing, 8);
        dtm_component_areas = dx * dy * cellfun(@numel, dtm_connected_components.PixelIdxList, 'UniformOutput', true)';
        
        idxn_nan = cell2mat(dtm_connected_components.PixelIdxList(dtm_component_areas > arg.Results.maxFillArea)');
        
        interpolant_dtm = scatteredInterpolant(grid_x(idxl_dtm_valid), grid_y(idxl_dtm_valid), dtm(idxl_dtm_valid), ...
            'linear', 'nearest');
        dtm = interpolant_dtm(grid_x, grid_y);
        
        % smooth terrain model values
        if ~isempty(arg.Results.smoothingFilter)
            
            dtm = imfilter(dtm, arg.Results.smoothingFilter, 'replicate', 'same');
            
        else
            
            fprintf('\nWarning: no smoothing filter was applied to the terrain model\n');
            
        end
        
        % NaN fill in missing areas
        dtm(idxn_nan) = NaN;
        
        % flip verticaly
        dtm = flipud(dtm);
        
        % store values and interpolant in structure
        models.terrain.values = dtm;
        models.terrain.interpolant = interpolant_dtm;

        % display terrain model plot
        if arg.Results.fig
            
            figure
            imagesc(models.terrain.values)
            colormap(gray)
            title('Terrain Model')
            axis equal tight
            colorbar
            
        end
        
        % compute terrain model point density
        if arg.Results.outputDensity
            
            %[nrows, ncols] = size(dtm);
            models.terrain.point_density = accumarray([dtm_cr(:,2), dtm_cr(:,1)], xyz_terrain(idxl_in,1), [nrows, ncols], @numel, 0);
            
            % display terrain model point density
            if arg.Results.fig
                
                figure
                imagesc(models.terrain.point_density)
                title('TerrainModel - Point density')
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

if any(ismember(arg.Results.outputModels, {'dsm', 'dhm'}));
    
    if arg.Results.verbose
        
        fprintf('computing surface model...');
        
    end
    
    xyz_surface = xyz(ismember(classification, [arg.Results.classTerrain, arg.Results.classSurface]), :);
    
    % rasterize
    [idxl_in, dsm_cr, dsm] = rasterize(xyz_surface, xv, yv, [], xyz_surface(:,3), @max, NaN);
    dsm = flipud(dsm); % flip verticaly
    idxl_dsm_missing = isnan(dsm);

    if ~isempty(arg.Results.smoothingFilter)
        
        dsm = imfilter(dsm, arg.Results.smoothingFilter, 'replicate', 'same');
        
    else
        
        fprintf('\nWarning: no smoothing filter was applied to surface model\n');
        
    end
    
    % interpolate missing areas with terrain model values
    dsm(idxl_dsm_missing) = dtm(idxl_dsm_missing);
    
    models.surface.values = dsm;
    
    % display surface model
    if arg.Results.fig
        
        figure
        imagesc(models.surface.values)
        colormap(gray)
        title('Surface Model')
        axis equal tight
        colorbar
        
    end
    
    % compute surface model point density
    if arg.Results.outputDensity
        
        models.surface.point_density = accumarray([dsm_cr(:,2), dsm_cr(:,1)], xyz_surface(idxl_in,1), [nrows, ncols], @numel, 0); 
        
        % display surface model point density
        if arg.Results.fig
            
            figure
            imagesc(models.surface.point_density)
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

if any(ismember(arg.Results.outputModels, 'dhm'));
    
    if arg.Results.verbose
        
        fprintf('computing height model...');
        
    end
    
    models.height.values = models.surface.values - models.terrain.values;
    
    % display height model
    if arg.Results.fig
        
        figure
        imagesc(models.height.values)
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

refmat = [0, -dy; dx, 0; grid_x(1,1) - dx/2, grid_y(end,1) + dy/2];
%refmat = makerefmat(grid_x(1,1) + dx/2, grid_y(end,1) - dy/2, dx, -dy)