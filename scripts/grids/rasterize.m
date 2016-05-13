function varargout = rasterize(xyz, xv, yv, varargin)
%RASTERIZE - creates a 2D or 3D raster from scattered points with coordinates X,Y,Z by binning
%them in the cells defined by XV, YV and optionaly ZV, applying the FUN
%function to the VAL values contained in each grid cell.
%
% Syntax:  [SUB, ...] = rasterize(xyz, xv, yv, ...)
%
% Inputs:
%    xyz - A nx3 vector containing the xyz coordinates
%
%    xv - A numeric vector containing the ordered x (column) coordinates of the grid cells
%
%    yv - A numeric vector containing the ordered y (row) coordinates of the grid cells
%
%    zv (optional, default: []) - A numeric vector containing the ordered z (stack) coordinates of the grid cells
%
%    val (optional, default: xyz(:,3)) - A nx1 vector containing the values which will be
%    aggregated in each cell (defaults to z)
%
%    fun (Optional) - An anonymous function used to aggregate the values contained in
%    each grid cell (e.g. max, min, std, mode, etc)
%
% Outputs:
%    sub - A Nx3 matrix containing the column, row, level subscripts of
%    each point within the raster
%
%    raster (Optional) - A 2D or 3D sparse matrix containing the gridded values aggregated
%    with the specified function handle
%
% Example:
%    pc = LASread('..\data\measurements\vector\als\zh_6995_2710_coniferous.las', false, true);
%    x = pc.record.x;
%    y = pc.record.y;
%    z = pc.record.z;
%      
%    dxy = 0.5;
%    xv = min(x)-dxy:dxy:max(x)+dxy;
%    yv = min(y)-dxy:dxy:max(y)+dxy;
%    zv = min(z)-dxy:dxy:max(z)+dxy;
%    [sub, raster] = rasterize([x, y, z], xv, yv, zv, z, @(x) numel(x));
%    
%    nrows = length(yv);
%    ncols = length(xv);
%    nlevels = length(zv);
%     
%    ind = sub2ind([nrows, ncols, nlevels], sub(:,2), sub(:,1), sub(:,3)); % 3d (row, col, lev) linear index for each point
%    [idxn_cells, ~, idxn_voxel_to_point] = unique(ind); % assign raster cell metrics to corresponding points
%    [row_voxel, col_voxel, lev_voxel] = ind2sub([nrows, ncols, nlevels], idxn_cells);
%    density = raster(idxn_cells);
%     
%    figure
%    scatter3(col_voxel, row_voxel, lev_voxel, 20, ...
%    density, ...
%    'Marker', '.');
%    xlabel('col')
%    ylabel('row')
%    zlabel('level')
%    title('point count per voxel')
%    colorbar
%    axis equal tight vis3d
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2016a
%
% See also:
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory
% Website: http://lasig.epfl.ch/
% Last revision: May 13, 2016
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund (WHFF, OFEV), project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'xyz', @(x) (size(x,2) == 3) && isnumeric(x));
addRequired(arg, 'xv', @(x) ismatrix(x) && isnumeric(x));
addRequired(arg, 'yv', @(x) ismatrix(x) && isnumeric(x));
addOptional(arg, 'zv', [], @(x) ismatrix(x) && isnumeric(x));
addOptional(arg, 'val', [], @(x) isnumeric(x) || islogical(x));
addOptional(arg, 'fun', [], @(x) isa(x, 'function_handle'));

parse(arg, xyz, xv, yv, varargin{:});

flag_raster = any(~isempty(arg.Results.val) || ~isempty(arg.Results.fun)));


%% reformat input

x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);


%% rasterize

if isempty(arg.Results.zv) % 2D raster
    
    % check if coordinates are located within the grid extent
    idx_extent = (x >= xv(1) & x <= xv(end) & y >= yv(1) & y <= yv(end));
    
    if ~all(idx_extent)
        
        error('some points are located outside the defined grid extent');
        
    end
    
    % pixel correspondance indices
    [~, ind_col] = histc(x, xv);
    [~, ind_row] = histc(y, yv);
    
    varargout{1} = [ind_col, ind_row]; % crl
    
    if flag_raster

        ncols = length(xv);
        nrows = length(yv);
        ind_cr = sub2ind([nrows, ncols], ind_row, ind_col);
        
        varargout{2} = accumarray(ind_cr, arg.Results.val, nrows*ncols, arg.Results.fun, 0, true); % raster
        
    end
    
else % 3D raster
    
    zv = arg.Results.zv;
    idxl_extent = (x >= xv(1) & x <= xv(end) & y >= yv(1) & y <= yv(end) & z >= zv(1) & z <= zv(end));
    
    if ~all(idxl_extent)
        
        error('some points are located outside the defined grid extent');
        
    end
    
    [~, ind_col] = histc(x, xv);
    [~, ind_row] = histc(y, yv);
    [~, ind_lev] = histc(z, zv);

    varargout{1} = uint32([ind_col, ind_row, ind_lev]);
    
    if flag_raster
        
        ncols = length(xv);
        nrows = length(yv);
        nlevs = length(zv);
        ind_crl = sub2ind([nrows, ncols, nlevs], ind_row, ind_col, ind_lev);
     
        varargout{2} = accumarray(ind_crl, arg.Results.val, [nrows*ncols*nlevs, 1], arg.Results.fun, 0, true); % raster

    end
    
end
