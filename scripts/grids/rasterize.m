function varargout = rasterize(xyz, xv, yv, varargin)
% RASTERIZE - creates a 2D or 3D raster from scattered 3D points with coordinates XYZ by binning
% them in the cells defined by edges XV, YV and optionaly ZV, applying the FUN
% function to the VAL values contained in each grid cell. Cells which do not
% contain any points are filled with the FILL value.
%
% Syntax:  [idxl_in, sub_crl, raster] = rasterize(xyz, xv, yv, ...)
%
% Inputs:
%    xyz - Nx3 numeric matrix, xyz coordinates of the 3D points
%
%    xv - Cx1 numeric vector, containing the ordered x coordinates of the
%    grid cells column centers
%
%    yv - Rx1 numeric vector, containing the ordered y coordinates of the grid
%    cells row centers
%
%    zv (optional, default: []) - Lx1 numeric vector, containing the ordered z coordinates of the grid
%    cells level (stack) centers
%
%    val (optional, default: xyz(:,3)) - Nx1 numeric vector, values which will be
%    aggregated in each cell (defaults to z)
%
%    fun (Optional) - anonymous function handle, function used to aggregate the values contained in
%    each grid cell (e.g. max, min, std, mode, etc)
%
%    fill (Optional) - numeric value, fill value attributed to cells containing no data
%
% Outputs:
%    idxl_in - Nx1 boolean vector, flag indicating if each input point is within
%    the defined grid extent 
%
%    sub_crl - Nx3 numeric matrix, column, row, level subscripts of each point
%
%    raster (optional) - RxC or RxCxL numeric matrix, gridded values aggregated
%    with the specified function handle
%
% Example:
%    pc = LASread('..\data\measurements\vector\als\zh_2014_coniferous.las', false, true);
%    xyz = [pc.record.x, pc.record.y, pc.record.z];
%    
%    dxy = 0.5;
%    XY = round(xyz(:,1:2) / dxy) * dxy;
%    xy_min = min(XY, [], 1);
%    xy_max = max(XY, [], 1);
%    xv = xy_min(1):dxy:xy_max(1);
%    yv = xy_min(2):dxy:xy_max(2);
%
%    [~, ~, raster] = rasterize(xyz, ...
%        xv, ...
%        yv, ...
%        [], ...
%        xyz(:,3), ...
%        @(x) numel(x));
%    
%    figure
%    imagesc(raster)
%    colorbar  
%    axis equal tight
%
% Other m-files required: none
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
addRequired(arg, 'xv', @(x) ismatrix(x) && isnumeric(x));
addRequired(arg, 'yv', @(x) ismatrix(x) && isnumeric(x));
addOptional(arg, 'zv', [], @(x) ismatrix(x) && isnumeric(x));
addOptional(arg, 'val', [], @(x) isnumeric(x) || islogical(x));
addOptional(arg, 'fun', [], @(x) isa(x, 'function_handle'));
addOptional(arg, 'fill', [], @(x) (numel(x) == 1) && (isnumeric(x) || islogical(x)));

parse(arg, xyz, xv, yv, varargin{:});

flag_raster = any(~isempty(arg.Results.val) || ~isempty(arg.Results.fun) || ~isempty(arg.Results.fill));


%% reformat input

x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);

dx = abs(xv(2) - xv(1));
dy = abs(yv(2) - yv(1));

xv = linspace(xv(1)-dx/2, xv(end)+dx/2, length(xv)+1);
yv = linspace(yv(1)-dy/2, yv(end)+dy/2, length(yv)+1);

    
%% rasterize

if isempty(arg.Results.zv) % 2D raster
    
    % check if coordinates are located within the grid extent
    idxl_in = (x >= xv(1)) & (x <= xv(end)) & (y >= yv(1)) & (y <= yv(end));
    
    varargout{1} = idxl_in;
    
    if ~all(idxl_in)
        
        fprintf('\nWarning: %u points are located outside the defined grid extent\n', nnz(~idxl_in));
        
    end
    
    % pixel correspondance indices
    [~, ind_col] = histc(x(idxl_in), xv);
    [~, ind_row] = histc(y(idxl_in), yv);
    
    varargout{2} = [ind_col, ind_row];
    
    if flag_raster

        ncols = length(xv)-1;
        nrows = length(yv)-1;
        varargout{3} = accumarray([ind_row, ind_col], ...
            arg.Results.val(idxl_in), ...
            [nrows, ncols], ...
            arg.Results.fun, ...
            arg.Results.fill); % raster
        
    end
    
else % 3D raster
    
    zv = arg.Results.zv;
    dz = abs(zv(2) - zv(1));
    zv = linspace(zv(1)-dz/2, zv(end)+dz/2, length(zv)+1);
    
    idxl_in = (x >= xv(1) & x <= xv(end) & y >= yv(1) & y <= yv(end) & z >= zv(1) & z <= zv(end));
    varargout{1} = idxl_in;
    
    if ~all(idxl_in)
        
        fprintf('\nWarning: %u points are located outside the defined grid extent\n', nnz(~idxl_in));
        
    end
    
    [~, ind_col] = histc(x(idxl_in), xv);
    [~, ind_row] = histc(y(idxl_in), yv);
    [~, ind_lev] = histc(z(idxl_in), zv);
    
    varargout{2} = uint32([ind_col, ind_row, ind_lev]);
    
    if flag_raster
        
        ncols = length(xv)-1;
        nrows = length(yv)-1;
        nlevs = length(zv)-1;
        varargout{3} = accumarray([ind_row, ind_col, ind_lev], ...
            arg.Results.val(idxl_in), ...
            [nrows, ncols, nlevs], ...
            arg.Results.fun, ...
            arg.Results.fill); % raster
        
    end
    
end
