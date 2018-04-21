function varargout = rasterize(xy, v, varargin)
% RASTERIZE - creates a 2D raster from scattered points with coordinates XY and values V by
% aggregating the values in each grid cell with function FUN. Cells which do not
% contain any points are filled with the FILL value. The dimensions of the
% raster are determined by the spatial reference matrix REFMAT or the pixel size CELLSIZE.
%
% Syntax:  [I, refmat] = rasterize(xy, v, ...)
%
% Inputs:
%    xy - Nx3 numeric matrix, xy coordinates of the points
%
%    v - Nx1 numeric vector, values which will be aggregated in each cell
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
%    fun (optional, default: @max) - anonymous function handle, function used to aggregate the values contained in
%    each grid cell (e.g. max, min, std, mode, etc)
%
%    fill (optional, default: nan) - numeric or logical value, fill value attributed to cells containing no data
%
%    fig (optional, default: true) - boolean value, switch to plot figures
%
% Outputs:
%    I - RxC matrix, rasterized values aggregated with the specified function handle
%
%    refmat - 3x2 numeric matrix, spatial referencing matrix, such that xy_map = [row, col, ones(nrows,1)] * refmat
%
%    ind - Nx1 numeric vector, pixel index of each xy point
%
% Example:
%    pc = LASread('..\data\measurements\vector\als\zh_2014_coniferous.las', false, true);
%    
%    % create a spatial reference matrix
%    dx = 0.5;
%    dy = -0.5;
%    refmat = [0, dy; dx, 0; min(pc.record.x)-dx, max(pc.record.y)-dy];
%
%    % find the maximum elevation in each 50cm x 50 cm grid cell
%    z_max = rasterize([pc.record.x, pc.record.y], ...
%         pc.record.z, ...
%         'refmat', refmat, ...
%         'fun', @(x) max(x), ...
%         'fill', nan, ...
%         'fig', true);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2017b, GNU Octave 4.2.2 (configured for "x86_64-w64-mingw32")
%
% See also:
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory (LASIG)
% Website: http://mparkan.github.io/Digital-Forestry-Toolbox/
% Last revision: April 21, 2018
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund (WHFF, OFEV), project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'xy', @(x) (size(x,2) == 2) && isnumeric(x));
addRequired(arg, 'v', @(x) (size(x,2) == 1) && isnumeric(x));
addParameter(arg, 'cellSize', 1, @(x) isnumeric(x) & (numel(x) == 1));
addParameter(arg, 'refmat', [], @(x) (all(size(x) == [3, 2]) && isnumeric(x)) || isempty(x));
addParameter(arg, 'rasterSize', [], @(x) (isnumeric(x) && (numel(x) == 2) && all(x > 0)) || isempty(x));
addParameter(arg, 'fun', @max, @(x) isa(x, 'function_handle'));
addParameter(arg, 'fill', nan, @(x) (isnumeric(x) || islogical(x)) & (numel(x) == 1));
addParameter(arg, 'fig', false, @(x) islogical(x) && (numel(x) == 1));
addParameter(arg, 'verbose', true, @(x) islogical(x) && (numel(x) == 1));

parse(arg, xy, v, varargin{:});

nargoutchk(0, 3);

%% setup spatial reference matrix

if isempty(arg.Results.refmat)
    
    dx = arg.Results.cellSize;
    dy = -arg.Results.cellSize;
    refmat = [0, dy; dx, 0; min(xy(:,1))-dx, max(xy(:,2))-dy];

else
    
    refmat = arg.Results.refmat;
    
end


%% convert map to image coordinates

P = [xy(:,1) - refmat(3,1), xy(:,2) - refmat(3,2)] / refmat(1:2,:); % center of upper left pixel -> [1,1]
row = round(P(:,1));
col = round(P(:,2));

if ~isempty(arg.Results.rasterSize)
    
    nrows = arg.Results.rasterSize(1);
    ncols = arg.Results.rasterSize(2);
    
else
    
    nrows = max(row);
    ncols = max(col);
    
end


%% clip values that fall outside the extent of the raster grid

idxl_roi = (row > 0) & (row <= nrows) & (col > 0) & (col <= ncols);
row = row(idxl_roi);
col = col(idxl_roi);
v = v(idxl_roi);

ind = nan(size(xy,1),1);
ind(idxl_roi) = sub2ind([nrows, ncols], row, col);


%% rasterize

if any(idxl_roi)
    
    varargout{1} = accumarray([row, col], ...
        v, ...
        [nrows, ncols], ...
        arg.Results.fun, ...
        cast(arg.Results.fill, class(arg.Results.fun(v(1)))));
    
    varargout{2} = refmat;
    varargout{3} = ind;
    
else
    
    fprintf('\nWarning: all the points are located outside of the defined raster extent\n');
    varargout{1} = [];
    varargout{2} = refmat;
    varargout{3} = ind;
    
end


%% plot

if arg.Results.fig && any(idxl_roi)
    
    figure
    imagesc(varargout{1})
    colorbar
    axis equal tight
    
end
