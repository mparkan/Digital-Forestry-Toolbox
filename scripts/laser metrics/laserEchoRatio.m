function echo_ratio = laserEchoRatio(xyz, varargin)
%LASERECHORATIO - Computes the echo ratio of a 3D point cloud. The echo
%ratio is a measure of local opacity with values between 0 (transparent) and 1 (opaque).
% ECHO_RATIO = LASERECHORATIO(XYZ, ...) with coordinates XYZ
%
% Syntax:  echo_ratio = laserEchoRatio(xyz, ...)
%
% Inputs:
%    xyz - nx1 float vector, xyz coordinates of the point cloud 
%        
%    rasterResolution (optional, default: 0.25) - numeric value, raster cell resolution used when converting the point cloud to a 3D raster
%
%    verbose (optional, default: true) - boolean value, verbosiy switch
%
%    fig (optional, default: false) - boolean value, switch to plot figures
%
% searchRadius
%
% Outputs:
%    echo_ratio - nx2 float vector, echo ratio
%
% Example:
%    pc = LASread('..\data\measurements\vector\als\zh_6850_2475.las', false, true);
%    x = pc.record.x;
%    y = pc.record.y;
%    z = pc.record.z;
%    
%    echo_ratio = laserEchoRatio([x, y, z], 'rasterResolution', 1, 'verbose', true, 'fig', true);
%    
% Other m-files required: rasterize.m
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2016a
%
% See also:
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Laboratory
% Website: http://lasig.epfl.ch/
% Last revision: May 14, 2016
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund (WHFF, OFEV), project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'xyz', @(x) isnumeric(x) && size(x,2) == 3);
addParamValue(arg, 'rasterResolution', 0.5, @(x) isnumeric(x) && (numel(x) == 1));
addParamValue(arg, 'fig', false, @(x) islogical(x) && (numel(x) == 1));
addParamValue(arg, 'verbose', true, @(x) islogical(x) && (numel(x) == 1));

parse(arg, xyz, varargin{:});


%% compute spatial extents

x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);

x_min = min(xyz(:,1));
x_max = max(xyz(:,1));
y_min = min(xyz(:,2));
y_max = max(xyz(:,2));
z_min = min(xyz(:,3));
z_max = max(xyz(:,3));


%% rasterize (3D) canopy points

if arg.Results.verbose
    
    tic
    fprintf('rasterizing point cloud...');
    
end

step = arg.Results.rasterResolution;
xv = x_min-step:step:x_max+step;
yv = y_min-step:step:y_max+step;
zv = z_min-step:step:z_max+step;

sub_crl = rasterize([x, y, z], xv, yv, zv);
nrows = length(yv);
ncols = length(xv);
nstacks = length(zv);

ind_crl = sub2ind([nrows, ncols, nstacks], sub_crl(:,2), sub_crl(:,1), sub_crl(:,3)); % linear index for 3D raster cells
ind_cr = sub2ind([nrows, ncols], sub_crl(:,2), sub_crl(:,1)); % linear index for 2D raster cells

if arg.Results.verbose
    
    fprintf('done!\n');
    toc
    
end


%% compute 3D point density

if arg.Results.verbose
    
    tic
    fprintf('computing 3D point density...');
    
end

ind_crl_unique = unique(ind_crl);
[~, idxn_crl] = ismember(ind_crl, ind_crl_unique);
[point_density_3D, ~] = histcounts(ind_crl, [ind_crl_unique; inf]);

if arg.Results.verbose
    
    fprintf('done!\n');
    toc
    
end


%% compute 2D point density

if arg.Results.verbose
    
    tic
    fprintf('computing 2D point density...');
    
end

ind_cr_unique = unique(ind_cr);
[~, idxn_cr] = ismember(ind_cr, ind_cr_unique);
[point_density_2D, ~] = histcounts(ind_cr, [ind_cr_unique; inf]);

if arg.Results.verbose
    
    fprintf('done!\n');
    toc
    
end


%% compute echo ratio

if arg.Results.verbose
    
    tic
    fprintf('computing 3D/2D echo ratio...');
    
end

echo_ratio = point_density_3D(idxn_crl) ./ point_density_2D(idxn_cr);

if arg.Results.fig
    
    % plot echo ratio
    figure
    scatter3(x, y, z, 12, ...
        echo_ratio, ...
        'Marker', '.');
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis equal tight vis3d
    title('echo ratio')
    colorbar
    
end

if arg.Results.verbose
    
    fprintf('done!\n');
    toc
    
end
