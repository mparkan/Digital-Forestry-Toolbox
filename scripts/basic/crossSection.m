function [index, footprint, profile] = crossSection(xyz, width, p0, p1, varargin)
% CROSSSECTION - computes a 2D vertical cross-section from a 3D point cloud
%
% [INDEX, FOOTPRINT, PROFILE] = CROSSSECTION(XYZ, WIDTH, P0, P1, ...) extract a 2D vertical cross-section/profile from a 3D point cloud
%
% Syntax:  [index, footprint, profile] = crossSection(xyz, width, p0, p1, ...)
%
% Inputs:
%    xyz - A Nx3 float matrix containing the xyz coordinates of the point
%    cloud
%
%    width - A float indicating the width of the buffer used when extracting the profile 
%
%    p0 - A 1x2 float matrix containing the xy coordinates of the profile
%    start point
%
%    p1 - A 1x2 float matrix containing the xy coordinates of the profile
%    end point
%
%    verbose (optional, default: false) - boolean value, verbosiy switch
%
%    fig (optional, default: false) - boolean value, switch to plot figures
%
% Outputs:
%    index - A Nx1 integer vector containg the indices of points located in
%    the cross-section buffer
%
%    footprint - A Nx2 float matrix containing the xy coordinates of the cross-section buffer  
%
%    profile - A Nx3 float matrix containing the range, elevation and distance to cross-section plane
%    of the points in the cross-section buffer
%    
%
% Example:
%    pc = LASread('..\data\measurements\vector\als\zh_6995_2710.las', false, true);
%    x = pc.record.x;
%    y = pc.record.y;
%    z = pc.record.z;
%    
%    width = 2;
%    p0 = [699841.5, 271000.2];
%    p1 = [699784, 271494.3];
%    [~, ~, profile] = crossSection([x, y, z], width, p0, p1, 'verbose', false, 'fig', true);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2016a
%
% See also: viewerCrossSection.m
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

addRequired(arg, 'xyz', @(x) isnumeric(x) && (size(x,2) == 3));
addRequired(arg, 'width', @(x) isnumeric(x) && (numel(x) == 1));
addRequired(arg, 'p0', @(x) isnumeric(x) && (numel(x) == 2));
addRequired(arg, 'p1', @(x) isnumeric(x) && (numel(x) == 2));
addParamValue(arg, 'verbose', true, @(x) islogical(x) && (numel(x) == 1));
addParamValue(arg, 'fig', true, @(x) islogical(x) && (numel(x) == 1));

parse(arg, xyz, width, p0, p1, varargin{:});


%% reformat input

x = xyz(:,1);
y = xyz(:,2);
w = width / 2;


%% create profile buffer

% compute vector between p0 and p1
v = p1 - p0;

% compute normals
n1 = [-v(2), v(1)];
n2 = [v(2), -v(1)];
n1 = n1 / norm(n1);
n2 = n2 / norm(n2);

% compute buffer corner points
footprint = [p0 + w * n1;
    p0 + w * n2;
    p1 + w * n2;
    p1 + w * n1];

idxl_buffer = inpolygon(x, y, footprint(:,1), footprint(:,2));
pq = xyz(idxl_buffer,:);
pb = xyz(~idxl_buffer,:);
index = find(idxl_buffer);


%% plot footprint extent

if arg.Results.fig
    
    % plot cross section footprint
    figure
    plot([footprint(:,1); footprint(1,1)], [footprint(:,2); footprint(1,2)], 'ro-')
    hold on
    plot(p0(:,1), p0(:,2), 'bx')
    hold on
    plot(p1(:,1), p1(:,2), 'bx')
    axis equal
    
    % plot 3D points with highlighted cross section 
    figure
    scatter3(pb(:,1), pb(:,2), pb(:,3), 8, pb(:,3), 'Marker', '.');
    hold on
    scatter3(pq(:,1), pq(:,2), pq(:,3), 8, [1, 0.65, 0], 'Marker', '.');
    colormap('gray')
    axis equal tight vis3d
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
end


%% project points on 2D vector

% direction vector of the line
line = [p0(:,1), p0(:,2), p1(:,1)-p0(:,1), p1(:,2)-p0(:,2)]; 
 
vx = line(:,3);
vy = line(:,4);

% difference of point with line origin
dx = pq(:,1) - line(:,1);
dy = pq(:,2) - line(:,2);

% position of projection on line, using dot product
tp = (dx .* vx + dy .* vy ) ./ (vx .* vx + vy .* vy);

% convert position on line to 2D cartesian coordinates
v_proj = [line(:,1) + tp .* vx, line(:,2) + tp .* vy];

% horizontal distance between 3D point and projection plane
d = sqrt(sum((v_proj - pq(:,1:2)).^2, 2));


%% compute linear coordinates

%rho = [v_proj(:,1)-p0(:,1) v_proj(:,2)-p0(:,2)];
rho = sqrt((v_proj(:,1) - p0(:,1)).^2 + (v_proj(:,2) - p0(:,2)).^2);

profile = [rho, pq(:,3), d];


%% plot profile

if arg.Results.fig
    
    % plot cross section
    figure
    plot(rho, pq(:,3), 'k.')
    ylabel('z')
    xlabel('rho')
    axis equal tight
    
end
