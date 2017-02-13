function [M, I] = geomedian(xyz, varargin)
% GEOMEDIAN - computes the geometric median of 3D points XYZ
% 
%
% Syntax:  [M, I] = geomedian(xyz, ...)
%
% Inputs:
%    xyz - Nx3 numeric matrix, 3D point cloud x,y,z coordinates
%
%    fig (optional, default: false) - boolean value, switch to plot figures
%
% Outputs:
%    M - 1x3 numeric vector, geometric median x,y,z coordinates
%
%    I - Nx1 numeric vector, index of the points that is the geometric median
%
% Example:
%
%    [M, I] = geomedian(xyz);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2016b
%
% See also:
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory (LASIG)
% Website: http://mparkan.github.io/Digital-Forestry-Toolbox/
% Last revision: February 13, 2017
% Acknowledgments: This work was supported by the Swiss Forestry and Wood
% Research Fund, WHFF (OFEV) - project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'xyz', @(x) (size(x,2) == 3) && isnumeric(x));
addParameter(arg, 'fig', false, @(x) islogical(x) && (numel(x) == 1));

parse(arg, xyz, varargin{:});


%% compute pairwise squared distances between all points

if size(xyz,1) >= 2
    
    % find unique combinations of all point pairs
    pairs = nchoosek(1:size(xyz, 1), 2);
    
    % compute pairwise squared distances between all unique points
    dist = sum((xyz(pairs(:, 1),:) - xyz(pairs(:, 2),:)) .^ 2, 2);
    PAIRS = [pairs; pairs(:,2:-1:1)];
    DIST = [dist; dist];
    
    % compute sum of squarred distances for each point
    SDS = accumarray(PAIRS(:,1), DIST, [], @(x) sum(x));
    
    % find point that minmizes sum of squarred distances to all others
    [SDS_min, I] = min(SDS);
    M = xyz(I,:);
    
    if arg.Results.fig
        
        figure
        scatter3(xyz(:,1), ...
            xyz(:,2), ...
            xyz(:,3), ...
            6, ...
            [0,0,0], ...
            'Marker', '.')
        hold on
        scatter3(xyz(I,1), ...
            xyz(I,2), ...
            xyz(I,3), ...
            12, ...
            [1,0,0], ...
            'Marker', 'o')
        axis equal tight vis3d
        xlabel('x')
        ylabel('y')
        zlabel('z')
        title('before majority vote')
        
    end
    
else
    
    I = 1;
    M = xyz;
    
end

