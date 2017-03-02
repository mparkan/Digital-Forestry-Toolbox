function [Xs, idxn_cell] = subsample(X, varargin)
% SUBSAMPLE - subsample a 2D or 3D point cloud 
% [XS, IDXN_CELL] = SUBSAMPLE(X) subsamples the 2D or 3D point cloud (X) by voxelizing the points
% and applying the aggregate function specified by METHOD to each voxel (e.g. mean, geometric median, grid, random).
% 
% Syntax:  [Xs, idxn_cell] = subsample(X, ...)
% 
% Inputs:
%    X - Nx3 or Nx2 numeric matrix, 3D or 2D point coordinates
% 
%    method (optional, default: 'mean') - any of the following 'geomedian', 'mean', 'grid', 'rand'
% 
%    resolution (optional, default: 0.5) - numeric value, voxel resolution
%
%    verbose (optional, default: true) - boolean value, verbosiy switch
%
%    fig (optional, default: false) - boolean value, switch to plot figures
%
% Outputs:
%    Xs - Nx3 or Nx2 numeric matrix, 3D or 2D point subsampled point coordinates
% 
%    idxn_cell - Nx1 integer vector, index indicating to which subsample voxel each
%    point was assigned to
% 
% Example:
%
% [xyz_s, idxn_cell] = subsample(xyz, ...
%     'method', 'mean', ...
%     'resolution', 0.5, ...
%     'fig', true, ...
%     'verbose', false);
% 
% Other m-files required: geomedian.m
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2016b
% 
% See also:
% 
% This code is part of the Matlab Digital Forestry Toolbox
% 
% Author: Matthew Parkan, EPFL - GIS Laboratory (LASIG)
% Website: http://mparkan.github.io/Digital-Forestry-Toolbox/
% Last revision: March 02, 2017
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund (WHFF, OFEV), project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'X', @(x) (size(x,2) >= 2) && (size(x,2) <= 3) && isnumeric(x));
addParameter(arg, 'method', 'mean', @(x) ismember(x, {'geomedian', 'mean', 'grid', 'rand'}));
addParameter(arg, 'resolution', 0.5, @(x) isnumeric(x) && (x >= 0) && (numel(x) == 1));
addParameter(arg, 'fig', true, @(x) islogical(x) && (numel(x) == 1));
addParameter(arg, 'verbose', true, @(x) islogical(x) && (numel(x) == 1));

parse(arg, X, varargin{:});


%% subsample point cloud (decimation)

if arg.Results.verbose
    
    fprintf('subsampling point cloud...');
    
end

switch arg.Results.method % 'geomedian', 'mean', 'grid', 'none'
    
    case 'geomedian' % geometric median
        
        [~, ~, idxn_cell] = unique(round(X / arg.Results.resolution) * arg.Results.resolution, 'rows');
        [Xs, ~] = splitapply(@geomedian, X, idxn_cell);
        
    case 'mean' % mean
        
        [~, ~, idxn_cell] = unique(round(X / arg.Results.resolution) * arg.Results.resolution, 'rows');
        Xs = splitapply(@(x) mean(x,1), X, idxn_cell);
        
    case 'grid' % grid
        
        [Xs, ~, idxn_cell] = unique(round(X / arg.Results.resolution) * arg.Results.resolution, 'rows');
        
    case 'rand' % random
        
        [~, ~, idxn_cell] = unique(round(X / arg.Results.resolution) * arg.Results.resolution, 'rows');
        Xs = splitapply(@(x) x(randperm(size(x,1),1),:), X, idxn_cell);
           
end

if arg.Results.verbose
    
    fprintf('done!\n');
    
end


%% plot figure

if arg.Results.fig
    
    switch size(X,2)
        
        case 2
        
            figure
            scatter(Xs(:,1), ...
                Xs(:,2), ...
                3, ...
                [0,0,0], ...
                'Marker', '.')
            axis equal tight
            xlabel('x')
            ylabel('y')
            title(sprintf('Subsampled points (resolution: %.2f)', arg.Results.resolution))
            
        case 3
            
            figure
            scatter3(Xs(:,1), ...
                Xs(:,2), ...
                Xs(:,3), ...
                3, ...
                [0,0,0], ...
                'Marker', '.')
            axis equal tight
            xlabel('x')
            ylabel('y')
            zlabel('z')
            title(sprintf('Subsampled points (resolution: %.2f)', arg.Results.resolution))
            
    end
    
end