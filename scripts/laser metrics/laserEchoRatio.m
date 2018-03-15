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
%    pc = LASread('..\data\measurements\vector\als\6850_2475.las', false, true);
%    x = pc.record.x;
%    y = pc.record.y;
%    z = pc.record.z;
%
%    echo_ratio = laserEchoRatio([pc.record.x, pc.record.y, pc.record.z], ...
%        'rasterResolution', 1, ...
%        'verbose', true, ...
%        'fig', true);
%
% Other m-files required: rasterize.m
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2017b, GNU Octave 4.2.1 (configured for "x86_64-w64-mingw32")
%
% See also:
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory (LASIG)
% Website: http://mparkan.github.io/Digital-Forestry-Toolbox/
% Last revision: March 15, 2018
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund (WHFF, OFEV), project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'xyz', @(x) isnumeric(x) && size(x,2) == 3);
addParameter(arg, 'rasterResolution', 0.5, @(x) isnumeric(x) && (numel(x) == 1));
addParameter(arg, 'method', 'standard', @(x) ismember(x, {'standard', 'cumulative'}));
addParameter(arg, 'fig', false, @(x) islogical(x) && (numel(x) == 1));
addParameter(arg, 'verbose', true, @(x) islogical(x) && (numel(x) == 1));

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
    
    fprintf('rasterizing point cloud...');
    pause(0.01);
    
end

step = arg.Results.rasterResolution;

xv = x_min:step:x_max; %x_min-scale:scale:x_max+scale;
yv = y_min:step:y_max;
zv = z_min:step:z_max;

xv = [xv, xv(end)+ceil(mod(x_max, step))*step];
yv = [yv, yv(end)+ceil(mod(y_max, step))*step];
zv = [zv, zv(end)+ceil(mod(z_max, step))*step];

[~, sub_crl] = rasterize([x, y, z], xv, yv, zv);
nrows = length(yv);
ncols = length(xv);
nstacks = length(zv);

ind_crl = sub2ind([nrows, ncols, nstacks], sub_crl(:,2), sub_crl(:,1), sub_crl(:,3)); % linear index for 3D raster cells
ind_cr = sub2ind([nrows, ncols], sub_crl(:,2), sub_crl(:,1)); % linear index for 2D raster cells

if arg.Results.verbose
    
    fprintf('done!\n');
    pause(0.01);
    
end


%% compute 3D point density

if arg.Results.verbose
    
    fprintf('computing 3D point density...');
    pause(0.01);
    
end

ind_crl_unique = unique(ind_crl);
[~, idxn_crl] = ismember(ind_crl, ind_crl_unique);

% [ind_crl_unique; inf];
% point_density_3D2 = accumarray(ind_crl, ind_crl, [], @length, 0);

% [point_density_3D, ~] = histcounts(ind_crl, [ind_crl_unique; inf]);
point_density_3D = histc(ind_crl, ind_crl_unique)';

if arg.Results.verbose
    
    fprintf('done!\n');
    pause(0.01);
    
end


%% compute cumulative 3D point density

switch arg.Results.method
    
    case 'cumulative'
        
        memstat = memory;
        memusage = nrows*ncols*nstacks*8/memstat.MaxPossibleArrayBytes;
        
        if memusage < 0.75
            
            tic
            % cumsum
            A = accumarray([sub_crl(:,2), sub_crl(:,1), sub_crl(:,3)], ones(size(x)), [nrows, ncols, nstacks], @(x) numel(x), 0); % raster
            % A = accumarray([sub_crl(:,2), sub_crl(:,1), sub_crl(:,3)], ones(size(x), 'uint16'), [nrows, ncols, nstacks], @(x) uint16(numel(x)), uint16(0)); % raster
            %A = accumarray(ind_crl, x, [nrows*ncols*nstacks, 1], @numel, 0, true); % raster
            toc
            
            lower_echo_ratio = A ./ cumsum(A,3);
            
            aa = lower_echo_ratio(ind_crl);
            
            
            if arg.Results.fig
                
                % plot cumulative echo ratio
                figure
                scatter3(x, y, z, 12, ...
                    aa, ...
                    'Marker', '.');
                xlabel('x')
                ylabel('y')
                zlabel('z')
                axis equal tight vis3d
                title('cumulative echo ratio')
                colorbar
                
            end
            
        end
        
    case 'standard'
        
        
        %% compute 2D point density
        
        if arg.Results.verbose
            
            fprintf('computing 2D point density...');
            pause(0.01)
            
        end
        
        ind_cr_unique = unique(ind_cr);
        [~, idxn_cr] = ismember(ind_cr, ind_cr_unique);
        %[point_density_2D, ~] = histcounts(ind_cr, [ind_cr_unique; inf]);
        point_density_2D = histc(ind_cr, ind_cr_unique)';
        
        if arg.Results.verbose
            
            fprintf('done!\n');
            pause(0.01);
            
        end
        
        
        %% compute echo ratio
        
        if arg.Results.verbose
            
            fprintf('computing 3D/2D echo ratio...');
            pause(0.01);
            
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
            axis equal tight
            title('standard echo ratio')
            colorbar
            
        end
        
        if arg.Results.verbose
            
            fprintf('done!\n');
            pause(0.01);
            
        end
        
end