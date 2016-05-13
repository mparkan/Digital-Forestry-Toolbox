function extent = LASextent(points, varargin)
%LASEXTENT - extracts the extent of LAS files from the header information
% EXTENT = LASEXTENT(POINTS, ...) extracts the extents (bounding boxes)
% of all the LAS files in POINTS and optionally export them to a single
% ESRI shapefile specified by OUTPUTFILEPATH
%
% Syntax:  extent = LASextent(folder, outputFilepath, ...)
%
% Inputs:
%    points - The path to the folder containing the LAS files
%
%    outputFilepath (optional) - The path to the output shapefile
%
%    verbose (optional, default: false) - boolean value, verbosiy switch
%
%    fig (optional, default: false) - boolean value, switch to plot figures
%
% Outputs:
%    extent - A structure containing the extent coordinates of each LAS file
%
% Example:
%    points = '..\data\measurements\vector\als\';
%    outputFilepath = '..\data\measurements\vector\extents\als_extents.shp';
%    extent = LASextent(points, outputFilepath, 'fig', true, 'verbose', false);
%
% Other m-files required: LASread.m, LASwrite.m
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2016a
%
% See also: LASCLIP, LASMERGE
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

addRequired(arg, 'points', @(x) iscell(x) || (ischar(x) && isdir(x)));
addOptional(arg, 'outputFilepath', [], @(x) ischar(x) || isempty(x));
addParamValue(arg, 'fig', false, @(x) islogical(x) && (numel(x) == 1));
addParamValue(arg, 'verbose', false, @(x) islogical(x) && (numel(x) == 1));

parse(arg, points, varargin{:});


%% list LAS input files

switch class(points)
    
    case 'cell'
        
        filelist = struct;
        
        for j = 1:length(points)
            
            [pathstr, name, ext] = fileparts(points{j});
            filelist(j).path = pathstr;
            filelist(j).name = name;
            filelist(j).ext = ext;
            
        end
        
    case 'char'
        
        if ~strcmpi(points(end), '\')
            
           points = [points, '\'];
            
        end
            
        filelist = dir([points, '*.las']);
        
        for j = 1:length(filelist)
            
            [pathstr, name, ext] = fileparts([points, filelist(j).name]);
            
            filelist(j).path = pathstr;
            filelist(j).name = name;
            filelist(j).ext = ext;
            
        end
        
end


%% check file extension

% idx_ext = strcmpi('.las', {filelist.ext})'; % index of LAS format files
% idx_file = ~[filelist.isdir]'; % index of non-directories
% idx_size = [filelist.bytes]' > 0; % index of non-empty files
% 
%idx_read = find(idx_ext & idx_file & idx_size);
% 
%n_tiles = length(idx_read);


%% read LAS header

n_tiles = length(filelist);
extent = struct;
reverse_str = '';

for j = 1:n_tiles
    
    extent(j).filepath = [filelist(j).path, '\', filelist(j).name, filelist(j).ext]; % [points, filelist(idx_read(j)).name];
    pc = LASread(extent(j).filepath, true, false);
    
    extent(j).name = filelist(j).name; % filelist(idx_read(j)).basename;
    extent(j).xmin = pc.header.min_x;
    extent(j).xmax = pc.header.max_x;
    extent(j).ymin = pc.header.min_y;
    extent(j).ymax = pc.header.max_y;
    extent(j).xpoly = [extent(j).xmin, extent(j).xmax, extent(j).xmax, extent(j).xmin, extent(j).xmin];
    extent(j).ypoly = [extent(j).ymax, extent(j).ymax, extent(j).ymin, extent(j).ymin, extent(j).ymax];
    
    if arg.Results.verbose
        
        msg = sprintf('processed %d/%d', j, n_tiles);
        fprintf([reverse_str, msg]);
        reverse_str = repmat(sprintf('\b'), 1, length(msg));
        
    end
    
end

if arg.Results.verbose
    
    fprintf('\n');
    
end


%% plot extents

if arg.Results.fig
    
    figure
    
    for j = 1:n_tiles
        
        plot([extent(j).xmin, extent(j).xmax, extent(j).xmax, extent(j).xmin, extent(j).xmin], ...
            [extent(j).ymax, extent(j).ymax, extent(j).ymin, extent(j).ymin, extent(j).ymax], 'r-')
        hold on
        
    end
    
    axis equal tight
    title('Extents')
    xlabel('x')
    ylabel('y')
    
end


%% export to shapefile

if ~isempty(arg.Results.outputFilepath)
    
    if arg.Results.verbose
        
        fprintf('writing extent to "%s"...', arg.Results.outputFilepath);
        
    end
    
    geom(1:n_tiles,:) = {'Polygon'};
    S = struct('Geometry', geom,...
        'ID', num2cell(1:n_tiles)',...
        'NAME', {extent.name}',...
        'FILE', {extent.filepath}',...
        'X', {extent.xpoly}',...
        'Y', {extent.ypoly}');
    
    shapewrite(S, arg.Results.outputFilepath);
    
    if arg.Results.verbose
        
        fprintf('done!\n');
        
    end
    
end
