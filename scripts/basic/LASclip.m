function s = LASclip(points, clipper, varargin)
%LASCLIP - Clips a LAS file with a polygon
% LASCLIP(POINTS, CLIPPER, ...) clips the point cloud stored in
% POINTS with CLIPPER. Optionally writes the clipped point cloud to OUTPUTFILEPATH
%
% Syntax:  LASclip(points, clipper, ...)
%
% Inputs:
%    points - Either a LAS structure or the path to an input LAS file
%
%    clipper - Either a Nx2 matrix containing the x,y coordinates of the
%    clipper polygon or a string indicating the path to an ESRI shapefile
%    containing the clipper polygon
%
%    outputFilepath (optional, default: []) - The path to the clipped output LAS file
%
%    verbose (optional, default: false) - boolean value, verbosiy switch
%
% Example:
%    points = '..\data\measurements\vector\als\6995_2710.las';
%    outputFilepath = '..\data\measurements\vector\als\zh_6995_2710_subset.las';
%    xv = [699673; 699729; 699500; 699503];
%    yv = [271341; 271036; 271042; 271346];
%    s = LASclip(points, [xv, yv], outputFilepath, 'verbose', true);
%
% Other m-files required: LASread.m, LASwrite.m, LASmerge.m, LASextent.m
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2016a
%
% See also: LASMERGE
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory
% Website: http://lasig.epfl.ch/
% Last revision: December 12, 2016
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund (WHFF, OFEV), project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% setup constants

OCTAVE_FLAG = (exist('OCTAVE_VERSION', 'builtin') ~= 0); % determine if system is Matlab or GNU Octave


%% check argument validity

arg = inputParser;

addRequired(arg, 'points', @(x) ischar(x) || isstruct(x) || iscell(x)); % tf = isdir('A') % if exist(Name, 'file') == 2
addRequired(arg, 'clipper', @(x) (isnumeric(x) && (size(x,2) == 2)) || ischar(x));
addOptional(arg, 'outputFilepath', [], @(x) ischar(x) || isempty(x));
addParamValue(arg, 'verbose', true, @(x) islogical(x) && (numel(x) == 1));

parse(arg, points, clipper, varargin{:});


%% check clipper is a shapefile or an array of coordinates

switch class(arg.Results.clipper)
    
    case 'char' % shapefile
        
        shape = shaperead(arg.Results.clipper);
        
        xc = shape.X;
        yc = shape.Y;
        idxl_nan = (isnan(xc) & isnan(yc)); 
        xc(idxl_nan) = [];
        yc(idxl_nan) = [];
        
    case 'double' % array of polygon coordinates
        
        xc = arg.Results.clipper(:,1);
        yc = arg.Results.clipper(:,2);
        
    case 'cell'
        
end


%% check clip polygon coordinates

if size(xc, 2) ~= 1
    
    xc = xc';
    
end

if size(yc, 2) ~= 1
    
    yc = yc';
    
end

% check polygon coordinates closure
if xc(1) ~= xc(end) || yc(1) ~= yc(end)
    
    xc = [xc; xc(1)];
    yc = [yc; yc(1)];
    
end


%% check input point format

switch class(arg.Results.points)
    
    case 'char'
        
        if isdir(arg.Results.points) % input is a directory path
            
            % get intersections between clipper and tile extents
            fileIntersections = getFileIntersections(arg.Results.points, xc, yc);
            n_intersections = length(fileIntersections);
            
            if n_intersections >= 1
                
                parts = cell(n_intersections,1);
                for j = 1:n_intersections
                    
                    % read LAS files
                    s = LASread(fileIntersections(j).filepath, false, false);
                    
                    % clip parts
                    parts{j,1} = clip(s);
                    
                end
                
                % merge parts
                r = LASmerge(parts, [], 'verbose', arg.Results.verbose);
                
            end
            
        else % input is a single file path
            
            s = LASread(arg.Results.points, false, false);
            r = clip(s);
            
        end
        
    case 'cell' % input is a cell array of file paths
        
        % get intersections between clipper and tile extents
        fileIntersections = getFileIntersections(arg.Results.points, xc, yc);
        n_intersections = length(fileIntersections);
        
        if n_intersections >= 1
            
            parts = cell(n_intersections,1);
            for j = 1:n_intersections
                
                % read LAS files
                s = LASread(fileIntersections(j).filepath, false, false);
                
                % clip parts
                parts{j,1} = clip(s);
                
            end
            
            % merge parts
            r = LASmerge(parts, [], 'verbose', false);
            
        end
        
    case 'struct' % input is a LAS structure (as returned by LASread)
        
        s = arg.Results.points;
        r = clip(s);
        
end


%% write clipped points to LAS file

if ~isempty(arg.Results.outputFilepath) && length(xc) >= 1
    
    outputFilepath = arg.Results.outputFilepath;
    [pathstr, name, ext] = fileparts(outputFilepath);
    
    % adjust output file name
    if exist(outputFilepath, 'file')
        
        outputFilepath =  [pathstr, filesep, name, '_v2.las'];
        
    else
        
        outputFilepath =  [pathstr, filesep, name, '.las'];
        
    end
    
    las_version = s.header.version_major * 10 + s.header.version_minor;
        
    LASwrite(r, outputFilepath, ...
        'version', las_version, ...
        'systemID', 'EXTRACTION', ...
        'recordFormat', s.header.point_data_format_id, ...
        'verbose', arg.Results.verbose);

end


%% get extent of tiles from LAS files
    function fileExtents = getFileIntersections(in, xc, yc)
        
        % get extents of all LAS files in specified input directory
        extent = LASextent(in, [], 'fig', false, 'verbose', false); % , 'fig', true, 'verbose', false
        
        n_tiles = length(extent);
        idxl_intersections = false(n_tiles,1);
        
        % check if clipper intersects LAS file extents
        for k = 1:n_tiles
            
            [xi, yi] = polybool('intersection', xc, yc, extent(k).xpoly, extent(k).ypoly);
            extent(k).x_intersection = xi;
            extent(k).y_intersection = yi;
            
            if ~isempty(extent(k).x_intersection)
                
                idxl_intersections(k,1) = true;
                
            end
            
        end
        
        fileExtents = extent(idxl_intersections);
        
    end


%% clip point cloud
    function r = clip(s)
        
        % find points in polygon
        idxl_in = inpolygon(s.record.x, s.record.y, xc, yc);

        pdr_skeys = fieldnames(s.record);
        
        % initialize point cloud structure
        r.header = s.header;
        
        for k = 1:length(pdr_skeys)
            
            r.record.(pdr_skeys{k}) = s.record.(pdr_skeys{k})(idxl_in);
            
        end
        
    end


end