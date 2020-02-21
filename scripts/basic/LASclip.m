function s = LASclip(points, clipper, varargin)
%LASCLIP - Clips a LAS file with a polygon
% LASCLIP(POINTS, CLIPPER, ...) clips the point cloud stored in
% POINTS with CLIPPER. Optionally writes the clipped point cloud to
% OUTPUTFILEPATH.
%
% Syntax:  LASclip(points, clipper, ...)
%
% Inputs:
%    points - Either a LAS structure, a path to a single LAS file, a cell array 
%    containing paths to multiple LAS files, a path to a directory
%    containing LAS files (Matlab only)
%
%    clipper - Either a Nx2 matrix containing the x,y coordinates of the
%    clipper polygon or a string indicating the path to an ESRI shapefile
%    containing the clipper polygon
%
%    outputFilepath (optional, default: []) - The path to the clipped output LAS file
%
%    UUID (optional, default: auto) - 32 character Univerally Unique
%    Identifier (e.g. '401db5f076cc277abc4aa217f593c48b').
%    If not specified a new UUID is automatically assigned.
%
%    verbose (optional, default: false) - boolean value, verbosiy switch
%
% Example:
%    points = 'E:\data\1143424.las';
%    outputFilepath = 'E:\data\1143424_clipped.las';
%    xv = [549042.0; 549295.0; 549295.0; 549042.0];
%    yv = [210236.5; 210236.5; 210032.5; 210032.5];
%    s = LASclip(points, [xv, yv], outputFilepath, 'verbose', true);
%
% Other m-files required: LASread.m, LASwrite.m, LASmerge.m, LASextent.m
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2019b, GNU Octave 5.2.0 (configured for "x86_64-w64-mingw32")
%
% See also: LASMERGE
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory (LASIG)
% Website: http://mparkan.github.io/Digital-Forestry-Toolbox/
% Last revision: February 20, 2020
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund (WHFF, OFEV), project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'points', @(x) ischar(x) || isstruct(x) || iscell(x)); % tf = isdir('A') % if exist(Name, 'file') == 2
addRequired(arg, 'clipper', @(x) (isnumeric(x) && (size(x,2) == 2)) || ischar(x));
addOptional(arg, 'outputFilepath', [], @(x) ischar(x) || isempty(x));
addParameter(arg, 'UUID', lower(strcat(dec2hex(randi(16,32,1)-1)')), @(x) ischar(x) && length(x) == 32);
addParameter(arg, 'verbose', true, @(x) islogical(x) && (numel(x) == 1));

parse(arg, points, clipper, varargin{:});


%% load additional packages in Octave

OCTAVE_FLAG = (exist('OCTAVE_VERSION', 'builtin') ~= 0); % determine if system is Matlab or GNU Octave

if OCTAVE_FLAG
    
    pkg load statistics

end

%% check if clipper is a shapefile or an array of coordinates

if arg.Results.verbose
    
     fprintf('reading clipper...');
    
end

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
        
end

if arg.Results.verbose
    
    fprintf('done!\n');
    
end


%% check clip polygon coordinates

if arg.Results.verbose
    
     fprintf('checking clipper coordinates...');
    
end

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

if arg.Results.verbose
    
    fprintf('done!\n');
    
end


%% check input point format

switch class(arg.Results.points)
    
    case 'char'
        
        if isfolder(arg.Results.points) % check if input is a directory path
            
            if OCTAVE_FLAG
                
                error('Clipping from a directory is not yet supported in Octave')
                
            end
            
            % get intersections between clipper and tile extents
            fileIntersections = getFileIntersections(arg.Results.points, xc, yc);
            n_intersections = length(fileIntersections);
            
            if n_intersections >= 1
                
                parts = cell(n_intersections,1);
                
                for j = 1:n_intersections
                    
                    % read LAS files
                    s = LASread(fileIntersections(j).FILEPATH, false, false);
                    
                    % clip parts
                    parts{j,1} = clip(s);
                    
                end
                
                % merge parts
                r = LASmerge(parts, [], ...
                    'verbose', arg.Results.verbose);
                
            else
                
                error('The clipper does not intersect the input file(s)')
                
            end
            
        else % input is a single file path
            
            s = LASread(arg.Results.points, false, false);
            r = clip(s);
            
        end
        
    case 'cell' % input is a cell array of file paths
        
        % get intersections between clipper and tile extents
        %fileIntersections = getFileIntersections(arg.Results.points, xc, yc);
        %n_intersections = height(fileIntersections);
        
        n_intersections = length(arg.Results.points);
        
        if n_intersections >= 1
            
            parts = cell(n_intersections,1);
            
            for j = 1:n_intersections
                
                % read LAS files
                %s = LASread(fileIntersections.FILEPATH{j}, false, false);
                s = LASread(arg.Results.points{j}, false, false);
                
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
    [pathstr, name, ~] = fileparts(outputFilepath);
    
    % adjust output file name
    if exist(outputFilepath, 'file')
        
        if ~isempty(pathstr)
            
            outputFilepath =  [pathstr, filesep, name, '_', datestr(now, 'ddmmyy_HHMM'), '.las'];
            
        else
            
            outputFilepath =  [name, '_', datestr(now, 'ddmmyy_HHMM'), '.las'];
            
        end
        
        
    else
        
        if ~isempty(pathstr)
            
            outputFilepath =  [pathstr, filesep, name, '.las'];
            
        else
            
            outputFilepath =  [name, '.las'];
            
        end
        
    end
    
    las_version = s.header.version_major * 10 + s.header.version_minor;
        
    LASwrite(r, outputFilepath, ...
        'version', las_version, ...
        'systemID', 'EXTRACTION', ...
        'guid', arg.Results.UUID, ...
        'recordFormat', s.header.point_data_format_id, ...
        'verbose', arg.Results.verbose);

end


%% get extent of tiles from LAS files

    function fileExtents = getFileIntersections(in, xc, yc)
        
        if arg.Results.verbose
            
            fprintf('computing clipper and file intersection...');
            
        end

        % get extents of all LAS files in specified input directory
        extent = LASextent(in, ...
            [], ...
            'method', 'header', ...
            'fig', false, ...
            'verbose', false);
        
        n_tiles = length(extent);
        idxl_intersections = false(n_tiles,1);
        
        % check if clipper intersects LAS file extents
        for k = 1:n_tiles
            
            % not supported in octave
            [xi, yi] = polybool('intersection', ...
                xc, ...
                yc, ...
                extent(k).X, ...
                extent(k).Y);
            
            extent(k).XI = xi;
            extent(k).YI = yi;
            
            if ~isempty(xi)
                
                idxl_intersections(k,1) = true;
                
            end
            
        end
        
        fileExtents = extent(idxl_intersections);
        
        if arg.Results.verbose
            
            fprintf('done!\n');
            
        end
        
    end


%% clip point cloud

    function r = clip(s)
        
        % find points in polygon
        idxl_in = inpolygon(s.record.x, s.record.y, xc, yc);

        pdr_skeys = fieldnames(s.record);
        
        % initialize point cloud structure
        r = s;
        
        for k = 1:length(pdr_skeys)
            
            r.record.(pdr_skeys{k}) = s.record.(pdr_skeys{k})(idxl_in);
            
        end
        
        
    end

end