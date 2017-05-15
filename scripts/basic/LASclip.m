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
%    points = 'E:\data\1143424.las';
%    outputFilepath = 'E:\data\1143424_clipped.las';
%    xv = [549042.0; 549295.0; 549295.0; 549042.0];
%    yv = [210236.5; 210236.5; 210032.5; 210032.5];
%    s = LASclip(points, [xv, yv], outputFilepath, 'verbose', true);
%
% Other m-files required: LASread.m, LASwrite.m, LASmerge.m, LASextent.m
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2016b
%
% See also: LASMERGE
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory (LASIG)
% Website: http://mparkan.github.io/Digital-Forestry-Toolbox/
% Last revision: May 15, 2017
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund (WHFF, OFEV), project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'points', @(x) ischar(x) || isstruct(x) || iscell(x)); % tf = isdir('A') % if exist(Name, 'file') == 2
addRequired(arg, 'clipper', @(x) (isnumeric(x) && (size(x,2) == 2)) || ischar(x));
addOptional(arg, 'outputFilepath', [], @(x) ischar(x) || isempty(x));
addParameter(arg, 'verbose', true, @(x) islogical(x) && (numel(x) == 1));

parse(arg, points, clipper, varargin{:});


%% check clipper is a shapefile or an array of coordinates

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
        
    case 'cell'
        
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
        
        if isdir(arg.Results.points) % input is a directory path
            
            % get intersections between clipper and tile extents
            fileIntersections = getFileIntersections(arg.Results.points, xc, yc);
            n_intersections = height(fileIntersections);
            
            if n_intersections >= 1
                
                parts = cell(n_intersections,1);
                
                for j = 1:n_intersections
                    
                    % read LAS files
                    s = LASread(fileIntersections.FILEPATH{j}, false, false);
                    
                    % clip parts
                    parts{j,1} = clip(s);
                    
                end
                
                % merge parts
                r = LASmerge(parts, [], ...
                    'verbose', arg.Results.verbose);
                
            end
            
        else % input is a single file path
            
            s = LASread(arg.Results.points, false, false);
            r = clip(s);
            
        end
        
    case 'cell' % input is a cell array of file paths
        
        % get intersections between clipper and tile extents
        fileIntersections = getFileIntersections(arg.Results.points, xc, yc);
        n_intersections = height(fileIntersections);
        
        if n_intersections >= 1
            
            parts = cell(n_intersections,1);
            
            for j = 1:n_intersections
                
                % read LAS files
                s = LASread(fileIntersections.FILEPATH{j}, false, false);
                
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
        
        outputFilepath =  [pathstr, '\', name, '_', datestr(now, 'ddmmyy_HHMM'), '.las'];
        
    else
        
        outputFilepath =  [pathstr, '\', name, '.las'];
        
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
        
        if arg.Results.verbose
            
            fprintf('computing clipper and file intersection...');
            
        end

        % get extents of all LAS files in specified input directory
        extent = LASextent(in, ...
            [], ...
            'method', 'header', ...
            'fig', false, ...
            'verbose', false);
        
        n_tiles = height(extent);
        idxl_intersections = false(n_tiles,1);
        
        % check if clipper intersects LAS file extents
        for k = 1:n_tiles
            
            [xi, yi] = polybool('intersection', ...
                xc, ...
                yc, ...
                extent.X{k}, ...
                extent.Y{k});
            
            extent.XI{k} = xi;
            extent.YI{k} = yi;
            
            if ~isempty(xi)
                
                idxl_intersections(k,1) = true;
                
            end
            
        end
        
        fileExtents = extent(idxl_intersections,:);
        
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
        r.header = s.header;
        
        for k = 1:length(pdr_skeys)
            
            r.record.(pdr_skeys{k}) = s.record.(pdr_skeys{k})(idxl_in);
            
        end
        
    end


end