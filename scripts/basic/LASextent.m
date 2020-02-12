function extent = LASextent(points, varargin)
%LASEXTENT - extracts the extent of LAS files from the header information
% EXTENT = LASEXTENT(POINTS, ...) extracts the extents (bounding boxes)
% of all the LAS files in POINTS and optionally (Matlab only) export them to a single
% ESRI shapefile specified by OUTPUTFILEPATH
%
% Syntax:  extent = LASextent(folder, outputFilepath, ...)
%
% Inputs:
%    points - The path to the folder containing the LAS files
%
%    outputFilepath (optional, default: []) - The path to the output shapefile
%
%    method (optional, default: 'bbox') - The type of extent ('header': use values in header, 
%    'bbox': compute rectangular bounding box, 'concavehull', computes the concave hull (Matlab only), 
%    'convexhull', computes the convex hull 
%
%    verbose (optional, default: false) - boolean value, verbosiy switch
%
%    fig (optional, default: false) - boolean value, switch to plot figures
%
% Outputs:
%    extent - A structure containing the extent coordinates of each LAS file
%
% Example:
%    points = '..\data\measurements\vector\als\zh_6995_2710_coniferous.las';
%    outputFilepath = '..\data\measurements\vector\extents\zh_6995_2710_coniferous.shp';
%    extent = LASextent(points, ...
%         outputFilepath, ...
%         'method', 'convexhull', ...
%         'fig', true, ...
%         'verbose', false);
%
% Other m-files required: LASread.m, LASwrite.m, subsample.m, gps2utc.m
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2017b, GNU Octave 4.2.1 (configured for "x86_64-w64-mingw32")
% 
% See also: LASCLIP, LASMERGE
% 
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory (LASIG)
% Website: http://mparkan.github.io/Digital-Forestry-Toolbox/
% Last revision: March 28, 2018
% Acknowledgments: This work was supported by the Swiss Forestry and Wood
% Research Fund, WHFF (OFEV) - project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'points', @(x) iscell(x) || ischar(x));
addOptional(arg, 'outputFilepath', [], @(x) ischar(x) || isempty(x));
addParameter(arg, 'method', 'bbox', @(x) ismember(x, {'header', 'bbox', 'concavehull', 'convexhull'}));
addParameter(arg, 'fig', false, @(x) islogical(x) && (numel(x) == 1));
addParameter(arg, 'verbose', false, @(x) islogical(x) && (numel(x) == 1));

parse(arg, points, varargin{:});

%% load additional packages in Octave

OCTAVE_FLAG = (exist('OCTAVE_VERSION', 'builtin') ~= 0); % determine if system is Matlab or GNU Octave

if OCTAVE_FLAG
    
    pkg load statistics
    pkg load image

end

%% list LAS input files

if arg.Results.verbose
   
    fprintf('creating file list...');
    pause(0.01);
    
end
    
switch class(points)
    
    case 'cell' % input is a collection of LAS file
        
        filelist = struct;
        
        for j = 1:length(points)
            
            [pathstr, name, ext] = fileparts(points{j});
            filelist(j).path = pathstr;
            filelist(j).name = name;
            filelist(j).ext = ext;
            
        end
        
    case 'char' % input is a single LAS file
        
        if isdir(points)
            
            if ~strcmpi(points(end), filesep)
                
                points = [points, filesep];
                
            end
        
            filelist = dir([points, '*.las']);
            for j = 1:length(filelist)
                
                [pathstr, name, ext] = fileparts([points, filelist(j).name]);
                
                filelist(j).path = pathstr;
                filelist(j).name = name;
                filelist(j).ext = ext;
                
            end
            
        else
            
            [pathstr, name, ext] = fileparts(points);
            filelist(1).path = pathstr;
            filelist(1).name = name;
            filelist(1).ext = ext;
            
        end
        
end

if arg.Results.verbose
    
    fprintf('done!\n');
    pause(0.01);
    
end
    

%% compute metadata

n_files = length(filelist);
extent = struct();
warning('off');

for j = 1:n_files
    
    if ~isempty(filelist(j).path)
        
        extent(j,1).FILEPATH = [filelist(j).path, filesep, filelist(j).name, filelist(j).ext];
        
    else
        
        extent(j,1).FILEPATH = [filelist(j).name, filelist(j).ext];
        
    end
    
    extent(j,1).FILENAME = filelist(j).name;
    
    % read LAS file
    switch arg.Results.method
        
        case 'header'
            
            pc = LASread(extent(j,1).FILEPATH, true, false);
            
            % read extrema
            extent(j,1).XMIN = pc.header.min_x;
            extent(j,1).XMAX = pc.header.max_x;
            extent(j,1).YMIN = pc.header.min_y;
            extent(j,1).YMAX = pc.header.max_y;
            extent(j,1).ZMIN = pc.header.min_z;
            extent(j,1).ZMAX = pc.header.max_z;
            
        otherwise
            
            if arg.Results.verbose
                
                fprintf('reading file %u/%u...', j, n_files);
                pause(0.01);
                
            end
    
            pc = LASread(extent(j,1).FILEPATH, false, false);
            
            if arg.Results.verbose
                
                fprintf('done!\n');
                pause(0.01);
                
            end
            
            
            if arg.Results.verbose
                
                fprintf('extracting metadata of file %u/%u...', j, n_files);
                pause(0.01);
                
            end
            
            % compute extrema
            extent(j,1).XMIN = min(pc.record.x);
            extent(j,1).XMAX = max(pc.record.x);
            extent(j,1).YMIN = min(pc.record.y);
            extent(j,1).YMAX = max(pc.record.y);
            extent(j,1).ZMIN = min(pc.record.z);
            extent(j,1).ZMAX = max(pc.record.z);
            
            % compute point density (per unit area)
            [~, idxn_unit] = subsample([pc.record.x, pc.record.y], ...
                'method', 'grid', ...
                'resolution', 1, ...
                'fig', false, ...
                'verbose', false);
            
            extent(j,1).PDENSITY = median(accumarray(idxn_unit, idxn_unit, [], @numel));
            
            % compute date metadata (if GPS time is in satellite format)
            if isfield(pc.header, 'global_encoding_gps_time_type')
                
                if pc.header.global_encoding_gps_time_type == 1
                    
                    utc_time = gps2utc(pc.record.gps_time + 10^9, ...
                        'verbose', false);
                    
                    [Y,M,D,~,~,~] = datevec(utc_time);
                    
                    unique_date = sort(unique(datenum(Y,M,D)));
                    start_date = min(unique_date);
                    end_date = max(unique_date);
                    
                    % number of surveys
                    extent(j,1).NSURVEYS = length(unique_date);
                    [extent(j,1).YSTART, extent(j,1).MSTART, extent(j,1).DSTART,~,~,~] = datevec(start_date);
                    [extent(j,1).YEND, extent(j,1).MEND, extent(j,1).DEND,~,~,~] = datevec(end_date);
                    extent(j,1).DATES = strjoin(cellstr(datestr(unique_date, 'yyyy/mm/dd')), ',');
                    
                end
                
            end
            
            % subsample point cloud
            xrange = extent(j,1).XMAX - extent(j,1).XMIN;
            yrange = extent(j,1).YMAX - extent(j,1).YMIN;
            max_range = max([xrange, yrange]);
            resolution = min(3, max_range / 100);
            
            [xy_s, idxn_cell] = subsample([pc.record.x, pc.record.y], ...
                'method', 'grid', ...
                'resolution', resolution, ...
                'fig', false, ...
                'verbose', false);
            
            if arg.Results.verbose
                
                fprintf('done!\n');
                pause(0.01);
                
            end

    end
    
    if arg.Results.verbose
        
        fprintf('computing spatial extent of file %u/%u...', j, n_files);
        pause(0.01);
        
    end
             
    % compute spatial extent
    switch arg.Results.method
        
        case 'header'
            
            extent(j,1).X = [extent(j,1).XMIN, extent(j,1).XMAX, extent(j,1).XMAX, extent(j,1).XMIN, extent(j,1).XMIN];
            extent(j,1).Y = [extent(j,1).YMAX, extent(j,1).YMAX, extent(j,1).YMIN, extent(j,1).YMIN, extent(j,1).YMAX];
        
        case 'bbox'
        
            extent(j,1).X = [extent(j,1).XMIN, extent(j,1).XMAX, extent(j,1).XMAX, extent(j,1).XMIN, extent(j,1).XMIN];
            extent(j,1).Y = [extent(j,1).YMAX, extent(j,1).YMAX, extent(j,1).YMIN, extent(j,1).YMIN, extent(j,1).YMAX];
        
        case 'concavehull'
            
            if OCTAVE_FLAG
                
              error('The concavehull option is not yet supported in Octave')
                
            end  
    
            % compute single region low resolution alpha hull
            shp = alphaShape(xy_s(:,1), xy_s(:,2), 10);
            ac = criticalAlpha(shp, 'one-region');
            shp_lr = alphaShape(xy_s(:,1), xy_s(:,2), ac);
            
            % find border points (near the alpha boundary)
            [~, d] = nearestNeighbor(shp_lr, xy_s);
            idxl_border_s = (d < resolution);
            idxl_inner_s = (d >= resolution);

            idxn_border = ismember(idxn_cell, find(idxl_border_s));
            
            % compute single region high resolution alpha hull
            x_hr = [xy_s(idxl_inner_s,1); pc.record.x(idxn_border)];
            y_hr = [xy_s(idxl_inner_s,2); pc.record.y(idxn_border)];
            
            shp_hr = alphaShape(x_hr, ...
                y_hr, ...
                ac, ...
                'HoleThreshold', shp.area);

            idxn_boundary = shp_hr.boundaryFacets;

            extent(j,1).X = shp_hr.Points(idxn_boundary,1);
            extent(j,1).Y = shp_hr.Points(idxn_boundary,2);

        case 'convexhull'
            
            idxn_boundary = convhull(pc.record.x, pc.record.y);

            extent(j,1).X = pc.record.x(idxn_boundary);
            extent(j,1).Y = pc.record.y(idxn_boundary);
            
    end
    
    if arg.Results.verbose
        
        fprintf('done!\n');
        pause(0.01);
        
    end
    
end

pause(0.1);


%% export to shapefile

if ~isempty(arg.Results.outputFilepath)
    
    if arg.Results.verbose
        
        fprintf('\nwriting result to "%s"...', arg.Results.outputFilepath);
        pause(0.01);
        
    end
    
        %[extent.Geometry] = deal('Polygon');
        
        for j = 1:length(extent)
            
           extent(j,1).Geometry = 'Polygon';
           extent(j,1).BoundingBox = [extent(j,1).XMIN, extent(j,1).YMIN; extent(j,1).XMAX, extent(j,1).YMAX];
            
        end
        
        figure
        plot(extent(j,1).X, extent(j,1).Y)
        axis equal tight
        
        j = 1;
        extent2 = struct;
        extent2(j,1).Geometry = 'Polygon';
        extent2(j,1).X = [597620 597630 597630 597620 597620];
        extent2(j,1).Y = [170230 170230 170220 170220 170230];
        extent2(j,1).BoundingBox = [min(extent2(j).X), min(extent2(j).Y); max(extent2(j).X), max(extent2(j).Y)];
        extent2(j,1).Custom1 = 'cdf';        
        
        
        extent2(j,1).BoundingBox = [extent(j,1).XMIN, extent(j,1).YMIN; extent(j,1).XMAX, extent(j,1).YMAX];
        extent2(j,1).X = extent(j,1).X';
        extent2(j,1).Y = extent(j,1).Y';
        extent2(j,1).Custom1 = 'cdf';

        shapewrite(extent, arg.Results.outputFilepath);
        
        
        s = struct;
        j = 1;
        s(j,1).Geometry = 'Polygon';
        s(j,1).X = [590610 590620 590620 590610 590610];
        s(j,1).Y = [178220 178220 178210 178210 178220];
        
        
        
        s(j,1).BoundingBox = [min(s(j,1).X), min(s(j,1).Y); max(s(j,1).X), max(s(j,1).Y)];
        s(j,1).Custom1 = 'abc';
        s(j,1).Custom2 = 42;
        
        %j = 2;
        %s(j,1).Geometry = 'Polygon';
        %s(j,1).X = [597620 597630 597630 597620 597620];
        %s(j,1).Y = [170230 170230 170220 170220 170230];
        %s(j,1).BoundingBox = [min(s(j).X), min(s(j).Y); max(s(j).X), max(s(j).Y)];
        %s(j,1).Custom1 = 'cdf';
        %s(j,1).Custom2 = 43;
        
        s2 = struct;
        j = 1;
        s2(j,1).Geometry = 'Point';
        s2(j,1).X = [590610];
        s2(j,1).Y = [178220];
        s2(j,1).id = 1;
        
       
        %s(j,1).BoundingBox = [min(s(j,1).X), min(s(j,1).Y); max(s(j,1).X), max(s(j,1).Y)];
        %s(j,1).Custom1 = 'abc';
        %s(j,1).Custom2 = 42;
        
        % test if input and output have the same field values
        shapewrite(s2, 'point_oct.shp')
        
        

    if arg.Results.verbose
        
        fprintf('done!\n');
        pause(0.01);
        
    end
    
end