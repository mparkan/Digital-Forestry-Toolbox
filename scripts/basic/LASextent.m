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
%    outputFilepath (optional, default: []) - The path to the output shapefile
%
%    method (optional, default: 'bbox') - The type of extent ('header': use values in header, 
%    'bbox': compute rectangular bounding box, 'concavehull', computes the concave hull, 
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
%         'method', 'concavehull', ...
%         'fig', true, ...
%         'verbose', false);
%
% Other m-files required: LASread.m, LASwrite.m, subsample.m, gps2utc.m
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2016b
% 
% See also: LASCLIP, LASMERGE, LASTILE
% 
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory (LASIG)
% Website: http://mparkan.github.io/Digital-Forestry-Toolbox/
% Last revision: March 13, 2017
% Acknowledgments: This work was supported by the Swiss Forestry and Wood
% Research Fund, WHFF (OFEV) - project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'points', @(x) iscell(x) || ischar(x));
addOptional(arg, 'outputFilepath', [], @(x) ischar(x) || isempty(x));
addParameter(arg, 'method', 'bbox', @(x) ismember(x, {'header', 'bbox', 'concavehull', 'convexhull'}));
% addParameter(arg, 'resolution', 1, @(x) isnumeric(x) && (numel(x) == 1) && (x >= 0));
addParameter(arg, 'fig', false, @(x) islogical(x) && (numel(x) == 1));
addParameter(arg, 'verbose', false, @(x) islogical(x) && (numel(x) == 1));

parse(arg, points, varargin{:});


%% list LAS input files

if arg.Results.verbose
   
    fprintf('creating file list...');
    
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
    
end
    

%% compute metadata

n_files = length(filelist);
extent = table();
warning('off');

for j = 1:n_files
    
    extent.FILEPATH{j,1} = [filelist(j).path, filesep, filelist(j).name, filelist(j).ext];
    extent.FILENAME{j,1} = filelist(j).name;
    
    % read LAS file
    switch arg.Results.method
        
        case 'header'
            
            pc = LASread(extent.FILEPATH{j,1}, true, false);
            
            % read extrema
            extent.XMIN(j,1) = pc.header.min_x;
            extent.XMAX(j,1) = pc.header.max_x;
            extent.YMIN(j,1) = pc.header.min_y;
            extent.YMAX(j,1) = pc.header.max_y;
            extent.ZMIN(j,1) = pc.header.min_z;
            extent.ZMAX(j,1) = pc.header.max_z;
            
        otherwise
            
            if arg.Results.verbose
                
                fprintf('reading file %u/%u...', j, n_files);
                
            end
    
            pc = LASread(extent.FILEPATH{j,1}, false, false);
            
            if arg.Results.verbose
                
                fprintf('done!\n');
                
            end
            
            
            if arg.Results.verbose
                
                fprintf('extracting metadata of file %u/%u...', j, n_files);
                
            end
            
            % compute extrema
            extent.XMIN(j,1) = min(pc.record.x);
            extent.XMAX(j,1) = max(pc.record.x);
            extent.YMIN(j,1) = min(pc.record.y);
            extent.YMAX(j,1) = max(pc.record.y);
            extent.ZMIN(j,1) = min(pc.record.z);
            extent.ZMAX(j,1) = max(pc.record.z);
            
            % compute point density (per unit area)
            [~, idxn_unit] = subsample([pc.record.x, pc.record.y], ...
                'method', 'grid', ...
                'resolution', 1, ...
                'fig', false, ...
                'verbose', false);
            
            extent.PDENSITY(j,1) = median(accumarray(idxn_unit, idxn_unit, [], @numel));
            
            % compute date metadata (if GPS time is in satellite format)
            if pc.header.global_encoding_gps_time_type == 1
                
                utc_time = gps2utc(pc.record.gps_time + 10^9, ...
                    'verbose', false);
                
                [Y,M,D,~,~,~] = datevec(utc_time);
                
                unique_date = sort(unique(datenum(Y,M,D)));
                start_date = min(unique_date);
                end_date = max(unique_date);
                
                % number of surveys
                extent.NSURVEYS(j,1) = length(unique_date);
                
                [extent.YSTART(j,1), extent.MSTART(j,1), extent.DSTART(j,1),~,~,~] = datevec(start_date);
                [extent.YEND(j,1), extent.MEND(j,1), extent.DEND(j,1),~,~,~] = datevec(end_date);
                
                extent.DATES{j,1} = strjoin(cellstr(datestr(unique_date, 'yyyy/mm/dd')), ',');

            end
            
            % subsample point cloud
            xrange = extent.XMAX(j) - extent.XMIN(j);
            yrange = extent.YMAX(j) - extent.YMIN(j);
            max_range = max([xrange, yrange]);
            resolution = min(3, max_range / 100);
            
            [xy_s, idxn_cell] = subsample([pc.record.x, pc.record.y], ...
                'method', 'grid', ...
                'resolution', resolution, ...
                'fig', false, ...
                'verbose', false);
            
            if arg.Results.verbose
                
                fprintf('done!\n');
                
            end

    end
    
    if arg.Results.verbose
        
        fprintf('computing spatial extent of file %u/%u...', j, n_files);
        
    end
            
    % compute spatial extent
    switch arg.Results.method
        
        case 'header'
            
            extent.X{j,1} = [extent.XMIN(j), extent.XMAX(j), extent.XMAX(j), extent.XMIN(j), extent.XMIN(j)];
            extent.Y{j,1} = [extent.YMAX(j), extent.YMAX(j), extent.YMIN(j), extent.YMIN(j), extent.YMAX(j)];
        
        case 'bbox'
            
            extent.X{j,1} = [extent.XMIN(j), extent.XMAX(j), extent.XMAX(j), extent.XMIN(j), extent.XMIN(j)];
            extent.Y{j,1} = [extent.YMAX(j), extent.YMAX(j), extent.YMIN(j), extent.YMIN(j), extent.YMAX(j)];
        
        case 'concavehull'
            
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

            extent.X{j,1} = shp_hr.Points(idxn_boundary,1);
            extent.Y{j,1} = shp_hr.Points(idxn_boundary,2);

        case 'convexhull'
            
            idxn_boundary = convhull(pc.record.x, pc.record.y);

            extent.X{j,1} = pc.record.x(idxn_boundary);
            extent.Y{j,1} = pc.record.y(idxn_boundary);
            
    end
    
    if arg.Results.verbose
        
        fprintf('done!\n');
        
    end
    
end

pause(0.1);


%% compute metadata fields


%% export to shapefile

if ~isempty(arg.Results.outputFilepath)
    
    if arg.Results.verbose
        
        fprintf('\nwriting result to "%s"...', arg.Results.outputFilepath);
        
    end
    
    extent.Geometry = repmat({'Polygon'}, height(extent), 1);
    shapewrite(table2struct(extent), arg.Results.outputFilepath);

    if arg.Results.verbose
        
        fprintf('done!\n');
        
    end
    
end


%% plot figures

if arg.Results.fig
    
    figure
    
    for j = 1:n_files
        
        plot(extent.X{j,1}, ...
            extent.Y{j,1}, ...
            'r-')
        hold on
        
    end
    
    axis equal tight
    title('Extents')
    xlabel('x')
    ylabel('y')
    
end
