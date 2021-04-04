% DIGITAL FORESTRY TOOLBOX - TUTORIAL 2
%
% Other m-files required: LASread.m, LASwrite.m, rasterize.m, elevationModels.m, canopyPeaks.m,
% treeWatershed.m, topoColor.m
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2020b, GNU Octave 6.2.0 (configured for "x86_64-w64-mingw32")
%
% See also:
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory
% Website: http://mparkan.github.io/Digital-Forestry-Toolbox/
% Last revision: April 4, 2021
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund (WHFF, OFEV), project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details

clc
clear
close all

OCTAVE_FLAG = (exist('OCTAVE_VERSION', 'builtin') ~= 0); % determine if system is Matlab or GNU Octave

if OCTAVE_FLAG
    
    pkg load statistics
    pkg load image
    more off
    
end


%% Step 1 - Reading the LAS file

% IMPORTANT: adjust the path to the input LAS file
pc = LASread('zh_2014_a.las');


%% Step 2 - Computing a raster Canopy Height Model (CHM)

cellSize = 0.8;
[models, refmat] = elevationModels([pc.record.x, pc.record.y, pc.record.z], ...
    pc.record.classification, ...
    'classTerrain', [2], ...
    'classSurface', [4,5], ...
    'cellSize', cellSize, ...
    'interpolation', 'idw', ... 
    'searchRadius', inf, ...
    'weightFunction', @(d) d^-3, ...
    'smoothingFilter', fspecial('gaussian', [3, 3], 0.8), ...
    'outputModels', {'terrain', 'surface', 'height'}, ...
    'fig', true, ...
    'verbose', true);


%% Step 3 - Exporting the elevation models to ARC/INFO ASCII grid file

% export the Digital Terrain Model (DTM) to an ARC/INFO ASCII grid file
ASCwrite('zh_2014_a_dtm.asc', ...
    models.terrain.values, ...
    refmat, ...
    'precision', 2, ...
    'noData', -99999, ...
    'verbose', true);
            
% export the Digital Surface Model (DSM) to an ARC/INFO ASCII grid file
ASCwrite('zh_2014_a_dsm.asc', ...
    models.surface.values, ...
    refmat, ...
    'precision', 2, ...
    'noData', -99999, ...
    'verbose', true);

% export the Digital Height Model (DHM) to an ARC/INFO ASCII grid file
ASCwrite('zh_2014_a_dhm.asc', ...
    models.height.values, ...
    refmat, ...
    'precision', 2, ...
    'noData', -99999, ...
    'verbose', true);


%% Step 4 - Tree top detection

[peaks_crh, ~] = canopyPeaks(models.height.values, ...
    refmat, ...
    'minTreeHeight', 2, ...
    'searchRadius', @(h) 0.28 * h^0.59, ...
    'fig', true, ...
    'verbose', true);


%% Step 5 - Marker controlled watershed segmentation

[label_2d, colors] = treeWatershed(models.height.values, ...
    'markers', peaks_crh(:,1:2), ...
    'minHeight', 1, ...
    'removeBorder', true, ...
    'fig', true, ...
    'verbose', true);


%% Step 6 - Computing segment metrics from the label matrix

% IMPORTANT: some of the metrics in regionprops() are currently only available in Matlab
metrics_2d = regionprops(label_2d, models.height.values, ...
    'Area', 'Centroid', 'MaxIntensity');


%% Step 7 - Transferring 2D labels to the 3D point cloud

idxl_veg = ismember(pc.record.classification, [4,5]);

% convert map coordinates (x,y) to image coordinates (column, row)
RC = [pc.record.x - refmat(3,1), pc.record.y - refmat(3,2)] / refmat(1:2,:);
RC(:,1) = round(RC(:,1)); % row
RC(:,2) = round(RC(:,2)); % column
ind = sub2ind(size(label_2d), RC(:,1), RC(:,2));

% transfer the label
label_3d = label_2d(ind);
label_3d(~idxl_veg) = 0;
[label_3d(label_3d ~= 0), ~] = grp2idx(label_3d(label_3d ~= 0));

% transfer the color index
color_3d = colors(ind);
color_3d(~idxl_veg) = 1;

% define a colormap
cmap = [0, 0, 0;
    166,206,227;
    31,120,180;
    178,223,138;
    51,160,44;
    251,154,153;
    227,26,28;
    253,191,111;
    255,127,0;
    202,178,214;
    106,61,154;
    255,255,153;
    177,89,40] ./ 255;


%% Step 8 - plot the colored points cloud (vegetation points only)

% Due to 3D plotting performance issues, the display of large points
% clouds is currently discouraged in Octave
if ~OCTAVE_FLAG
    
    % plot all the segments (Matlab only)
    figure('Color', [1,1,1])
    scatter3(pc.record.x(idxl_veg), ...
        pc.record.y(idxl_veg), ...
        pc.record.z(idxl_veg), ...
        6, ...
        color_3d(idxl_veg), ...
        'Marker', '.')
    axis equal tight
    colormap(cmap)
    caxis([1, size(cmap,1)])
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    % plot a single segment (Matlab only)
    idxl_sample = (label_3d == 42);
    
    figure
    scatter3(pc.record.x(idxl_sample), ...
        pc.record.y(idxl_sample), ...
        pc.record.z(idxl_sample), ...
        12, ...
        pc.record.intensity(idxl_sample), ...
        'Marker', '.');
    colorbar
    axis equal tight
    title('Return intensity')
    xlabel('x')
    ylabel('y')
    ylabel('z')
    
end


%% Step 9 - Computing segment metrics from the labelled point cloud

% compute metrics
[metrics_3d, fmt, idxl_scalar] = treeMetrics(label_3d, ...
    [pc.record.x, pc.record.y, pc.record.z], ...
    pc.record.intensity, ...
    pc.record.return_number, ...
    pc.record.number_of_returns, ...
    nan(length(pc.record.x), 3), ...
    models.terrain.values, ...
    refmat, ...
    'metrics', {'UUID', 'XPos', 'YPos', 'ZPos', 'H','BBOX2D', 'XCVH2D', 'YCVH2D', 'CVH2DArea', 'IQ50'}, ...
    'intensityScaling', true, ...
    'alphaMin', 1.5, ...
    'verbose', true);

% list field names
sfields = fieldnames(metrics_3d);


%% Step 10 - Exporting the segment metrics to a CSV file

% IMPORTANT: adjust the path to the output CSV file (otherwise it will be created in the current folder)
fid = fopen('zh_2014_a_seg_metrics.csv', 'w+'); % open file
fprintf(fid, '%s\n', strjoin(sfields(idxl_scalar), ',')); % write header line
C = struct2cell(rmfield(metrics_3d, sfields(~idxl_scalar))); % convert structure to cell array
fprintf(fid, [strjoin(fmt(idxl_scalar), ','), '\n'], C{:}); % write cell array to CSV file
fclose(fid); % close file


%% Step 11 - Exporting the segment points (and metrics) to a SHP file

% duplicate the metrics structure (scalar fields only)
S1 = rmfield(metrics_3d, sfields(~idxl_scalar));

% add the geometry type
[S1.Geometry] = deal('Point');

% add the X coordinates of the polygons
[S1.X] = metrics_3d.XPos;

% add the Y coordinates of the polygons
[S1.Y] = metrics_3d.YPos;

% write non-scalar structure to SHP file
% IMPORTANT: the shapewrite() function included here is currently 
% not compatible with Matlab. Matlab users should use the shapewrite()
% function from the offical Matlab mapping toolbox instead.
shapewrite(S1, 'zh_2014_a_seg_points.shp'); % IMPORTANT: adjust the path to the output SHP file (otherwise it will be created in the current folder)
clear S1


%% Step 12 - Exporting the segment polygons (and metrics) to a SHP file

% duplicate the metrics structure (scalar fields only)
S2 = rmfield(metrics_3d, sfields(~idxl_scalar));

% add the geometry type
[S2.Geometry] = deal('Polygon');

% add the X coordinates of the polygons
[S2.X] = metrics_3d.XCVH2D;

% add the Y coordinates of the polygons
[S2.Y] = metrics_3d.YCVH2D;

% add the bounding boxes [minX, minY; maxX, maxY]
[S2.BoundingBox] = metrics_3d.BBOX2D;

% write non-scalar structure to SHP file
% IMPORTANT: the shapewrite() function included here is currently 
% not compatible with Matlab. Matlab users should use the shapewrite() 
% function from the offical Matlab mapping toolbox instead.
shapewrite(S2, 'zh_2014_a_seg_polygons.shp'); % IMPORTANT: adjust the path to the output SHP file (otherwise it will be created in the current folder)
clear S2


%% Step 13 - Exporting the labelled and colored point cloud to a LAS file

% duplicate the source file
r = pc;

% rescale the RGB colors to 16 bit range and add them to the point record
rgb = uint16(cmap(color_3d,:) * 65535);
r.record.red = rgb(:,1);
r.record.green = rgb(:,2);
r.record.blue = rgb(:,3);

% add the "label" field to the point record (as an uint32 field)
r.record.label = uint32(label_3d);

% add the "label" uint32 field metadata in the variable length records
% check the ASPRS LAS 1.4 specification for details about the meaning of the fields
% https://www.asprs.org/a/society/committees/standards/LAS_1_4_r13.pdf
vlr = struct;
vlr.value.reserved = 0;
vlr.value.data_type = 5;
vlr.value.options.no_data_bit = 0;
vlr.value.options.min_bit = 0;
vlr.value.options.max_bit = 0;
vlr.value.options.scale_bit = 0;
vlr.value.options.offset_bit = 0;
vlr.value.name = 'label';
vlr.value.unused = 0;
vlr.value.no_data = [0; 0; 0];
vlr.value.min = [0; 0; 0];
vlr.value.max = [0; 0; 0];
vlr.value.scale = [0; 0; 0];
vlr.value.offset = [0; 0; 0];
vlr.value.description = 'LABEL';

vlr.reserved = 43707;
vlr.user_id = 'LASF_Spec';
vlr.record_id = 4;
vlr.record_length_after_header = length(vlr.value) * 192;
vlr.description = 'Extra bytes';

% append the new VLR to the existing VLR
if isfield(r, 'variable_length_records')
    
    r.variable_length_records(length(r.variable_length_records)+1) = vlr;
    
else
    
    r.variable_length_records = vlr;
    
end

% if necessary, adapt the output record format to add the RGB channel
switch pc.header.point_data_format_id
    
    case 1 % 1 -> 3
        
        recordFormat = 3;
        
    case 4 % 4 -> 5
        
        recordFormat = 5;
        
    case 6 % 6 -> 7
        
        recordFormat = 7;
        
    case 9 % 9 -> 10
        
        recordFormat = 10;
        
    otherwise % 2,3,5,7,8,10
        
        recordFormat = pc.header.point_data_format_id;
        
end

% write the LAS 1.4 file
% IMPORTANT: adjust the path to the output LAS file
LASwrite(r, ...
    'zh_2014_a_ws_seg.las', ...
    'version', 14, ...
    'guid', lower(strcat(dec2hex(randi(16,32,1)-1)')), ...
    'systemID', 'SEGMENTATION', ...
    'recordFormat', recordFormat, ...
    'verbose', true);

% you can optionally read the exported file and check it has the 
% RGB color and label records
% IMPORTANT: adjust the path to the input LAS file
r2 = LASread('zh_2014_a_ws_seg.las');
