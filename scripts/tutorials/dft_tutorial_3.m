% DIGITAL FORESTRY TOOLBOX - TUTORIAL 3
%
% Other m-files required: LASread.m, elevationModels.m, treeStems.m
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2017b, GNU Octave 4.4.1 (configured for "x86_64-w64-mingw32")
%
% See also:
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory
% Website: http://mparkan.github.io/Digital-Forestry-Toolbox/
% Last revision: December 7, 2018
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund (WHFF, OFEV), project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details

clc
clear
close all

OCTAVE_FLAG = (exist('OCTAVE_VERSION', 'builtin') ~= 0); % determine if system is Matlab or GNU Octave

if OCTAVE_FLAG
    
    pkg load statistics
    pkg load image
    pkg load io
    pkg load mapping

end
    
%% Step 1 - Read the LAS file

% IMPORTANT: adjust the path to the input LAS file
pc = LASread('ge_2017_a.las');


%% Step 2 - Normalize the point cloud elevation

% compute the terrain model
cellSize = 0.5;
[models, refmat] = elevationModels([pc.record.x, pc.record.y, pc.record.z], ...
    pc.record.classification, ...
    'classTerrain', [2], ...
    'classSurface', [4,5], ...
    'cellSize', cellSize, ...
    'closing', 5, ...
    'interpolation', 'idw', ...
    'searchRadius', inf, ...
    'weightFunction', @(d) d^-3, ...
    'smoothingFilter', fspecial('gaussian', [2, 2], 0.8), ...
    'outputModels', {'terrain'}, ...
    'fig', false, ...
    'verbose', true);

[nrows, ncols] = size(models.terrain.values);

% subtract the terrain elevation from the point cloud elevation
P = round([pc.record.x - refmat(3,1), pc.record.y - refmat(3,2)] / refmat(1:2,:));
ind = sub2ind([nrows, ncols], P(:,1), P(:,2));
xyh = [pc.record.x, pc.record.y, pc.record.z - models.terrain.values(ind)];


%% Step 3 - Filter points by classification and return index

idxl_last = pc.record.return_number == pc.record.number_of_returns; % derniers echos
idxl_veg = ismember(pc.record.classification, [4, 5]); % classe végétation haute
idxl_filter = idxl_veg & idxl_last; % combiner les filtres


%% Step 4 - Detect stems

[label, xyh_stem] = treeStems(xyh, ...
    idxl_filter, ...
    'cellSize', 0.4, ...
    'bandWidth', 1.5, ...
    'verticalStep', 0.25, ...
    'searchRadius', 2, ...
    'minLength', 5, ...
    'verbose', true, ...
    'fig', true);


%% Step 5 - Export stem attributes to a CSV file

% IMPORTANT: adjust the path to the output CSV file
fid = fopen('ge_2017_a_stems.csv', 'w+'); % open file
fprintf(fid, 'X, Y, H\n'); % write header line
fprintf(fid, '%.2f, %.2f, %.2f\n', xyh_stem'); % write records
fclose(fid); % close file


%% Step 6 - Export stem attributes to an ESRI shapefile

S = struct;
for j = 1:size(xyh_stem,1)
    
    S(j,1).Geometry = 'Point';
    S(j,1).BoundingBox = [];
    S(j,1).ID = j;
    S(j,1).X = xyh_stem(j,1);
    S(j,1).Y = xyh_stem(j,2);
    S(j,1).H = xyh_stem(j,3);
     
end

% IMPORTANT: Octave users, please make sure you are using the latest versions of 
% the 'io' (2.4.12 or above) and 'mapping' (1.4.0 or above) packages. 
% Previous versions contain critical issues in the shapewrite function.
shapewrite(S, 'ge_2017_a_stems.shp'); % IMPORTANT: adjust the path to the output SHP file

