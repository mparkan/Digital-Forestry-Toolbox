% DIGITAL FORESTRY TOOLBOX - TUTORIAL 1 
%
% Other m-files required: LASread.m, gps2utc.m, crossSection.m, LASclip.m,
% LASextent.m
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
% Last revision: March 23, 2018
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund (WHFF, OFEV), project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details

clc
clear
close all

OCTAVE_FLAG = (exist('OCTAVE_VERSION', 'builtin') ~= 0); % determine if system is Matlab or GNU Octave

if OCTAVE_FLAG
    
    pkg load statistics
    pkg load image
    pkg load mapping

end

    
%% Step 1 - Reading the LAS file

% IMPORTANT: adjust the path to the input LAS file
pc = LASread('zh_2014_a.las');


%% Step 2 - Subsetting

% define a Region Of Interest (ROI)
x_roi = [2699597, 2699584, 2699778, 2699804];
y_roi = [1271166, 1271307, 1271341, 1271170];

% create ROI logical index (with inclusion test)
idxl_roi = inpolygon(pc.record.x, pc.record.y, x_roi, y_roi);

% create class subset logical index
idxl_class = ismember(pc.record.classification, [4,5]);

% combine subsetting indices
idxl_sample = idxl_roi & idxl_class;

% apply subsetting indices
xyz_s = [pc.record.x(idxl_sample), pc.record.y(idxl_sample), pc.record.z(idxl_sample)];
intensity_s = single(pc.record.intensity(idxl_sample));


%% Step 3 - Visualization

figure
scatter3(xyz_s(:,1), ...
xyz_s(:,2), ...
xyz_s(:,3), ...
6, ...
intensity_s, ...
'Marker', '.');
colorbar
caxis(quantile(intensity_s, [0.01, 0.99]))
axis equal tight
title('Return intensity')
xlabel('x')
ylabel('y')
ylabel('z')


%% Step 4 - Checking the acquisition dates

utc_time = gps2utc(pc.record.gps_time + 10^9, 'verbose', true);


%% Step 5 - Clipping directly to a LAS file

% define Region Of Interest (ROI)
x_roi = [2699597, 2699584, 2699778, 2699804];
y_roi = [1271166, 1271307, 1271341, 1271170];

% IMPORTANT: adjust the path to the input SHP file and output LAS file
LASclip('zh_2014_a.las', ...
    [x_roi', y_roi'], ...
    'zh_2014_a_clip.las', ...
    'verbose', true);


%% Step 6 - Computing the spatial extent of a LAS file (Matlab only)

% IMPORTANT: adjust the path to the input LAS file and output SHP file
extent = LASextent('zh_2014_a.las', ...
    'zh_2014_a_extent.shp', ...
    'method', 'convexhull', ...
    'fig', false, ...
    'verbose', true);


%% Step 7 - Extracting a cross-section from a 3D point cloud

width = 2;
p0 = [2699582, 1271120];
p1 = [2699884, 1271333];
[~, ~, profile] = crossSection([pc.record.x, pc.record.y, pc.record.z], ...
    width, ...
    p0, ...
    p1, ...
    'verbose', true, ...
    'fig', true);
