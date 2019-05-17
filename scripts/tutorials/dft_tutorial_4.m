% DIGITAL FORESTRY TOOLBOX - TUTORIAL 4
%
% Other files required: LASread.m, elevationModels.m, canopyPeaks.m, treeWatershed.m, ASCwrite.m, slicmex.mex
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2018b, GNU Octave 5.1.0 (configured for "x86_64-w64-mingw32")
%
% See also:
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan
% Website: http://mparkan.github.io/Digital-Forestry-Toolbox/
% Last revision: May 17, 2019
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


%% Step 2 - Compute a raster Canopy Height Model (CHM)

cellSize = 0.5;
[models, refmat] = elevationModels([pc.record.x, pc.record.y, pc.record.z], ...
    pc.record.classification, ...
    'classTerrain', [2, 9, 16], ...
    'classSurface', [3, 4, 5], ...
    'cellSize', cellSize, ...
    'closing', 3, ...
    'smoothingFilter', [], ...
    'outputModels', {'terrain', 'surface', 'height'}, ...
    'fig', false, ...
    'verbose', true);

% convert map coordinates to image coordinates
[nrows, ncols] = size(models.terrain.values);
P = [pc.record.x - refmat(3,1), pc.record.y - refmat(3,2)] / refmat(1:2,:);
row = round(P(:,1));
col = round(P(:,2));
ind = sub2ind([nrows, ncols], row, col);


%% Step 3 - Filter the canopy points

idxl_veg = ismember(pc.record.classification, [4,5]);

z_margin = 1;
idxl_canopy = pc.record.z >= (models.surface.values(ind) - z_margin);
idxl_filter = idxl_veg & idxl_canopy;


%% Step 4 - Compute an intensity raster

I = accumarray([row(idxl_filter), col(idxl_filter)], double(pc.record.intensity(idxl_filter)), [nrows, ncols], @mean, 0);

% plot the intensity map
figure('Color', [1,1,1])
imagesc(I);
axis equal tight
colorbar
caxis(quantile(I(:), [0.02, 0.98]))

% export the intensity map to an ARC/INFO ASCII grid file
ASCwrite('ge_2017_a_intensity.asc', ...
    I, ...
    refmat, ...
    'precision', 2, ...
    'noData', -99999, ...
    'verbose', true);


%% Step 5 - Segmentation

segmentation = 'slic'; % specify 'slic' or 'watershed'

switch segmentation
    
    case 'slic' % Superpixel segmentation
        
        % rescale values to an 8 bit range (0-255)
        I_n = (I - quantile(I(:), 0.02)) ./ diff(quantile(I(:), [0.02, 0.98]));
        I_n = uint8(I_n * 255);
        
        % before running slicmex, you have to compile it with the following command:
        % mex slicmex.c
        r = 2.5; % target radius of one superpixel (in map units)
        a = pi*r^2; % target area of one superpixel (in map units)
        n_superpixels = round(nrows * ncols * cellSize^2 / a); % number of superpixels
        [L, nlabels] = slicmex(I_n, n_superpixels, 100);
        
        % Note: Matlab users who have the "image processing" toolbox can also use the superpixels() function
        % [L, nlabels] = superpixels(I_n, n_superpixels, 'Compactness', 75);
        
    case 'watershed' % Marker controlled watershed segmentation
        
        [peaks_crh, ~] = canopyPeaks(models.height.values, ...
            refmat, ...
            'minTreeHeight', 2, ...
            'smoothingFilter', fspecial('gaussian', [3 3], 1), ...
            'searchRadius', @(h) 0.28 * h^0.59, ...
            'fig', false, ...
            'verbose', true);
        
        [L, ~] = treeWatershed(models.height.values, ...
            'markers', peaks_crh(:,1:2), ...
            'minHeight', 1, ...
            'removeBorder', false, ...
            'fig', false, ...
            'verbose', true);
        
end


%% Step 6 - compute segment level mean intensity and mean height

I_s = accumarray(L(:)+1, I(:), [nrows*ncols,1], @nanmean, nan);
I_s = I_s(L+1);

H_s = accumarray(L(:)+1, models.height.values(:), [nrows*ncols,1], @nanmax, single(nan));
H_s = H_s(L+1);

% threshold minimum segment height
h_min = 3;
I_s(H_s <= h_min) = nan;

% remove null values
I_s(I_s == 0) = nan;

% plot
figure
imagesc(I_s)
axis equal tight
colorbar
caxis(quantile(I_s(:), [0.02, 0.98]))


%% Step 7 - convert the intensity image to a classification map

% plot the intensity histogram
figure
hist(I_s(~isnan(I_s)), 100)
xlabel('intensity')
ylabel('count')

% notice the bimodal distribution, with a threshold around 17'000 for the
% persistant foliage

% threshold the segment scale intensity map
[~, C] = histc(I_s(:), [7500, 17000, inf]);

% plot the classification map (0 = non-forest, 1 = deciduous, 2 = persistant)
C = reshape(C, nrows, ncols);
figure
imagesc(C)
axis equal tight


%% Step 8 - filter small patches with by applying morphological reconstruction on the eroded

se = strel("disk", 3, 0); % circular structuring element with a radius of 3 pixels 
J = imreconstruct(imerode(C == 2, se), C == 2);
C(J ~= (C == 2)) = 1;

% plot the classification map
figure
imagesc(C)
axis equal tight

% export the classification map to an ARC/INFO ASCII grid file
ASCwrite('ge_2017_a_class.asc', ...
    C, ...
    refmat, ...
    'precision', 0, ...
    'noData', -99999, ...
    'verbose', true);
