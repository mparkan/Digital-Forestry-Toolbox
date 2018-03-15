function intensity_corr = laserIntensityCorrection(intensity_target, xyzt_target, xyzt_sensor, varargin)
%LASERINTENSITYCORRECTION - Corrects return intensity for the sensor to target range using the theoretical 
% model proposed in [1,2,3]: intensity_corr = intensity_target * (range^exponent / referenceRange^exponent)
%
% INTENSITY = LASERINTENSITYCORRECTION(INTENSITY_TARGET, XYZT_TARGET,
% XYZT_SENSOR, ...) corrects the raw target intensity values specified in
% INTENSITY_TARGET for the range. The range is computed from the position and acquision
% time of the target points XYZT_TARGET and the sensor trajectory XYZT_SENSOR.
% The correction uses the simple theoretical model proposed in [1,2,3] which requires 
% the sensor to target range, a reference distance and a chosen exponent (typically = 2)
% 
% [1] A. G. Kashani, M. J. Olsen, C. E. Parrish, et N. Wilson, A Review of LIDAR Radiometric Processing: 
% From Ad Hoc Intensity Correction to Rigorous Radiometric Calibration, 
% Sensors, vol. 15, n°11, p. 28099-28128, nov. 2015.
% 
% [2] M. Starek, B. Luzum, R. Kumar, K.C Slatton, Normalizing LIDAR Intensities, GEM Center
% Report No. Rep_2006-12-001, University of Florida, Gainesville, FL, USA, 2006.
% 
% [3] I. Korpela, H. O. Ørka, J. Hyyppä, V. Heikkinen, et T. Tokola, 
% Range and AGC normalization in airborne discrete-return LiDAR intensity data for forest canopies, 
% ISPRS Journal of Photogrammetry and Remote Sensing, vol. 65, n°4, p. 369?379, jul. 2010.
%
% Syntax: intensity_corr = laserIntensityCorrection(intensity_target, xyzt_target, xyzt_sensor, ...)
%
% Inputs:
%    intensity_target - Nx1 numeric vector, raw intensity of the target
%    points
%       
%    xyzt_target - Nx4 numeric matrix, xyz coordinates and acquisition time
%    t of the target points. The time format must be the same as for the sensor trajectory 
%    (preferablly GPS satellite time or adjusted GPS time).  
%    
%    xyzt_sensor - Mx4 numeric matrix, xyz coordinates and acquisition time
%    t of the sensor trajectory. The time format must be the same as for the target points.  
%    (preferablly GPS satellite time or adjusted GPS time).  
%    
%    referenceRange (optional, default: minimum sensor-target range) - numeric value, reference
%    range of the theoretical correction model
%    
%    exponent (optional, default: 2) - numeric value, exponent of the
%    theoretical correction model
%    
%    verbose (optional, default: true) - boolean value, verbosiy switch
%    
%    fig (optional, default: false) - boolean value, switch to plot figures
%    
% Outputs:
%    intensity_corr - Nx1 numeric vector, range corrected intensity of the target
%    points
%
% Example:
%
%    pc = LASread('E:\als\mypoints.las', false, true);
%    xyzt_target = [pc.record.x, pc.record.y, pc.record.z, pc.record.gps_time];
%
%    traj = TRJread('E:\trajectories\4242.trj', ...
%    'headerOnly', false, ...
%    'verbose', true);
%    xyzt_sensor = [traj.record.X, traj.record.Y, traj.record.Z, traj.record.Time];
%
%    intensity_corr = laserIntensityCorrection(pc.record.intensity, ...
%    xyzt_target, ...
%    xyzt_sensor, ...
%    'referenceRange', 1000, ...
%    'exponent', 2, ...
%    'verbose', true, ...
%    'fig', true);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2017b, GNU Octave 4.2.1 (configured for "x86_64-w64-mingw32")
%
% See also:
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory (LASIG)
% Website: http://mparkan.github.io/Digital-Forestry-Toolbox/
% Last revision: March 15, 2018
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund (WHFF, OFEV), project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'intensity_target', @isnumeric);
addRequired(arg, 'xyzt_target', @isnumeric);
addRequired(arg, 'xyzt_sensor', @isnumeric);
addParameter(arg, 'referenceRange', [], @(x) isnumeric(x) && (numel(x) == 1));
addParameter(arg, 'exponent', 2, @(x) isnumeric(x) && (numel(x) == 1));
addParameter(arg, 'fig', false, @(x) islogical(x) && (numel(x) == 1));
addParameter(arg, 'verbose', true, @(x) islogical(x) && (numel(x) == 1));

parse(arg, intensity_target, xyzt_target, xyzt_sensor, varargin{:});


%% interpolate trajectories

if arg.Results.verbose
    
    tic
    fprintf('interpolating sensor trajectory...');
    
end

xi_sensor = interp1(xyzt_sensor(:,4), xyzt_sensor(:,1), xyzt_target(:,4), 'linear');
yi_sensor = interp1(xyzt_sensor(:,4), xyzt_sensor(:,2), xyzt_target(:,4), 'linear');
zi_sensor = interp1(xyzt_sensor(:,4), xyzt_sensor(:,3), xyzt_target(:,4), 'linear');

if arg.Results.verbose
    
    fprintf('done!\n');
    toc
    
end


%% compute sensor-target range

if arg.Results.verbose
    
    tic
    fprintf('computing sensor-target range...');
    
end

d = sqrt(sum(([xi_sensor, yi_sensor, zi_sensor] - xyzt_target(:,1:3)).^2,2));

if arg.Results.verbose
    
    fprintf('done!\n');
    toc
    
end


%% correct intensity for range

if arg.Results.verbose
    
    tic
    fprintf('correcting intensity for range...');
    
end

intensity_target = double(intensity_target);

if isempty(arg.Results.referenceRange)
    
    %referenceRange = median(d);
    referenceRange = min(d);
    
else
    
    referenceRange = arg.Results.referenceRange;
    
end

intensity_corr = intensity_target .* (d.^arg.Results.exponent / referenceRange^arg.Results.exponent); 

if arg.Results.verbose
    
    fprintf('done!\n');
    toc
    
end


%% plot figures

if arg.Results.fig
    
    figure
    scatter3(xyzt_sensor(:,1), ...
        xyzt_sensor(:,2), ...
        xyzt_sensor(:,3), ...
        6, ...
        [0, 0, 0], ...
        'Marker', '.');
    hold on
    scatter3(xi_sensor, ...
        yi_sensor, ...
        zi_sensor, ...
        6, ...
        [1, 0, 0], ...
        'Marker', '.');
    axis equal tight
    title('Interpolated sensor trajectory')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
end
