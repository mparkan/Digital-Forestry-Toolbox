function t = laserTimeFormat(gps_week_number, gps_week_time, format)
%LASERTIMEFORMAT - Convert GPS week time to GPS satellite time (i.e. number of seconds since 1980-1-6), GPS adjusted time
% (i.e. GPS satellite time - 10^9) or UTC time.
% T = LASERTIMEFORMAT(GPS_WEEK_NUMBER, GPS_WEEK_TIME, FORMAT, ...)
% 
% For more information, check:
% U.S. Naval Observatory (USNO), Time Service - GPS week, http://tycho.usno.navy.mil/gps_week.html
% NOAA, Continuously Operating Reference Station (CORS) - GPS calendar, http://www.ngs.noaa.gov/CORS/Gpscal.shtml
% QPS, UTC to GPS Time Correction, https://confluence.qps.nl/display/KBE/UTC+to+GPS+Time+Correction
%
% Syntax:  time = laserTimeFormat(gps_week_number, gps_week_time, flightline_id, ...)
%
% Inputs:
%    gps_week_number - Nx1 numeric vector, GPS week number for each target point
%
%    gps_week_time - Nx1 numeric vector, GPS week time for each target point
%
%    format (optional, default: 'adjusted') - a string with value 'satellite', 'adjusted' or
%    'utc', indicating the output time format 
%
% Outputs:
%    t - Nx1 float vector
%
% Example:
%    pc = LASread(pointcloud.las);
%    trajectory = shaperead(trajectory.shp);
%
%    trajectory_flightline_id = [trajectory.Line_ID]';
%    trajectory_gps_week_number = [trajectory.GPS_WN]';
        
%    [~, idxn_flightline_id] = ismember(pc.record.point_source_id, trajectory_flightline_id);
%    gps_week_number = trajectory_gps_week_number(idxn_flightline_id);
%    gps_week_time = pc.record.gps_time;
%
%    time = laserTimeFormat(gps_week_number, ...
%        gps_week_time, ...
%        'adjusted');
%
% Other m-files required: gps2utc.m
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2016b
%
% See also: gps2utc.m
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory (LASIG)
% Website: http://mparkan.github.io/Digital-Forestry-Toolbox/
% Last revision: January 18, 2017
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund (WHFF, OFEV), project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'gps_week_time', @isnumeric);
addRequired(arg, 'gps_week_number', @isnumeric);
% addRequired(arg, 'flightline_id', @isnumeric);
% addRequired(arg, 'gps_week_number', @(x) (size(x,2) == 1) && isnumeric(x));
% addRequired(arg, 'flightline_id', @(x) (size(x,2) == 1) && isnumeric(x));
addOptional(arg, 'format', 'adjusted', @(x) ismember(x, {'satellite', 'adjusted', 'utc'}));

parse(arg, gps_week_time, gps_week_number, format);


%% reformat time

% update LAS GPS time
switch format
    
    case 'satellite' % GPS satellite time
        
        % GPS satellite time (i.e. number of seconds since 1980-1-6)
        t = gps_week_number * 7 * 24 * 3600 + gps_week_time;
        
    case 'adjusted' % GPS adjusted time
        
        % GPS adjusted time (i.e. GPS satellite time - 10^9, the offset moves the time back to
        % near zero to improve floating point resolution.)
        t = gps_week_number * 7 * 24 * 3600 - 10^9 + gps_week_time;
        
    case 'utc' % Universal Time Count (UTC)
        
        % GPS, Global Positioning System time, is the atomic time scale implemented by the atomic clocks in the GPS ground 
        % control stations and the GPS satellites themselves. GPS time was zero at 0h 6-Jan-1980 and since it is not 
        % perturbed by leap seconds GPS is now ahead of UTC.  
        gps_sat_time = gps_week_number * 7 * 24 * 3600 + gps_week_time;
        t = gps2utc(gps_sat_time, 'verbose', true);
        
end

