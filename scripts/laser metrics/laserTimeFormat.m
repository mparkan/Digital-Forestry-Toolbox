function time = laserTimeFormat(gps_week_time, point_source_id, gps_week_number, flightline_id, format)
%LASERTIMEFORMAT - Convert GPS week time to GPS satellite time (i.e. number of seconds since 1980-1-6), GPS adjusted time
% (i.e. GPS satellite time - 10^9) or UTC time.
% TIME = LASERTIMEFORMAT(GPS_WEEK_TIME, POINT_SOURCE_ID, GPS_WEEK_NUMBER, FLIGHTLINE_ID, FORMAT)
% 
% For more information, check:
% U.S. Naval Observatory (USNO), Time Service - GPS week, http://tycho.usno.navy.mil/gps_week.html
% NOAA, Continuously Operating Reference Station (CORS) - GPS calendar, http://www.ngs.noaa.gov/CORS/Gpscal.shtml
% QPS, UTC to GPS Time Correction, https://confluence.qps.nl/display/KBE/UTC+to+GPS+Time+Correction
%
% Syntax:  time = laserTimeFormat(gps_week_time, point_source_id, gps_week_number, flightline_id, format)
%
% Inputs:
%    gps_week_time - Nx1 numeric vector, GPS week time for each target point
%        
%    point_source_id - Nx1 numeric vector, flight line identifier for each target point
%
%    gps_week_number - Mx1 numeric vector, GPS week number for each unique flight line identifier 
%
%    flightline_id - Mx1 numeric vector, unique flight line identifier for each target point
%
%    format (optional, default: 'satellite') - a string with value 'satellite', 'adjusted' or
%    'utc', indicating the output time format 
%
% Outputs:
%    time - Nx1 float vector
%
% Example:
%   time = laserTimeFormat(gps_week_time, point_source_id, gps_week_number, flightline_id, 'adjusted');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2016a, GNU Octave 4.0
%
% See also:
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Laboratory
% Website: http://lasig.epfl.ch/
% Last revision: May 19, 2016
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund (WHFF, OFEV), project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'gps_week_time', @isnumeric);
addRequired(arg, 'point_source_id', @isnumeric);
addRequired(arg, 'gps_week_number', @(x) (size(x,2) == 1) && isnumeric(x));
addRequired(arg, 'flightline_id', @(x) (size(x,2) == 1) && isnumeric(x));
addOptional(arg, 'format', 'satellite', @(x) ismember(x, {'satellite','adjusted','utc'}));

parse(arg, gps_week_time, point_source_id, gps_week_number, flightline_id, format);


%% reformat time

% gps_week_number
[~, idxn_flid] = ismember(point_source_id, flightline_id);
gps_week_number_target = gps_week_number(idxn_flid);

% update LAS GPS time
switch format
    
    case 'satellite' % GPS satellite time
        
        % GPS satellite time (i.e. number of seconds since 1980-1-6)
        time = gps_week_number_target * 7 * 24 * 3600 + gps_week_time;
        
    case 'adjusted' % GPS adjusted time
        
        % GPS adjusted time (i.e. GPS satellite time - 10^9, the offset moves the time back to
        % near zero to improve floating point resolution.)
        time = gps_week_number_target * 7 * 24 * 3600 - 10^9 + gps_week_time;
        
    case 'utc' % Universal Time Count (UTC)
        
        % GPS, Global Positioning System time, is the atomic time scale implemented by the atomic clocks in the GPS ground 
        % control stations and the GPS satellites themselves. GPS time was zero at 0h 6-Jan-1980 and since it is not 
        % perturbed by leap seconds GPS is now ahead of UTC by 17 seconds.
        datum = datenum(1980,1,6,0,0,0); % * 24 * 3600  % start of GPS time number of days from January 0, 0000.
        week = datum + gps_week_number_target * 7; % * 24 * 3600; % add weeks
        time = week + gps_week_time / (24*3600); % add seconds since beginning of week
        
        % https://en.wikipedia.org/wiki/Leap_second
        % https://confluence.qps.nl/display/KBE/UTC+to+GPS+Time+Correction
        % https://www.andrews.edu/~tzs/timeconv/timealgorithm.html
        
        leap_dates = [datenum(1980,1,1,0,0,0),...
            datenum(1981,7,1,0,0,0),...
            datenum(1982,7,1,0,0,0),...
            datenum(1983,7,1,0,0,0),...
            datenum(1985,7,1,0,0,0),...
            datenum(1988,1,1,0,0,0),...
            datenum(1990,1,1,0,0,0),...
            datenum(1991,1,1,0,0,0),...
            datenum(1992,7,1,0,0,0),...
            datenum(1993,7,1,0,0,0),...
            datenum(1994,7,1,0,0,0),...
            datenum(1996,1,1,0,0,0),...
            datenum(1997,7,1,0,0,0),...
            datenum(1999,1,1,0,0,0),...
            datenum(2006,1,1,0,0,0),...
            datenum(2009,1,1,0,0,0),...
            datenum(2012,7,1,0,0,0),...
            datenum(2015,7,1,0,0,0),...
            ];
       
        leap_seconds = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17]';
        
        % add leap seconds
        [~, ~, idxn_bin] = histcounts(time, [leap_dates, inf]);
        time = time + leap_seconds(idxn_bin) / (24*3600);
        
end

