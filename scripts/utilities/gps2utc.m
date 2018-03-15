function utc_time = gps2utc(gps_time, varargin)
% GPS2UTC - Convert GPS satellite time (i.e. number of seconds since 1980-1-6) to UTC time.
% UTC_TIME = GPS2UTC(GPS_TIME, ...)
% 
% GPS, Global Positioning System time, is the atomic time scale implemented by the atomic clocks in the GPS ground 
% control stations and the GPS satellites themselves. GPS time was zero at 0h 6-Jan-1980 and since it is not 
% perturbed by leap seconds GPS time is now ahead of UTC time by several seconds.
%
% For more information, check:
% U.S. Naval Observatory (USNO), Time Service - GPS week, http://tycho.usno.navy.mil/gps_week.html
% NOAA, Continuously Operating Reference Station (CORS) - GPS calendar, http://www.ngs.noaa.gov/CORS/Gpscal.shtml
% QPS, UTC to GPS Time Correction, https://confluence.qps.nl/display/KBE/UTC+to+GPS+Time+Correction
% The Internet Engineering Task Force (IETF), Leap second list, https://www.ietf.org/timezones/data/leap-seconds.list
%
% Syntax:  utc_time = gps2utc(gps_time, ...)
%
% Inputs:
%    gps_time - Nx1 numeric vector, GPS satellite time (non-adjusted)
%
%    verbose (optional, default: true) - boolean value, verbosiy switch
%
% Outputs:
%    utc_time - Nx1 numeric vector, UTC time (serial date number)
%
% Example:
%
%    pc = LASread(zh_2014_coniferous.las);
%
%    utc_time = gps2utc(pc.record.gps_time + 10^9, ...
%                'verbose', true);
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

addRequired(arg, 'gps_time', @(x) (size(x,2) == 1) && isnumeric(x));
addParameter(arg, 'verbose', true, @(x) islogical(x) && (numel(x) == 1));

parse(arg, gps_time, varargin{:});


%% convert gps time to UTC time
% For update, see: https://www.ietf.org/timezones/data/leap-seconds.list

% define leap seconds since January 1, 1980
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
    datenum(2017,1,1,0,0,0),...
    ];

leap_seconds = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18]';


%% remove leap seconds

t = datenum(1980,1,6,0,0,0) + (gps_time / (24 * 3600));
[~, idxn_bin] = histc(t, leap_dates);
utc_time = t - (leap_seconds(idxn_bin) / (24 * 3600));


%% convert to vector and string format 

if arg.Results.verbose
    
    [Y,M,D,~,~,~] = datevec(utc_time);
    
    unique_date = sort(unique(datenum(Y,M,D)));
    
    % print unique days
    for j = 1:length(unique_date)
        
        fprintf('day %u: %s\n', j, datestr(unique_date(j)));
        
    end
    
end
