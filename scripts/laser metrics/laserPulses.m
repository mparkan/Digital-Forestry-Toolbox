function [pulse_number, pulse_returns, ind_pulse] = laserPulses(return_number)
% LASERPULSES - determines the pulse number (index) associated with the individual
% return numbers (sorted by acquisition GPS time).
%
% [PULSE_NUMBER, PULSE_RETURNS, IND_PULSE] = LASERPULSES(RETURN_NUMBER)
% returns the PULSE_NUMBER, number of return (PULSE_RETURNS) and pulse index (IND_PULSE) 
% associated with the individual return number specified in RETURN_NUMBER.
%
% Syntax:  [pulse_number, pulse_returns, ind_pulse] = laserPulses(return_number)
%
% Inputs:
%    return_number - nx1 uint8 vector containing the individual return
%    numbers (sorted by acquisition GPS time)
%
% Outputs:
%    pulse_number - nx1 float vector containing the laser pulse number
%    associated with each individual return
%
%    pulse_returns - nx1 uint8 containing the number of returns per pulse
%
%    ind_pulse - float vector containing the pulse index when
%    considering a pulse matrix with dimensions NxM:
%    
%            | return 1 | return 2 | return 3 | ... | return M |
%            ---------------------------------------------------
%    pulse 1 | return 1 | return 2 | return 3 | ... | return M |  
%    pulse 2 | return 1 | return 2 | return 3 | ... | return M |
%    pulse 3 | return 1 | return 2 | return 3 | ... | return M |
%    ...     |   ....   |   ....   |   ....   | ... |   ....   |
%    pulse N | return 1 | return 2 | return 3 | ... | return M |
%
%
% Example:
%
% pc = LASread('E:\als\mypoints.las', false, true);
% [pulse_number, pulse_returns, ind_pulse] = laserPulses(pc.record.return_number);
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

addRequired(arg, 'return_number', @(x) (size(x,2) == 1) && isnumeric(x));

parse(arg, return_number);


%% compute pulse number

rdiff = [true; logical(diff(return_number) <= 0)];
pulse_number = cumsum(rdiff);


%% compute number of return in pulse

pulse_returns = uint8(diff(find([rdiff; true])));


%% build pulse/index index matrix

n_pulses = pulse_number(end);
max_return_number = double(max(unique(return_number)));

pulse_return_number = reshape((repmat(pulse_returns, 1, max_return_number) - repmat(uint8(max_return_number-1:-1:0), n_pulses, 1))', [], 1);
pulse_return_number = pulse_return_number(pulse_return_number > 0);

ind_pulse = sub2ind([n_pulses max_return_number], pulse_number, double(pulse_return_number));
