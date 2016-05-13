function trajectory = SBETread(filepath, verbose)
%SBETREAD - Reads a Smoothed Best Estimated Trajectory (SBET) Applanix file.
% TRAJECTORY = SBETread(FILEPATH, VERBOSE) reads the SBET file specified in FILEPATH
% into the structure TRAJECTORY.
%
% The SBET format contains the 8-byte double-precision following fields (in order).
%
% 1: GPS time (s)
% 2: Latitude (rad)
% 3: Longitude (rad)
% 4: Altitude (m)
% 5: X (E-W) velocity (m/s)
% 6: Y (N-S) velocity (m/s)
% 7: Z (Vertical) velocity (m/s)
% 8: Attitude roll component (rad)
% 9: Attitude pitch component (rad)
% 10: Attitude heading component (rad)
% 11: Wander angle (rad)
% 12: E-W acceleration (m/s^2)
% 13: N-S acceleration (m/s^2)
% 14: Vertical acceleration (m/s^2)
% 15: X-axis angular rate (rad/s)
% 16: Y-axis angular rate (rad/s)
% 17: Z-axis angular rate (rad/s)
%
% For more information see:
% http://vislab-ccom.unh.edu/~schwehr/Classes/2011/esci895-researchtools/21-python-binary-files.html
%
% Syntax:  trajectory = SBETread(filepath, verbose)
%
% Inputs:
%    filepath - The path to the input SBET file
%    verbose - If set to true, will display status information in the console
%
% Outputs:
%    trajectory - A structure containing the trajectory records
%
% Example:
%    trajectory = SBETread('E:\trajectories\sbet_2014041001.out', true);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2013b, GNU Octave 4.0.0 (configured for "i686-w64-mingw32")
%
% See also:
%
% Author: Matthew Parkan, EPFL - GIS Laboratory
% Website: http://lasig.epfl.ch/
% Last revision: September 4, 2015
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund, WHFF (OFEV) - project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details

%% check argument validity

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; % determine if system is Matlab or GNU Octave

if ~isOctave
    
    arg = inputParser;
    
    addRequired(arg, 'filepath', @checkFilepath);
    addRequired(arg, 'verbose', @islogical);
    
    parse(arg, filepath, verbose);
    
end

fclose('all');

    function checkFilepath(filepath)
        
        fid = fopen(filepath, 'r');
        
        if fid == -1
            
            error('Could not open file');
            
        end
        
        fclose(fid);
        
    end

%% record block format definition

% GPS time (s), 8 bytes
k = 1;
r.record(k).full_name = 'GPS time (s)';
r.record(k).short_name = 'gps_time';
r.record(k).type = 'double';
r.record(k).byte_length = 8;
r.record(k).n_values = 1;
r.record(k).print_format = '%g';
r.record(k).value = [];
k = k + 1;

% Latitude (rad), 8 bytes
r.record(k).full_name = 'Latitude (rad)';
r.record(k).short_name = 'latitude';
r.record(k).type = 'double';
r.record(k).byte_length = 8;
r.record(k).n_values = 1;
r.record(k).print_format = '%g';
r.record(k).value = [];
k = k + 1;

% Longitude (rad), 8 bytes
r.record(k).full_name = 'Longitude (rad)';
r.record(k).short_name = 'longitude';
r.record(k).type = 'double';
r.record(k).byte_length = 8;
r.record(k).n_values = 1;
r.record(k).print_format = '%g';
r.record(k).value = [];
k = k + 1;

% Altitude (m), 8 bytes
r.record(k).full_name = 'Altitude (m)';
r.record(k).short_name = 'altitude';
r.record(k).type = 'double';
r.record(k).byte_length = 8;
r.record(k).n_values = 1;
r.record(k).print_format = '%g';
r.record(k).value = [];
k = k + 1;

% E-W velocity (m/s), 8 bytes
r.record(k).full_name = 'X (E-W) velocity (m/s)';
r.record(k).short_name = 'x_velocity';
r.record(k).type = 'double';
r.record(k).byte_length = 8;
r.record(k).n_values = 1;
r.record(k).print_format = '%g';
r.record(k).value = [];
k = k + 1;

% N-S velocity (m/s), 8 bytes
r.record(k).full_name = 'Y (N-S) velocity (m/s)';
r.record(k).short_name = 'y_velocity';
r.record(k).type = 'double';
r.record(k).byte_length = 8;
r.record(k).n_values = 1;
r.record(k).print_format = '%g';
r.record(k).value = [];
k = k + 1;

% Vertical velocity (m/s), 8 bytes
r.record(k).full_name = 'Z (vertical) velocity (m/s)';
r.record(k).short_name = 'z_velocity';
r.record(k).type = 'double';
r.record(k).byte_length = 8;
r.record(k).n_values = 1;
r.record(k).print_format = '%g';
r.record(k).value = [];
k = k + 1;

% Attitude roll component (rad), 8 bytes
r.record(k).full_name = 'Attitude roll component (rad)';
r.record(k).short_name = 'roll';
r.record(k).type = 'double';
r.record(k).byte_length = 8;
r.record(k).n_values = 1;
r.record(k).print_format = '%g';
r.record(k).value = [];
k = k + 1;

% Attitude pitch component (rad), 8 bytes
r.record(k).full_name = 'Attitude pitch component (rad)';
r.record(k).short_name = 'pitch';
r.record(k).type = 'double';
r.record(k).byte_length = 8;
r.record(k).n_values = 1;
r.record(k).print_format = '%g';
r.record(k).value = [];
k = k + 1;

% Attitude heading component (rad), 8 bytes
r.record(k).full_name = 'Attitude heading component (rad)';
r.record(k).short_name = 'heading';
r.record(k).type = 'double';
r.record(k).byte_length = 8;
r.record(k).n_values = 1;
r.record(k).print_format = '%g';
r.record(k).value = [];
k = k + 1;

% Wander angle (rad), 8 bytes
r.record(k).full_name = 'Wander angle (rad)';
r.record(k).short_name = 'wander';
r.record(k).type = 'double';
r.record(k).byte_length = 8;
r.record(k).n_values = 1;
r.record(k).print_format = '%g';
r.record(k).value = [];
k = k + 1;

% E-W acceleration (m/s^2), 8 bytes
r.record(k).full_name = 'X (E-W) acceleration (m/s^2)';
r.record(k).short_name = 'x_acceleration';
r.record(k).type = 'double';
r.record(k).byte_length = 8;
r.record(k).n_values = 1;
r.record(k).print_format = '%g';
r.record(k).value = [];
k = k + 1;

% N-S acceleration (m/s^2), 8 bytes
r.record(k).full_name = 'N-S acceleration (m/s^2)';
r.record(k).short_name = 'y_acceleration';
r.record(k).type = 'double';
r.record(k).byte_length = 8;
r.record(k).n_values = 1;
r.record(k).print_format = '%g';
r.record(k).value = [];
k = k + 1;

% Vertical acceleration (m/s^2), 8 bytes
r.record(k).full_name = 'Z (vertical) acceleration (m/s^2)';
r.record(k).short_name = 'z_acceleration';
r.record(k).type = 'double';
r.record(k).byte_length = 8;
r.record(k).n_values = 1;
r.record(k).print_format = '%g';
r.record(k).value = [];
k = k + 1;

% X-axis angular rate (rad/s), 8 bytes
r.record(k).full_name = 'X-axis angular rate (rad/s)';
r.record(k).short_name = 'x_angular_rate';
r.record(k).type = 'double';
r.record(k).byte_length = 8;
r.record(k).n_values = 1;
r.record(k).print_format = '%g';
r.record(k).value = [];
k = k + 1;

% Y-axis angular rate (rad/s), 8 bytes
r.record(k).full_name = 'Y-axis angular rate (rad/s)';
r.record(k).short_name = 'y_angular_rate';
r.record(k).type = 'double';
r.record(k).byte_length = 8;
r.record(k).n_values = 1;
r.record(k).print_format = '%g';
r.record(k).value = [];
k = k + 1;

% Z-axis angular rate (rad/s), 8 bytes
r.record(k).full_name = 'Z-axis angular rate (rad/s)';
r.record(k).short_name = 'z_angular_rate';
r.record(k).type = 'double';
r.record(k).byte_length = 8;
r.record(k).n_values = 1;
r.record(k).print_format = '%g';
r.record(k).value = [];

%% open file

fid = fopen(filepath);

if fid == -1
    error('Could not open file')
end


%% read trajectory position record


if verbose
    
    fprintf('reading trajectory position record...');
    
end

% read blob
blob = fread(fid, '*double');

fclose(fid); % close file

% parse blob

step_size = length(r.record);

for j = 1:step_size
    
    r.record(j).value = blob(j:step_size:end); % cast uint32/uint64 types to single/double precision
    
end

%% format data for export

% record block

for j = 1:length(r.record)
    
    trajectory.(r.record(j).short_name) = r.record(j).value;
    
end


if verbose
    
    fprintf('done!\n');
    
end


end
