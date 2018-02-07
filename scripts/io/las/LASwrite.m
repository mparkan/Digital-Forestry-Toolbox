function varargout = LASwrite(s, filepath, varargin)
%LASWRITE - Writes to the ASPRS LAS point cloud format (LAS versions 1.0 to 1.4 are supported)
% LASWRITE(S, FILEPATH, ...) writes the point cloud stored in the structure S to
% the LAS file specified in FILEPATH.
%
% The las format is a binary format defined by the American Society for
% Photogrammetry and Remote Sensing (ASPRS) to store point clouds. The
% specification has gone through several revisions since its creation:
%
% LAS 1.0, 09/05/2003
% LAS 1.1, 07/05/2005
% LAS 1.2, 02/09/2008
% LAS 1.3, 24/08/2010
% LAS 1.4, 15/07/2013
%
% For more information see: ASPRS, LASer (LAS) File Format Exchange
% Activities,
% http://www.asprs.org/Committee-General/LASer-LAS-File-Format-Exchange-Activities.html
%
% Syntax:  LASwrite(s, filepath, ...)
%
% Inputs:
%    s - A structure containing the point cloud
%    filepath - The path to the input LAS file
%    version (optional, default: 14) - The output LAS version (see LAS specification)
%    systemID (optional, default: 'OTHER') - The output LAS version (see LAS specification)
%    recordFormat (optional, default: 3) - The output LAS record format (see LAS specification)
%    verbose (optional, default: true) - If set to true, will display information
%              in the console
%
% Example:
%    s = LASwrite(s, 'E:\las_data\my_points.las', 'version', 12, 'verbose', false)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2017b, GNU Octave 4.2.1 (configured for "x86_64-w64-mingw32")
%
% See also: LASREAD
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory (LASIG)
% Website: http://mparkan.github.io/Digital-Forestry-Toolbox/
% Last revision: February 7, 2018
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund, WHFF (OFEV) - project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details

fclose('all'); % close any open files

%% setup constants

MACHINE_FORMAT = 'ieee-le'; % all data is in little-endian format
OCTAVE_FLAG = (exist('OCTAVE_VERSION', 'builtin') ~= 0); % determine if system is Matlab or GNU Octave
CURRENT_DATE = clock; % current date
AUTH_RECORD_FORMAT_ID = 0:10; % authorized point record format IDs
AUTH_LAS_MAJOR_VERSIONS = [1]; % authorized LAS major versions
AUTH_LAS_MINOR_VERSIONS = [0, 1, 2, 3, 4]; % authorized LAS minor versions
AUTH_LAS_VERSIONS = [10, 11, 12, 13, 14]; % authorized LAS versions
PUBLIC_HEADER_SIZES = [227, 227, 227, 235, 375]; % public header sizes for LAS 1.0-1.4


%% check argument validity

arg = inputParser;

addRequired(arg, 's', @isstruct);
addOptional(arg, 'filepath', [], @ischar);
addParameter(arg, 'version', 14, @(x) (ismember(x, AUTH_LAS_VERSIONS) && numel(x) == 1));
addParameter(arg, 'systemID', 'OTHER', @(x) ischar(x) && (length(x) <= 32));
addParameter(arg, 'recordFormat', 3, @(x) (ismember(x, AUTH_RECORD_FORMAT_ID) && numel(x) == 1));
addParameter(arg, 'verbose', true, @(x) islogical(x) && (numel(x) == 1));

parse(arg, s, filepath, varargin{:});


%% determine LAS version

las_version = fix(arg.Results.version / 10)*10 + mod(arg.Results.version, 10);
[~, ind_las_version] = ismember(las_version, AUTH_LAS_VERSIONS);


%% compute spatial extent

x_min = min(s.record.x);
x_max = max(s.record.x);
y_min = min(s.record.y);
y_max = max(s.record.y);
z_min = min(s.record.z);
z_max = max(s.record.z);


%% public header block format definition

k = 1;

% File Signature (LASF), char[4], 4 bytes, *
r.header(k).compatibility = [10, 11, 12, 13, 14];
r.header(k).full_name = 'File Signature (LASF)';
r.header(k).short_name = 'file_signature';
r.header(k).type = {'char', 'char', 'char', 'char', 'char'};
r.header(k).storage_type = {'char', 'char', 'char', 'char', 'char'};
r.header(k).byte_length = [1, 1, 1, 1, 1]; % byte length per value
r.header(k).n_values = [4, 4, 4, 4, 4];
r.header(k).flag_bit_field = [false, false, false, false, false];
r.header(k).print_format = {'%s', '%s', '%s', '%s', '%s'};
r.header(k).default_value = 'LASF';
r.header(k).value = r.header(k).default_value;
r.header(k).validation = []; % @(x) strcmp(x, 'LASF');
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = [];
k = k + 1;

% Reserved (LAS 1.0 only), unsigned long, 4 bytes
r.header(k).compatibility = [10];
r.header(k).full_name = 'Reserved (LAS 1.0 only)';
r.header(k).short_name = 'reserved';
r.header(k).type = {'uint32'};
r.header(k).storage_type = {'double'};
r.header(k).byte_length = [4]; % byte length per value
r.header(k).n_values = [1];
r.header(k).flag_bit_field = [false];
r.header(k).print_format = {'%u'};
r.header(k).default_value = 0;
r.header(k).value = r.header(k).default_value;
r.header(k).validation = []; % @(x) (x == 0) && (numel(x) == 1);
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = [];
k = k + 1;

% File Source ID, unsigned short, 2 bytes, *
r.header(k).compatibility = [11, 12, 13, 14];
r.header(k).full_name = 'File Source ID';
r.header(k).short_name = 'file_source_id';
r.header(k).type = {'uint16', 'uint16', 'uint16', 'uint16'};
r.header(k).storage_type = {'double', 'double', 'double', 'double'};
r.header(k).byte_length = [2, 2, 2, 2];
r.header(k).n_values = [1, 1, 1, 1];
r.header(k).flag_bit_field = [false, false, false, false];
r.header(k).print_format = {'%u', '%u', '%u', '%u'};
r.header(k).default_value = 0;
r.header(k).value = [];
r.header(k).validation = @(x) (x >= 0) && (x <= intmax('uint16')) && (numel(x) == 1);
r.header(k).error_id = 'LASwrite:fileSourceID';
r.header(k).error_message = 'The file source ID field was not specified';
r.header(k).warning_message = 'The file source ID field was not specified';
k = k + 1;

% Reserved (LAS 1.1 only), unsigned short, 2 bytes
r.header(k).compatibility = [11];
r.header(k).full_name = 'Reserved (LAS 1.1 only)';
r.header(k).short_name = 'reserved';
r.header(k).type = {'uint16'};
r.header(k).storage_type = {'double'};
r.header(k).byte_length = [2]; % byte length per value
r.header(k).n_values = [1];
r.header(k).flag_bit_field = [false];
r.header(k).print_format = {'%u'};
r.header(k).default_value = 0;
r.header(k).value = r.header(k).default_value;
r.header(k).validation = []; % @(x) (x == 0) && (numel(x) == 1);
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = 'The Reserved field (LAS 1.1 only) was not specified';
k = k + 1;

% Global Encoding - GPS Time Type, bit 0, 1 bit, *
% The meaning of GPS Time in the Point Records
% 0 (not set) -> GPS time in the point record fields is GPS Week Time (the same as previous versions of LAS)
% 1 (set) -> GPS Time is standard GPS Time (satellite GPS Time) minus 1 x 109.
% The offset moves the time back to near zero to improve floating point resolution.
r.header(k).compatibility = [12, 13, 14];
r.header(k).full_name = 'Global Encoding - GPS Time Type';
r.header(k).short_name = 'global_encoding_gps_time_type';
r.header(k).type = {'ubit1', 'ubit1', 'ubit1'};
r.header(k).storage_type = {'logical', 'logical', 'logical'};
r.header(k).parent_type = {'uint16', 'uint16', 'uint16'};
r.header(k).byte_length = [1/8, 1/8, 1/8];
r.header(k).bit_position = [0, 0, 0];
r.header(k).n_values = [1, 1, 1];
r.header(k).flag_bit_field = [true, true, true];
r.header(k).print_format = {'%u', '%u', '%u'};
r.header(k).default_value = false;
r.header(k).value = [];
r.header(k).validation = @(x) islogical(x) && (numel(x) == 1);
r.header(k).error_id = 'LASwrite:gpsTimeType';
r.header(k).error_message = 'The global encoding GPS Time Type (LAS >= 1.2 only) field was not specified';
r.header(k).warning_message = [];
k = k + 1;

% Global Encoding (LAS >= 1.3 only) - Waveform Data Packets Internal, bit 1, 1 bit, *
r.header(k).compatibility = [13, 14];
r.header(k).full_name = 'Global Encoding - Waveform Data Packets Internal (LAS >= 1.3 only)';
r.header(k).short_name = 'global_encoding_waveform_data_packets_internal';
r.header(k).type = {'ubit1', 'ubit1'};
r.header(k).parent_type = {'uint16', 'uint16'};
r.header(k).storage_type = {'logical', 'logical'};
r.header(k).byte_length = [1/8, 1/8];
r.header(k).bit_position = [1, 1];
r.header(k).n_values = [1, 1];
r.header(k).flag_bit_field = [true, true];
r.header(k).print_format = {'%u', '%u'};
r.header(k).default_value = false;
r.header(k).value = [];
r.header(k).validation = @(x) islogical(x) && (numel(x) == 1);
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = 'The Global Encoding - Waveform Data Packets Internal (LAS >= 1.3 only) field was not specified';
k = k + 1;

% Global Encoding (LAS >= 1.3 only) - Waveform Data Packets External, bit 2, 1 bit, *
r.header(k).compatibility = [13, 14];
r.header(k).full_name = 'Global Encoding - Waveform Data Packets External (LAS >= 1.3 only)';
r.header(k).short_name = 'global_encoding_waveform_data_packets_external';
r.header(k).type = {'ubit1', 'ubit1'};
r.header(k).parent_type = {'uint16', 'uint16'};
r.header(k).storage_type = {'logical', 'logical'};
r.header(k).byte_length = [1/8, 1/8];
r.header(k).bit_position = [2, 2];
r.header(k).n_values = [1, 1];
r.header(k).flag_bit_field = [true, true];
r.header(k).print_format = {'%u', '%u'};
r.header(k).default_value = false;
r.header(k).value = [];
r.header(k).validation = @(x) islogical(x) && (numel(x) == 1);
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = 'The Global Encoding - Waveform Data Packets External (LAS >= 1.3 only) field was not specified';
k = k + 1;

% Global Encoding (LAS >= 1.3 only) - Return numbers have been synthetically generated, bit 3, 1 bit, *
r.header(k).compatibility = [13, 14];
r.header(k).full_name = 'Global Encoding - Return numbers have been synthetically generated (LAS >= 1.3 only)';
r.header(k).short_name = 'global_encoding_synthetic_return_numbers';
r.header(k).type = {'ubit1', 'ubit1'};
r.header(k).parent_type = {'uint16', 'uint16'};
r.header(k).storage_type = {'logical', 'logical'};
r.header(k).byte_length = [1/8, 1/8];
r.header(k).bit_position = [3, 3];
r.header(k).n_values = [1, 1];
r.header(k).flag_bit_field = [true, true];
r.header(k).print_format = {'%u', '%u'};
r.header(k).default_value = false;
r.header(k).value = [];
r.header(k).validation = @(x) islogical(x) && (numel(x) == 1);
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = 'The Global Encoding - Return numbers have been synthetically generated (LAS >= 1.3 only) field was not specified';
k = k + 1;

% WKT (LAS >= 1.4 only) - If set, the Coordinate Reference System (CRS) is WKT, bit 4, 1 bit, *
r.header(k).compatibility = [14];
r.header(k).full_name = 'Global Encoding - WKT (LAS >= 1.4 only)';
r.header(k).short_name = 'global_encoding_wkt';
r.header(k).type = {'ubit1'};
r.header(k).parent_type = {'uint16'};
r.header(k).storage_type = {'logical'};
r.header(k).byte_length = [1/8];
r.header(k).bit_position = [4];
r.header(k).n_values = [1];
r.header(k).flag_bit_field = [true];
r.header(k).print_format = {'%u'};
r.header(k).default_value = false;
r.header(k).value = [];
r.header(k).validation = @(x) islogical(x) && (numel(x) == 1);
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = 'The Global Encoding - WKT (LAS >= 1.4 only) field was not specified';
k = k + 1;

% Global Encoding - Reserved (LAS >= 1.2 only), *
r.header(k).compatibility = [12, 13, 14];
r.header(k).full_name = 'Global Encoding - Reserved (LAS >= 1.2 only)';
r.header(k).short_name = 'global_encoding_reserved';
r.header(k).type = {'ubit15', 'ubit12', 'ubit11'};
r.header(k).parent_type = {'uint16', 'uint16', 'uint16'};
r.header(k).storage_type = {'double', 'double', 'double'};
r.header(k).byte_length = [15/8, 12/8, 11/8];
r.header(k).bit_position = [1, 4, 5];
r.header(k).n_values = [1, 1, 1];
r.header(k).flag_bit_field = [true, true, true];
r.header(k).print_format = {'%u', '%u', '%u'};
r.header(k).default_value = 0;
r.header(k).value = r.header(k).default_value;
r.header(k).validation = [];
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = 'The Global Encoding - Reserved (LAS >= 1.2 only) field was not specified';
k = k + 1;

% Project ID - GUID data 1, unsigned long, 4 bytes
r.header(k).compatibility = [10, 11, 12, 13, 14];
r.header(k).full_name = 'Project ID - GUID data 1';
r.header(k).short_name = 'project_id_1';
r.header(k).type = {'uint32', 'uint32', 'uint32', 'uint32', 'uint32'};
r.header(k).storage_type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).byte_length = [4, 4, 4, 4, 4];
r.header(k).n_values = [1, 1 ,1 ,1 ,1];
r.header(k).flag_bit_field = [false, false, false, false, false];
r.header(k).print_format = {'%u', '%u', '%u', '%u', '%u'};
r.header(k).default_value = 0;
r.header(k).value = [];
r.header(k).validation = @(x) isnumeric(x) && (numel(x) == 1);
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = 'The Project ID - GUID data 1 field was not specified';
k = k + 1;

% Project ID - GUID data 2, unsigned short, 2 bytes
r.header(k).compatibility = [10, 11, 12, 13, 14];
r.header(k).full_name = 'Project ID - GUID data 2';
r.header(k).short_name = 'project_id_2';
r.header(k).type = {'uint16', 'uint16', 'uint16', 'uint16', 'uint16'};
r.header(k).storage_type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).byte_length = [2, 2, 2, 2, 2];
r.header(k).n_values = [1, 1 ,1 ,1 ,1];
r.header(k).flag_bit_field = [false, false, false, false, false];
r.header(k).print_format = {'%u', '%u', '%u', '%u', '%u'};
r.header(k).default_value = 0;
r.header(k).value = [];
r.header(k).validation = @(x) isnumeric(x) && (numel(x) == 1);
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = 'The Project ID - GUID data 2 field was not specified';
k = k + 1;

% Project ID - GUID data 3, unsigned short, 2 bytes
r.header(k).compatibility = [10, 11, 12, 13, 14];
r.header(k).full_name = 'Project ID - GUID data 3';
r.header(k).short_name = 'project_id_3';
r.header(k).type = {'uint16', 'uint16', 'uint16', 'uint16', 'uint16'};
r.header(k).storage_type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).byte_length = [2, 2, 2, 2, 2];
r.header(k).n_values = [1, 1 ,1 ,1 ,1];
r.header(k).flag_bit_field = [false, false, false, false, false];
r.header(k).print_format = {'%u', '%u', '%u', '%u', '%u'};
r.header(k).default_value = 0;
r.header(k).value = [];
r.header(k).validation = @(x) isnumeric(x) && (numel(x) == 1);
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = 'The Project ID - GUID data 3 field was not specified';
k = k + 1;

% Project ID - GUID data 4, unsigned char[8], 8 bytes
r.header(k).compatibility = [10, 11, 12, 13, 14];
r.header(k).full_name = 'Project ID - GUID data 4';
r.header(k).short_name = 'project_id_4';
r.header(k).type = {'uint64', 'uint64', 'uint64', 'uint64', 'uint64'};
r.header(k).storage_type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).byte_length = [8, 8, 8, 8 ,8];
r.header(k).n_values = [1, 1 ,1 ,1 ,1];
r.header(k).flag_bit_field = [false, false, false, false, false];
r.header(k).print_format = {'%u', '%u', '%u', '%u', '%u'};
r.header(k).default_value = 0;
r.header(k).value = [];
r.header(k).validation = @(x) isnumeric(x) && (numel(x) == 1);
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = 'The Project ID - GUID data 4 field was not specified';
k = k + 1;

% Version Major, unsigned char, 1 byte, *
r.header(k).compatibility = [10, 11, 12, 13, 14];
r.header(k).full_name = 'Version Major';
r.header(k).short_name = 'version_major';
r.header(k).type = {'uint8', 'uint8', 'uint8', 'uint8', 'uint8'};
r.header(k).storage_type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).byte_length = [1, 1, 1, 1, 1];
r.header(k).n_values = [1, 1 ,1 ,1 ,1];
r.header(k).flag_bit_field = [false, false, false, false, false];
r.header(k).print_format = {'%u', '%u', '%u', '%u', '%u'};
r.header(k).default_value = fix(arg.Results.version / 10);
r.header(k).value = r.header(k).default_value;
r.header(k).validation = []; %@(x) ismember(x, AUTH_LAS_MAJOR_VERSIONS) && (numel(x) == 1);
r.header(k).error_id = 'LASwrite:versionMajor';
r.header(k).error_message = 'The major LAS version field was not specified or is invalid';
r.header(k).warning_message = [];
k = k + 1;

% Version Minor, unsigned char, 1 byte, *
r.header(k).compatibility = [10, 11, 12, 13, 14];
r.header(k).full_name = 'Version Minor';
r.header(k).short_name = 'version_minor';
r.header(k).type = {'uint8', 'uint8', 'uint8', 'uint8', 'uint8'};
r.header(k).storage_type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).byte_length = [1, 1, 1, 1, 1];
r.header(k).n_values = [1, 1 ,1 ,1 ,1];
r.header(k).flag_bit_field = [false, false, false, false, false];
r.header(k).print_format = {'%u', '%u', '%u', '%u', '%u'};
r.header(k).default_value = mod(arg.Results.version, 10);
r.header(k).value = r.header(k).default_value;
r.header(k).validation = []; %@(x) ismember(x, AUTH_LAS_MINOR_VERSIONS) && (numel(x) == 1);
r.header(k).error_id = 'LASwrite:minorVersion';
r.header(k).error_message = 'The minor LAS version field was not specified or is invalid';
r.header(k).warning_message = [];
k = k + 1;

% System Identifier, char[32], 32 bytes, *
r.header(k).compatibility = [10, 11, 12, 13, 14];
r.header(k).full_name = 'System Identifier';
r.header(k).short_name = 'system_identifier';
r.header(k).type = {'char', 'char', 'char', 'char', 'char'};
r.header(k).storage_type = {'char', 'char', 'char', 'char', 'char'};
r.header(k).byte_length = [1, 1, 1, 1, 1];
r.header(k).n_values = [32, 32, 32, 32, 32];
r.header(k).flag_bit_field = [false, false, false, false, false];
r.header(k).print_format = {'%s', '%s', '%s', '%s', '%s'};
r.header(k).default_value = arg.Results.systemID;
r.header(k).value = r.header(k).default_value;
r.header(k).validation = []; % @(x) ischar(x) && (length(x) <= 32);
r.header(k).error_id = 'LASwrite:systemIdentifier';
r.header(k).error_message = 'The System identifier field was not specified';
r.header(k).warning_message = [];
k = k + 1;

% Generating Software, char[32], 32 bytes, *
r.header(k).compatibility = [10, 11, 12, 13, 14];
r.header(k).full_name = 'Generating Software';
r.header(k).short_name = 'generating_software';
r.header(k).type = {'char', 'char', 'char', 'char', 'char'};
r.header(k).storage_type = {'char', 'char', 'char', 'char', 'char'};
r.header(k).byte_length = [1, 1, 1, 1, 1];
r.header(k).n_values = [32, 32, 32, 32, 32];
r.header(k).flag_bit_field = [false, false, false, false, false];
r.header(k).print_format = {'%s', '%s', '%s', '%s', '%s'};
r.header(k).default_value = 'Matlab / GNU Octave - LASwrite';
r.header(k).value = r.header(k).default_value;
r.header(k).validation = []; % @(x) ischar(x) && (length(x) <= 32);
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = [];
k = k + 1;

% File Creation Day of Year, unsigned short, 2 bytes
r.header(k).compatibility = [10, 11, 12, 13, 14];
r.header(k).full_name = 'File Creation Day of Year';
r.header(k).short_name = 'file_creation_doy';
r.header(k).type = {'uint16', 'uint16', 'uint16', 'uint16', 'uint16'};
r.header(k).storage_type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).byte_length = [2, 2, 2, 2, 2];
r.header(k).n_values = [1, 1, 1, 1, 1];
r.header(k).flag_bit_field = [false, false, false, false, false];
r.header(k).print_format = {'%u', '%u', '%u', '%u', '%u'};
r.header(k).default_value = floor(now-datenum(CURRENT_DATE(1),1,1) + 1);
r.header(k).value = r.header(k).default_value;
r.header(k).validation = []; %@(x) isa(x, 'uint16') && (numel(x) == 1) && (x <= 366);
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = [];
k = k + 1;

% File Creation Year, unsigned short, 2 bytes
r.header(k).compatibility = [10, 11, 12, 13, 14];
r.header(k).full_name = 'File Creation Year';
r.header(k).short_name = 'file_creation_year';
r.header(k).type = {'uint16', 'uint16', 'uint16', 'uint16', 'uint16'};
r.header(k).storage_type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).byte_length = [2, 2, 2, 2, 2];
r.header(k).n_values = [1, 1, 1, 1, 1];
r.header(k).flag_bit_field = [false, false, false, false, false];
r.header(k).print_format = {'%u', '%u', '%u', '%u', '%u'};
r.header(k).default_value = CURRENT_DATE(1);
r.header(k).value = r.header(k).default_value;
r.header(k).validation = []; %@(x) isa(x, 'uint16') && (numel(x) == 1) && (x <= 366);
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = [];
k = k + 1;

% Header Size, unsigned short, 2 bytes, *
r.header(k).compatibility = [10, 11, 12, 13, 14];
r.header(k).full_name = 'Header Size';
r.header(k).short_name = 'header_size';
r.header(k).type = {'uint16', 'uint16', 'uint16', 'uint16', 'uint16'};
r.header(k).storage_type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).byte_length = [2, 2, 2, 2, 2];
r.header(k).n_values = [1, 1, 1, 1, 1];
r.header(k).flag_bit_field = [false, false, false, false, false];
r.header(k).print_format = {'%u', '%u', '%u', '%u', '%u'};
r.header(k).default_value = PUBLIC_HEADER_SIZES(ind_las_version);
r.header(k).value = r.header(k).default_value;
r.header(k).validation = [];
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = [];
k = k + 1;

% Offset to point data, unsigned long, 4 bytes, *
r.header(k).compatibility = [10, 11, 12, 13, 14];
r.header(k).full_name = 'Offset to point data';
r.header(k).short_name = 'offset_to_data';
r.header(k).type = {'uint32', 'uint32', 'uint32', 'uint32', 'uint32'};
r.header(k).storage_type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).byte_length = [4, 4, 4, 4, 4];
r.header(k).n_values = [1, 1, 1, 1, 1];
r.header(k).flag_bit_field = [false, false, false, false, false];
r.header(k).print_format = {'%u', '%u', '%u', '%u', '%u'};
r.header(k).default_value = [];
r.header(k).value = [];
r.header(k).validation = @(x) (x >= intmin('uint32')) && (x <= intmax('uint32')) && (numel(x) == 1);
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = [];
k = k + 1;

% Number of Variable Length Records, unsigned long, 4 bytes, *
r.header(k).compatibility = [10, 11, 12, 13, 14];
r.header(k).full_name = 'Number of Variable Length Records';
r.header(k).short_name = 'n_variable_length_records';
r.header(k).type = {'uint32', 'uint32', 'uint32', 'uint32', 'uint32'};
r.header(k).storage_type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).byte_length = [4, 4, 4, 4, 4];
r.header(k).n_values = [1, 1, 1, 1, 1];
r.header(k).flag_bit_field = [false, false, false, false, false];
r.header(k).print_format = {'%u', '%u', '%u', '%u', '%u'};
r.header(k).default_value = [];
r.header(k).value = [];
r.header(k).validation = @(x) (x >= intmin('uint32')) && (x <= intmax('uint32')) && (numel(x) == 1);
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = 'WARNING: The mandatory GeoKeyDirectoryTag was not specified in the variable length header\n';
k = k + 1;

% Point Data Format ID (0-99 for spec), unsigned char, 1 byte, *
r.header(k).compatibility = [10, 11, 12, 13, 14];
r.header(k).full_name = 'Point Data Format ID (0-99 for spec)';
r.header(k).short_name = 'point_data_format_id';
r.header(k).type = {'uint8', 'uint8', 'uint8', 'uint8', 'uint8'};
r.header(k).storage_type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).byte_length = [1, 1, 1, 1, 1];
r.header(k).n_values = [1, 1, 1, 1, 1];
r.header(k).flag_bit_field = [false, false, false, false, false];
r.header(k).print_format = {'%u', '%u', '%u', '%u', '%u'};
r.header(k).default_value = uint8(arg.Results.recordFormat);
r.header(k).value = r.header(k).default_value;
r.header(k).validation = [];
r.header(k).error_id = 'LASwrite:formatID';
r.header(k).error_message = 'The Point Data Format ID field was not specified';
r.header(k).warning_message = [];
k = k + 1;

% Point Data Record Length, unsigned short, 2 bytes, *
r.header(k).compatibility = [10, 11, 12, 13, 14];
r.header(k).full_name = 'Point Data Record Length';
r.header(k).short_name = 'point_data_record_length';
r.header(k).type = {'uint16', 'uint16', 'uint16', 'uint16', 'uint16'};
r.header(k).storage_type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).byte_length = [2, 2, 2, 2, 2];
r.header(k).n_values = [1, 1, 1, 1, 1];
r.header(k).flag_bit_field = [false, false, false, false, false];
r.header(k).print_format = {'%u', '%u', '%u', '%u', '%u'};
r.header(k).default_value = [];
r.header(k).value = [];
r.header(k).validation = @(x) isa(x, 'uint16') && (numel(x) == 1);
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = [];
k = k + 1;

% Number of point records (Legacy for LAS >= 1.4), unsigned long, 4 bytes, *
r.header(k).compatibility = [10, 11, 12, 13, 14];
r.header(k).full_name = 'Number of point records (Legacy for LAS >= 1.4)';
r.header(k).short_name = 'n_point_records';
r.header(k).type = {'uint32', 'uint32', 'uint32', 'uint32', 'uint32'};
r.header(k).storage_type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).byte_length = [4, 4, 4, 4, 4];
r.header(k).n_values = [1, 1, 1, 1, 1];
r.header(k).flag_bit_field = [false, false, false, false, false];
r.header(k).print_format = {'%u', '%u', '%u', '%u', '%u'};
r.header(k).default_value = length(s.record.x);
r.header(k).value = r.header(k).default_value;
r.header(k).validation = [];
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = [];
k = k + 1;

% Number of points by return (Legacy for LAS >= 1.4), unsigned long[5], 20 bytes, *
r.header(k).compatibility = [10, 11, 12, 13, 14];
r.header(k).full_name = 'Number of points by return (1-5) (Legacy for LAS >= 1.4)';
r.header(k).short_name = 'n_points_by_return';
r.header(k).type = {'uint32', 'uint32', 'uint32', 'uint32', 'uint32'};
r.header(k).storage_type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).byte_length = [4, 4, 4, 4, 4];
r.header(k).n_values = [5, 5, 5, 5, 5];
r.header(k).flag_bit_field = [false, false, false, false, false];
r.header(k).print_format = {'%u', '%u', '%u', '%u', '%u'};
r.header(k).default_value = [];
r.header(k).value = [];
r.header(k).validation = @(x) (x >= intmin('uint32')) && (x <= intmax('uint32')) && (numel(x) == 5);
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = [];
k = k + 1;

% X scale factor, double, 8 bytes, *
r.header(k).compatibility = [10, 11, 12, 13, 14];
r.header(k).full_name = 'X scale factor';
r.header(k).short_name = 'x_scale_factor';
r.header(k).type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).storage_type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).byte_length = [8, 8, 8, 8, 8];
r.header(k).n_values = [1, 1, 1, 1, 1];
r.header(k).flag_bit_field = [false, false, false, false, false];
r.header(k).print_format = {'%g', '%g', '%g', '%g', '%g'};
r.header(k).default_value = 0.001;
r.header(k).value = [];
r.header(k).validation = @(x) isa(x, 'double') && (x <= 1) && (x > 0) && (numel(x) == 1);
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = [];
k = k + 1;

% Y scale factor, double, 8 bytes, *
r.header(k).compatibility = [10, 11, 12, 13, 14];
r.header(k).full_name = 'Y scale factor';
r.header(k).short_name = 'y_scale_factor';
r.header(k).type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).storage_type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).byte_length = [8, 8, 8, 8, 8];
r.header(k).n_values = [1, 1, 1, 1, 1];
r.header(k).flag_bit_field = [false, false, false, false, false];
r.header(k).print_format = {'%g', '%g', '%g', '%g', '%g'};
r.header(k).default_value = 0.001;
r.header(k).value = [];
r.header(k).validation = @(x) isa(x, 'double') && (x <= 1) && (x > 0) && (numel(x) == 1);
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = [];
k = k + 1;

% Z scale factor, double, 8 bytes, *
r.header(k).compatibility = [10, 11, 12, 13, 14];
r.header(k).full_name = 'Z scale factor';
r.header(k).short_name = 'z_scale_factor';
r.header(k).type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).storage_type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).byte_length = [8, 8, 8, 8, 8];
r.header(k).n_values = [1, 1, 1, 1, 1];
r.header(k).flag_bit_field = [false, false, false, false, false];
r.header(k).print_format = {'%g', '%g', '%g', '%g', '%g'};
r.header(k).default_value = 0.001;
r.header(k).value = [];
r.header(k).validation = @(x) isa(x, 'double') && (x <= 1) && (x > 0) && (numel(x) == 1);
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = [];
k = k + 1;

% X offset, double, 8 bytes, *
r.header(k).compatibility = [10, 11, 12, 13, 14];
r.header(k).full_name = 'X offset';
r.header(k).short_name = 'x_offset';
r.header(k).type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).storage_type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).byte_length = [8, 8, 8, 8, 8];
r.header(k).n_values = [1, 1, 1, 1, 1];
r.header(k).flag_bit_field = [false, false, false, false, false];
r.header(k).print_format = {'%3.3f', '%3.3f', '%3.3f', '%3.3f', '%3.3f'};
r.header(k).default_value = floor(x_min);
r.header(k).value = [];
r.header(k).validation = @(x) isa(x, 'double') && (numel(x) == 1);
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = [];
k = k + 1;

% Y offset, double, 8 bytes, *
r.header(k).compatibility = [10, 11, 12, 13, 14];
r.header(k).full_name = 'Y offset';
r.header(k).short_name = 'y_offset';
r.header(k).type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).storage_type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).byte_length = [8, 8, 8, 8, 8];
r.header(k).n_values = [1, 1, 1, 1, 1];
r.header(k).flag_bit_field = [false, false, false, false, false];
r.header(k).print_format = {'%3.3f', '%3.3f', '%3.3f', '%3.3f', '%3.3f'};
r.header(k).default_value = floor(y_min);
r.header(k).value = [];
r.header(k).validation = @(x) isa(x, 'double') && (numel(x) == 1);
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = [];
k = k + 1;

% Z offset, double, 8 bytes, *
r.header(k).compatibility = [10, 11, 12, 13, 14];
r.header(k).full_name = 'Z offset';
r.header(k).short_name = 'z_offset';
r.header(k).type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).storage_type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).byte_length = [8, 8, 8, 8, 8];
r.header(k).n_values = [1, 1, 1, 1, 1];
r.header(k).flag_bit_field = [false, false, false, false, false];
r.header(k).print_format = {'%3.3f', '%3.3f', '%3.3f', '%3.3f', '%3.3f'};
r.header(k).default_value = floor(z_min);
r.header(k).value = [];
r.header(k).validation = @(x) isa(x, 'double') && (numel(x) == 1);
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = [];
k = k + 1;

% Max X, double, 8 bytes, *
r.header(k).compatibility = [10, 11, 12, 13, 14];
r.header(k).full_name = 'Max X';
r.header(k).short_name = 'max_x';
r.header(k).type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).storage_type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).byte_length = [8, 8, 8, 8, 8];
r.header(k).n_values = [1, 1, 1, 1, 1];
r.header(k).flag_bit_field = [false, false, false, false, false];
r.header(k).print_format = {'%3.3f', '%3.3f', '%3.3f', '%3.3f', '%3.3f'};
r.header(k).default_value = x_max;
r.header(k).value = r.header(k).default_value;
r.header(k).validation = [];
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = [];
k = k + 1;

% Min X, double, 8 bytes, *
r.header(k).compatibility = [10, 11, 12, 13, 14];
r.header(k).full_name = 'Min X';
r.header(k).short_name = 'min_x';
r.header(k).type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).storage_type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).byte_length = [8, 8, 8, 8, 8];
r.header(k).n_values = [1, 1, 1, 1, 1];
r.header(k).flag_bit_field = [false, false, false, false, false];
r.header(k).print_format = {'%3.3f', '%3.3f', '%3.3f', '%3.3f', '%3.3f'};
r.header(k).default_value = x_min;
r.header(k).value = r.header(k).default_value;
r.header(k).validation = [];
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = [];
k = k + 1;

% Max Y, double, 8 bytes, *
r.header(k).compatibility = [10, 11, 12, 13, 14];
r.header(k).full_name = 'Max Y';
r.header(k).short_name = 'max_y';
r.header(k).type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).storage_type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).byte_length = [8, 8, 8, 8, 8];
r.header(k).n_values = [1, 1, 1, 1, 1];
r.header(k).flag_bit_field = [false, false, false, false, false];
r.header(k).print_format = {'%3.3f', '%3.3f', '%3.3f', '%3.3f', '%3.3f'};
r.header(k).default_value = y_max;
r.header(k).value = r.header(k).default_value;
r.header(k).validation = [];
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = [];
k = k + 1;

% Min Y, double, 8 bytes, *
r.header(k).compatibility = [10, 11, 12, 13, 14];
r.header(k).full_name = 'Min Y';
r.header(k).short_name = 'min_y';
r.header(k).type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).storage_type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).byte_length = [8, 8, 8, 8, 8];
r.header(k).n_values = [1, 1, 1, 1, 1];
r.header(k).flag_bit_field = [false, false, false, false, false];
r.header(k).print_format = {'%3.3f', '%3.3f', '%3.3f', '%3.3f', '%3.3f'};
r.header(k).default_value = y_min;
r.header(k).value = r.header(k).default_value;
r.header(k).validation = [];
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = [];
k = k + 1;

% Max Z, double, 8 bytes, *
r.header(k).compatibility = [10, 11, 12, 13, 14];
r.header(k).full_name = 'Max Z';
r.header(k).short_name = 'max_z';
r.header(k).type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).storage_type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).byte_length = [8, 8, 8, 8, 8];
r.header(k).n_values = [1, 1, 1, 1, 1];
r.header(k).flag_bit_field = [false, false, false, false, false];
r.header(k).print_format = {'%3.3f', '%3.3f', '%3.3f', '%3.3f', '%3.3f'};
r.header(k).default_value = z_max;
r.header(k).value = r.header(k).default_value;
r.header(k).validation = [];
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = [];
k = k + 1;

% Min Z, double, 8 bytes, *
r.header(k).compatibility = [10, 11, 12, 13, 14];
r.header(k).full_name = 'Min Z';
r.header(k).short_name = 'min_z';
r.header(k).type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).storage_type = {'double', 'double', 'double', 'double', 'double'};
r.header(k).byte_length = [8, 8, 8, 8, 8];
r.header(k).n_values = [1, 1, 1, 1, 1];
r.header(k).flag_bit_field = [false, false, false, false, false];
r.header(k).print_format = {'%3.3f', '%3.3f', '%3.3f', '%3.3f', '%3.3f'};
r.header(k).default_value = z_min;
r.header(k).value = r.header(k).default_value;
r.header(k).validation = [];
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = [];
k = k + 1;

% Start of Waveform Data Packet Record, Unsigned long long, 8 bytes, *
r.header(k).compatibility = [13, 14];
r.header(k).full_name = 'Start of Waveform Data Packet Record (LAS >= 1.3 only)';
r.header(k).short_name = 'offset_to_waveform_data';
r.header(k).type = {'uint64', 'uint64'};
r.header(k).storage_type = {'double', 'double'};
r.header(k).byte_length = [8, 8];
r.header(k).n_values = [1, 1];
r.header(k).flag_bit_field = [false, false];
r.header(k).print_format = {'%u', '%u'};
r.header(k).default_value = [];
r.header(k).value = [];
r.header(k).validation = @(x) isa(x, 'uint64') && (numel(x) == 1);
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = [];
k = k + 1;

% Start of first Extended Variable Length Record, unsigned long long, 8 bytes, *
r.header(k).compatibility = [14];
r.header(k).full_name = 'Start of first Extended Variable Length Record (LAS >= 1.4 only)';
r.header(k).short_name = 'offset_to_evlr';
r.header(k).type = {'uint64'};
r.header(k).storage_type = {'double'};
r.header(k).byte_length = [8];
r.header(k).n_values = [1];
r.header(k).flag_bit_field = [false];
r.header(k).print_format = {'%u'};
r.header(k).default_value = [];
r.header(k).value = [];
r.header(k).validation = @(x) isa(x, 'uint64') && (numel(x) == 1);
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = [];
k = k + 1;

% Number of Extended Variable Length Records, unsigned long, 4 bytes, *
r.header(k).compatibility = [14];
r.header(k).full_name = 'Number of Extended Variable Length Records (LAS >= 1.4 only)';
r.header(k).short_name = 'number_of_evlr';
r.header(k).type = {'uint32'};
r.header(k).storage_type = {'double'};
r.header(k).byte_length = [4];
r.header(k).n_values = [1];
r.header(k).flag_bit_field = [false];
r.header(k).print_format = {'%u'};
r.header(k).default_value = [];
r.header(k).value = [];
r.header(k).validation = @(x) isa(x, 'uint32') && (numel(x) == 1);
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = [];
k = k + 1;

% Number of point records, unsigned long long 8 bytes *
r.header(k).compatibility = [14];
r.header(k).full_name = 'Number of point records (LAS >= 1.4 only)';
r.header(k).short_name = 'n_point_records_extended';
r.header(k).type = {'uint64'};
r.header(k).storage_type = {'double'};
r.header(k).byte_length = [8];
r.header(k).n_values = [1];
r.header(k).flag_bit_field = [false];
r.header(k).print_format = {'%u'};
r.header(k).default_value = length(s.record.x);
r.header(k).value = r.header(k).default_value;
r.header(k).validation = [];
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = [];
k = k + 1;

% Number of points by return, unsigned long long [15], 120 bytes, *
r.header(k).compatibility = [14];
r.header(k).full_name = 'Number of points by return (LAS >= 1.4 only)';
r.header(k).short_name = 'n_points_by_return_extended'; % number_of_points_by_return
r.header(k).type = {'uint64'};
r.header(k).storage_type = {'double'};
r.header(k).byte_length = [8];
r.header(k).n_values = [15];
r.header(k).flag_bit_field = [false];
r.header(k).print_format = {'%u'};
r.header(k).default_value = [];
r.header(k).value = [];
r.header(k).validation = @(x) isa(x, 'uint64') && (numel(x) == 15);
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = [];


%% point data record format definition

k = 1;

% X, long, 4 bytes, *
r.record(k).full_name = 'X';
r.record(k).short_name = 'x';
r.record(k).version_compatibility = [10, 11, 12, 13, 14];
r.record(k).format_compatibility = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
r.record(k).bit_position = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32'};
r.record(k).storage_type = {'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32'};
r.record(k).bit_offset = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1 1];
r.record(k).bit_length = [32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = true;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% Y, long, 4 bytes, *
r.record(k).full_name = 'Y';
r.record(k).short_name = 'y';
r.record(k).version_compatibility = [10, 11, 12, 13, 14];
r.record(k).format_compatibility = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
r.record(k).bit_position = [32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32'};
r.record(k).storage_type = {'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32'};
r.record(k).bit_offset = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1 1];
r.record(k).bit_length = [32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = true;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% Z, long, 4 bytes, *
r.record(k).full_name = 'Z';
r.record(k).short_name = 'z';
r.record(k).version_compatibility = [10, 11, 12, 13, 14];
r.record(k).format_compatibility = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
r.record(k).bit_position = [64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32'};
r.record(k).storage_type = {'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32'};
r.record(k).bit_offset = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1 1];
r.record(k).bit_length = [32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = true;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% Intensity, unsigned short, 2 bytes
r.record(k).full_name = 'Intensity';
r.record(k).short_name = 'intensity';
r.record(k).version_compatibility = [10, 11, 12, 13, 14];
r.record(k).format_compatibility = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
r.record(k).bit_position = [96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16'};
r.record(k).storage_type = {'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16'};
r.record(k).bit_offset = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1 1];
r.record(k).bit_length = [16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = false;
r.record(k).flag_transform = false;
r.record(k).scale = 1;
r.record(k).offset = 0;
r.record(k).value = [];
k = k + 1;

% Return Number, 3 bits (bits 0, 1, 2), 3 bits, *
r.record(k).full_name = 'Return Number';
r.record(k).short_name = 'return_number';
r.record(k).version_compatibility = [10, 11, 12, 13, 14];
r.record(k).format_compatibility = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
r.record(k).bit_position = [112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'ubit3', 'ubit3', 'ubit3', 'ubit3', 'ubit3', 'ubit3', 'ubit4', 'ubit4', 'ubit4', 'ubit4', 'ubit4'};
r.record(k).storage_type = {'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8'};
r.record(k).bit_offset = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1 1];
r.record(k).bit_length = [3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = true;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = 1;
r.record(k).offset = 0;
r.record(k).value = [];
k = k + 1;

% Number of Returns (given pulse), 3 bits (bits 3, 4, 5), 3 bits, *
r.record(k).full_name = 'Number of Returns (given pulse)';
r.record(k).short_name = 'number_of_returns';
r.record(k).version_compatibility = [10, 11, 12, 13, 14];
r.record(k).format_compatibility = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
r.record(k).bit_position = [115, 115, 115, 115, 115, 115, 116, 116, 116, 116, 116];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'ubit3', 'ubit3', 'ubit3', 'ubit3', 'ubit3', 'ubit3', 'ubit4', 'ubit4', 'ubit4', 'ubit4', 'ubit4'};
r.record(k).storage_type = {'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8'};
r.record(k).bit_offset = [4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5];
r.record(k).bit_length = [3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4];
r.record(k).byte_length = r.record(k).bit_length / 8;r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = true;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = 1;
r.record(k).offset = 0;
r.record(k).value = [];
k = k + 1;

% Classification flag 1 - Synthetic, 1 bit (bit 0), 1 bit, *
r.record(k).full_name = 'Classification flag 1 - Synthetic';
r.record(k).short_name = 'classification_flag_1';
r.record(k).version_compatibility = [14];
r.record(k).format_compatibility = [6, 7, 8, 9, 10];
r.record(k).bit_position = [120, 120, 120, 120, 120];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1'};
r.record(k).storage_type = {'uint8', 'uint8', 'uint8', 'uint8', 'uint8'};
r.record(k).bit_offset = [1, 1, 1, 1, 1];
r.record(k).bit_length = [1, 1, 1, 1, 1];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = true;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = 1;
r.record(k).offset = 0;
r.record(k).value = [];
k = k + 1;

% Classification flag 2 - Key-point, 1 bit (bit 1), 1 bit, *
r.record(k).full_name = 'Classification flag 2 - Key-point';
r.record(k).short_name = 'classification_flag_2';
r.record(k).version_compatibility = [14];
r.record(k).format_compatibility = [6, 7, 8, 9, 10];
r.record(k).bit_position = [121, 121, 121, 121, 121];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1'};
r.record(k).storage_type = {'uint8', 'uint8', 'uint8', 'uint8', 'uint8'};
r.record(k).bit_offset = [2, 2, 2, 2, 2];
r.record(k).bit_length = [1, 1, 1, 1, 1];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = true;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = 1;
r.record(k).offset = 0;
r.record(k).value = [];
k = k + 1;

% Classification flag 3 - Withheld, 1 bit (bit 2), 1 bit, *
r.record(k).full_name = 'Classification flag 3 - Withheld';
r.record(k).short_name = 'classification_flag_3';
r.record(k).version_compatibility = [14];
r.record(k).format_compatibility = [6, 7, 8, 9, 10];
r.record(k).bit_position = [122, 122, 122, 122, 122];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1'};
r.record(k).storage_type = {'uint8', 'uint8', 'uint8', 'uint8', 'uint8'};
r.record(k).bit_offset = [3, 3, 3, 3, 3];
r.record(k).bit_length = [1, 1, 1, 1, 1];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = true;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = 1;
r.record(k).offset = 0;
r.record(k).value = [];
k = k + 1;

% Classification flag 4 - Overlap, 1 bit (bit 3), 1 bit, *
r.record(k).full_name = 'Classification flag 4 - Overlap';
r.record(k).short_name = 'classification_flag_4';
r.record(k).version_compatibility = [14];
r.record(k).format_compatibility = [6, 7, 8, 9, 10];
r.record(k).bit_position = [123, 123, 123, 123, 123];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1'};
r.record(k).storage_type = {'uint8', 'uint8', 'uint8', 'uint8', 'uint8'};
r.record(k).bit_offset = [4, 4, 4, 4, 4];
r.record(k).bit_length = [1, 1, 1, 1, 1];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = true;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = 1;
r.record(k).offset = 0;
r.record(k).value = [];
k = k + 1;

% Scanner channel, 2 bits (bit 4-5), 1 bit, *
r.record(k).full_name = 'Scanner channel';
r.record(k).short_name = 'scanner_channel';
r.record(k).version_compatibility = [14];
r.record(k).format_compatibility = [6, 7, 8, 9, 10];
r.record(k).bit_position = [124, 124, 124, 124, 124];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'ubit2', 'ubit2', 'ubit2', 'ubit2', 'ubit2'};
r.record(k).storage_type = {'uint8', 'uint8', 'uint8', 'uint8', 'uint8'};
r.record(k).bit_offset = [5, 5, 5, 5, 5];
r.record(k).bit_length = [2, 2, 2, 2, 2];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = true;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = 1;
r.record(k).offset = 0;
r.record(k).value = [];
k = k + 1;

% Scan Direction Flag, 1 bit (bit 6), 1 bit, *
r.record(k).full_name = 'Scan Direction Flag';
r.record(k).short_name = 'scan_direction_flag';
r.record(k).version_compatibility = [10, 11, 12, 13, 14];
r.record(k).format_compatibility = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
r.record(k).bit_position = [118, 118, 118, 118, 118, 118, 126, 126, 126, 126, 126];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1'};
r.record(k).storage_type = {'logical', 'logical', 'logical', 'logical', 'logical', 'logical', 'logical', 'logical', 'logical', 'logical', 'logical'};
r.record(k).bit_offset = [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7];
r.record(k).bit_length = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = true;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = 1;
r.record(k).offset = 0;
r.record(k).value = [];
k = k + 1;

% Edge of Flight Line, 1 bit (bit 7), 1 bit, *
r.record(k).full_name = 'Edge of Flight Line';
r.record(k).short_name = 'flightline_edge_flag';
r.record(k).version_compatibility = [10, 11, 12, 13, 14];
r.record(k).format_compatibility = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
r.record(k).bit_position = [119, 119, 119, 119, 119, 127, 127, 127, 127, 127];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1'};
r.record(k).storage_type = {'logical', 'logical', 'logical', 'logical', 'logical', 'logical', 'logical', 'logical', 'logical', 'logical', 'logical'};
r.record(k).bit_offset = [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8];
r.record(k).bit_length = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = true;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = 1;
r.record(k).offset = 0;
r.record(k).value = [];
k = k + 1;

% ASPRS classification number (0-31), 5 bits (bits 0, 1, 2, 3, 4), *
r.record(k).full_name = 'Classification (number)';
r.record(k).short_name = 'classification';
r.record(k).version_compatibility = [11, 12, 13, 14];
r.record(k).format_compatibility = [0, 1, 2, 3, 4, 5];
r.record(k).bit_position = [120, 120, 120, 120, 120, 120];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'ubit5', 'ubit5', 'ubit5', 'ubit5', 'ubit5', 'ubit5'};
r.record(k).storage_type = {'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8'};
r.record(k).bit_offset = [1, 1, 1, 1, 1, 1];
r.record(k).bit_length = [5, 5, 5, 5, 5, 5];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = true;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = 1;
r.record(k).offset = 0;
r.record(k).value = [];
k = k + 1;

% Synthetic - If set, then this point was created by a technique other than LIDAR collection such as digitized from a photogrammetric stereo model, 1 bit, *
r.record(k).full_name = 'Classification (synthetic flag)';
r.record(k).short_name = 'classification_synthetic';
r.record(k).version_compatibility = [11, 12, 13, 14];
r.record(k).format_compatibility = [0, 1, 2, 3, 4, 5];
r.record(k).bit_position = [125, 125, 125, 125, 125, 125];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1'};
r.record(k).storage_type = {'logical', 'logical', 'logical', 'logical', 'logical', 'logical'};
r.record(k).bit_offset = [6, 6, 6, 6, 6, 6];
r.record(k).bit_length = [1, 1, 1, 1, 1, 1];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = true;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = 1;
r.record(k).offset = 0;
r.record(k).value = [];
k = k + 1;

% Keypoint - If set, this point is considered to be a model key-point and thus generally should not be withheld in a thining algorithm, 1 bit, *
r.record(k).full_name = 'Classification (key-point flag)';
r.record(k).short_name = 'classification_keypoint';
r.record(k).version_compatibility = [11, 12, 13, 14];
r.record(k).format_compatibility = [0, 1, 2, 3, 4, 5];
r.record(k).bit_position = [126, 126, 126, 126, 126, 126];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1'};
r.record(k).storage_type = {'logical', 'logical', 'logical', 'logical', 'logical', 'logical'};
r.record(k).bit_offset = [7, 7, 7, 7, 7, 7];
r.record(k).bit_length = [1, 1, 1, 1, 1, 1];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = true;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = 1;
r.record(k).offset = 0;
r.record(k).value = [];
k = k + 1;

% Withheld - If set, this point should not be included in processing (synonymous with Deleted), 1 bit, *
r.record(k).full_name = 'Classification (withheld flag)';
r.record(k).short_name = 'classification_withheld';
r.record(k).version_compatibility = [11, 12, 13, 14];
r.record(k).format_compatibility = [0, 1, 2, 3, 4, 5];
r.record(k).bit_position = [127, 127, 127, 127, 127, 127];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1'};
r.record(k).storage_type = {'logical', 'logical', 'logical', 'logical', 'logical', 'logical'};
r.record(k).bit_offset = [8, 8, 8, 8, 8, 8];
r.record(k).bit_length = [1, 1, 1, 1, 1, 1];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = true;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = 1;
r.record(k).offset = 0;
r.record(k).value = [];
k = k + 1;

% Scan Angle, unsigned char, 1 byte, *
r.record(k).full_name = 'Scan Angle';
r.record(k).short_name = 'scan_angle';
r.record(k).version_compatibility = [10, 11, 12, 13, 14];
r.record(k).format_compatibility = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
r.record(k).bit_position = [128, 128, 128, 128, 128, 128, 144, 144, 144, 144, 144];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'int16', 'int16', 'int16', 'int16', 'int16'};
r.record(k).storage_type = {'uchar', 'uchar', 'uchar', 'uchar', 'uchar', 'uchar', 'single', 'single', 'single', 'single', 'single'};
r.record(k).bit_offset = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
r.record(k).bit_length = [8, 8, 8, 8, 8, 8, 16, 16, 16, 16, 16];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = 1;
r.record(k).offset = 0;
r.record(k).value = [];
k = k + 1;

% Classification, unsigned char, 1 byte, * (LAS 1.4)
r.record(k).full_name = 'Classification';
r.record(k).short_name = 'classification';
r.record(k).version_compatibility = [14];
r.record(k).format_compatibility = [6, 7, 8, 9, 10];
r.record(k).bit_position = [128, 128, 128, 128, 128];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'uint8', 'uint8', 'uint8', 'uint8', 'uint8'};
r.record(k).storage_type = {'uint8', 'uint8', 'uint8', 'uint8', 'uint8'};
r.record(k).bit_offset = [1, 1, 1, 1, 1];
r.record(k).bit_length = [8, 8, 8, 8, 8];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = 1;
r.record(k).offset = 0;
r.record(k).value = [];
k = k + 1;

% Classification, unsigned char, 1 byte, * (LAS 1.0)
r.record(k).full_name = 'Classification';
r.record(k).short_name = 'classification';
r.record(k).version_compatibility = [10];
r.record(k).format_compatibility = [0,1];
r.record(k).bit_position = [120, 120];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'uint8', 'uint8'};
r.record(k).storage_type = {'uint8', 'uint8'};
r.record(k).bit_offset = [1, 1];
r.record(k).bit_length = [8, 8];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = 1;
r.record(k).offset = 0;
r.record(k).value = [];
k = k + 1;

% User Data, unsigned char, 1 byte
r.record(k).full_name = 'User Data';
r.record(k).short_name = 'user_data';
r.record(k).version_compatibility = [10, 11, 12, 13, 14];
r.record(k).format_compatibility = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
r.record(k).bit_position = [136, 136, 136, 136, 136, 136, 136, 136, 136, 136, 136];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8'};
r.record(k).storage_type = {'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8'};
r.record(k).bit_offset = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
r.record(k).bit_length = [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = false;
r.record(k).flag_transform = false;
r.record(k).scale = 1;
r.record(k).offset = 0;
r.record(k).value = [];
k = k + 1;

% Point Source ID, unsigned short, 2 bytes, *
r.record(k).full_name = 'Point Source ID';
r.record(k).short_name = 'point_source_id';
r.record(k).version_compatibility = [10, 11, 12, 13, 14];
r.record(k).format_compatibility = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
r.record(k).bit_position = [144, 144, 144, 144, 144, 144, 160, 160, 160, 160, 160];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16'};
r.record(k).storage_type = {'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16'};
r.record(k).bit_offset = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
r.record(k).bit_length = [16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = 1;
r.record(k).offset = 0;
r.record(k).value = [];
k = k + 1;

% GPS Time, double, 8 bytes, *
r.record(k).full_name = 'GPS Time';
r.record(k).short_name = 'gps_time';
r.record(k).version_compatibility = [10, 11, 12, 13, 14];
r.record(k).format_compatibility = [1, 3, 4, 5, 6, 7, 8, 9, 10];
r.record(k).bit_position = [160, 160, 160, 160, 176, 176, 176, 176, 176];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'};
r.record(k).storage_type = {'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'};
r.record(k).bit_offset = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1,];
r.record(k).bit_length = [64, 64, 64, 64, 64, 64, 64, 64, 64, 64];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = 1;
r.record(k).offset = 0;
r.record(k).value = [];
k = k + 1;

% Red, unsigned short, 2 bytes, *
r.record(k).full_name = 'Red';
r.record(k).short_name = 'red';
r.record(k).version_compatibility = [12, 13, 14];
r.record(k).format_compatibility = [2, 3, 5, 7, 8, 10];
r.record(k).bit_position = [160, 224, 224, 240, 240, 240];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16'};
r.record(k).storage_type = {'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16'};
r.record(k).bit_offset = [1, 1, 1, 1, 1, 1];
r.record(k).bit_length = [16, 16, 16, 16, 16, 16];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = 1;
r.record(k).offset = 0;
r.record(k).value = [];
k = k + 1;

% Green, unsigned short, 2 bytes, *
r.record(k).full_name = 'Green';
r.record(k).short_name = 'green';
r.record(k).version_compatibility = [12, 13, 14];
r.record(k).format_compatibility = [2, 3, 5, 7, 8, 10];
r.record(k).bit_position = [176, 240, 240, 256, 256, 256];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16'};
r.record(k).storage_type = {'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16'};
r.record(k).bit_offset = [1, 1, 1, 1, 1, 1];
r.record(k).bit_length = [16, 16, 16, 16, 16, 16];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = 1;
r.record(k).offset = 0;
r.record(k).value = [];
k = k + 1;

% Blue, unsigned short, 2 bytes, *
r.record(k).full_name = 'Blue';
r.record(k).short_name = 'blue';
r.record(k).version_compatibility = [12, 13, 14];
r.record(k).format_compatibility = [2, 3, 5, 7, 8, 10];
r.record(k).bit_position = [192, 256, 256, 272, 272, 272];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16'};
r.record(k).storage_type = {'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16'};
r.record(k).bit_offset = [1, 1, 1, 1, 1, 1];
r.record(k).bit_length = [16, 16, 16, 16, 16, 16];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = 1;
r.record(k).offset = 0;
r.record(k).value = [];
k = k + 1;

% NIR, unsigned short, 2 bytes, *
r.record(k).full_name = 'Near infrared';
r.record(k).short_name = 'nir';
r.record(k).version_compatibility = [14];
r.record(k).format_compatibility = [8, 10];
r.record(k).bit_position = [288, 288];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'uint16', 'uint16'};
r.record(k).storage_type = {'uint16', 'uint16'};
r.record(k).bit_offset = [1, 1];
r.record(k).bit_length = [16, 16];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = 1;
r.record(k).offset = 0;
r.record(k).value = [];
k = k + 1;

% Wave Packet Descriptor Index, Unsigned char, 1 byte, *
r.record(k).full_name = 'Wave Packet Descriptor Index';
r.record(k).short_name = 'wave_packet_descriptor_index';
r.record(k).version_compatibility = [13, 14];
r.record(k).format_compatibility = [4, 5, 9, 10];
r.record(k).bit_position = [224, 272, 240, 304];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'uint8', 'uint8', 'uint8', 'uint8'};
r.record(k).storage_type = {'uint8', 'uint8', 'uint8', 'uint8'};
r.record(k).bit_offset = [1, 1, 1, 1];
r.record(k).bit_length = [8, 8, 8, 8];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = 1;
r.record(k).offset = 0;
r.record(k).value = [];
k = k + 1;

% Byte offset to waveform data, Unsigned long long, 8 bytes, *
r.record(k).full_name = 'Byte offset to waveform data';
r.record(k).short_name = 'offset_to_waveform_data';
r.record(k).version_compatibility = [13, 14];
r.record(k).format_compatibility = [4, 5, 9, 10];
r.record(k).bit_position = [232, 280, 248, 312];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'uint64', 'uint64', 'uint64', 'uint64'};
r.record(k).storage_type = {'uint64', 'uint64', 'uint64', 'uint64'};
r.record(k).bit_offset = [1, 1, 1, 1];
r.record(k).bit_length = [64, 64, 64, 64];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = 1;
r.record(k).offset = 0;
r.record(k).value = [];
k = k + 1;

% Waveform packet size in bytes, Unsigned long, 4 bytes, *
r.record(k).full_name = 'Waveform packet size in bytes';
r.record(k).short_name = 'waveform_packet_size';
r.record(k).version_compatibility = [13, 14];
r.record(k).format_compatibility = [4, 5, 9, 10];
r.record(k).bit_position = [296, 344, 312, 376];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'uint32', 'uint32', 'uint32', 'uint32'};
r.record(k).storage_type = {'uint32', 'uint32', 'uint32', 'uint32'};
r.record(k).bit_offset = [1, 1, 1, 1];
r.record(k).bit_length = [32, 32, 32, 32];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = 1;
r.record(k).offset = 0;
r.record(k).value = [];
k = k + 1;

% Return Point Waveform Location, float, 4 bytes, *
r.record(k).full_name = 'Return Point Waveform Location';
r.record(k).short_name = 'return_point_waveform_location';
r.record(k).version_compatibility = [13, 14];
r.record(k).format_compatibility = [4, 5, 9, 10];
r.record(k).bit_position = [328, 376, 344, 408];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'uint32', 'uint32', 'uint32', 'uint32'};
r.record(k).storage_type = {'uint32', 'uint32', 'uint32', 'uint32'};
r.record(k).bit_offset = [1, 1, 1, 1];
r.record(k).bit_length = [32, 32, 32, 32];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = 1;
r.record(k).offset = 0;
r.record(k).value = [];
k = k + 1;

% X(t), float, 4 bytes, *
r.record(k).full_name = 'X(t)';
r.record(k).short_name = 'x_t';
r.record(k).version_compatibility = [13, 14];
r.record(k).format_compatibility = [4, 5, 9, 10];
r.record(k).bit_position = [360, 408, 376, 440];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'single', 'single', 'single', 'single'};
r.record(k).storage_type = {'single', 'single', 'single', 'single'};
r.record(k).bit_offset = [1, 1, 1, 1];
r.record(k).bit_length = [32, 32, 32, 32];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = 1;
r.record(k).offset = 0;
r.record(k).value = [];
k = k + 1;

% Y(t), float, 4 bytes, *
r.record(k).full_name = 'Y(t)';
r.record(k).short_name = 'y_t';
r.record(k).version_compatibility = [13, 14];
r.record(k).format_compatibility = [4, 5, 9, 10];
r.record(k).bit_position = [392, 440, 408, 472];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'single', 'single', 'single', 'single'};
r.record(k).storage_type = {'single', 'single', 'single', 'single'};
r.record(k).bit_offset = [1, 1, 1, 1];
r.record(k).bit_length = [32, 32, 32, 32];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = 1;
r.record(k).offset = 0;
r.record(k).value = [];
k = k + 1;

% Z(t), float, 4 bytes, *
r.record(k).full_name = 'Z(t)';
r.record(k).short_name = 'z_t';
r.record(k).version_compatibility = [13, 14];
r.record(k).format_compatibility = [4, 5, 9, 10];
r.record(k).bit_position = [424, 472, 440, 504];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'single', 'single', 'single', 'single'};
r.record(k).storage_type = {'single', 'single', 'single', 'single'};
r.record(k).bit_offset = [1, 1, 1, 1];
r.record(k).bit_length = [32, 32, 32, 32];
r.record(k).byte_length = r.record(k).bit_length / 8;
r.record(k).flag_extra_field = false;
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = 1;
r.record(k).offset = 0;
r.record(k).value = [];


%% extra bytes format definition

k = 1;
r.extra_bytes(k).id = 0;
r.extra_bytes(k).type = 'NaN';
r.extra_bytes(k).n_values = 1;
r.extra_bytes(k).byte_length = 'NaN';
k = k + 1;

r.extra_bytes(k).id = 1;
r.extra_bytes(k).type = 'uint8'; % uchar
r.extra_bytes(k).n_values = 1;
r.extra_bytes(k).byte_length = 1;
k = k + 1;

r.extra_bytes(k).id = 2;
r.extra_bytes(k).type = 'int8'; % char
r.extra_bytes(k).n_values = 1;
r.extra_bytes(k).byte_length = 1;
k = k + 1;

r.extra_bytes(k).id = 3;
r.extra_bytes(k).type = 'uint16';
r.extra_bytes(k).n_values = 1;
r.extra_bytes(k).byte_length = 2;
k = k + 1;

r.extra_bytes(k).id = 4;
r.extra_bytes(k).type = 'int16';
r.extra_bytes(k).n_values = 1;
r.extra_bytes(k).byte_length = 2;
k = k + 1;

r.extra_bytes(k).id = 5;
r.extra_bytes(k).type = 'uint32';
r.extra_bytes(k).n_values = 1;
r.extra_bytes(k).byte_length = 4;
k = k + 1;

r.extra_bytes(k).id = 6;
r.extra_bytes(k).type = 'int32';
r.extra_bytes(k).n_values = 1;
r.extra_bytes(k).byte_length = 4;
k = k + 1;

r.extra_bytes(k).id = 7;
r.extra_bytes(k).type = 'uint64';
r.extra_bytes(k).n_values = 1;
r.extra_bytes(k).byte_length = 8;
k = k + 1;

r.extra_bytes(k).id = 8;
r.extra_bytes(k).type = 'int64';
r.extra_bytes(k).n_values = 1;
r.extra_bytes(k).byte_length = 8;
k = k + 1;

r.extra_bytes(k).id = 9;
r.extra_bytes(k).type = 'single';
r.extra_bytes(k).n_values = 1;
r.extra_bytes(k).byte_length = 4;
k = k + 1;

r.extra_bytes(k).id = 10;
r.extra_bytes(k).type = 'double';
r.extra_bytes(k).n_values = 1;
r.extra_bytes(k).byte_length = 8;
k = k + 1;

r.extra_bytes(k).id = 11;
r.extra_bytes(k).type = 'uchar';
r.extra_bytes(k).n_values = 2;
r.extra_bytes(k).byte_length = 1;
k = k + 1;

r.extra_bytes(k).id = 12;
r.extra_bytes(k).type = 'char';
r.extra_bytes(k).n_values = 2;
r.extra_bytes(k).byte_length = 1;
k = k + 1;

r.extra_bytes(k).id = 13;
r.extra_bytes(k).type = 'uint16';
r.extra_bytes(k).n_values = 2;
r.extra_bytes(k).byte_length = 2;
k = k + 1;

r.extra_bytes(k).id = 14;
r.extra_bytes(k).type = 'int16';
r.extra_bytes(k).n_values = 2;
r.extra_bytes(k).byte_length = 2;
k = k + 1;

r.extra_bytes(k).id = 15;
r.extra_bytes(k).type = 'uint32';
r.extra_bytes(k).n_values = 2;
r.extra_bytes(k).byte_length = 4;
k = k + 1;

r.extra_bytes(k).id = 16;
r.extra_bytes(k).type = 'int32';
r.extra_bytes(k).n_values = 2;
r.extra_bytes(k).byte_length = 4;
k = k + 1;

r.extra_bytes(k).id = 17;
r.extra_bytes(k).type = 'uint8';
r.extra_bytes(k).n_values = 2;
r.extra_bytes(k).byte_length = 8;
k = k + 1;

r.extra_bytes(k).id = 18;
r.extra_bytes(k).type = 'int64';
r.extra_bytes(k).n_values = 2;
r.extra_bytes(k).byte_length = 8;
k = k + 1;

r.extra_bytes(k).id = 19;
r.extra_bytes(k).type = 'single';
r.extra_bytes(k).n_values = 2;
r.extra_bytes(k).byte_length = 4;
k = k + 1;

r.extra_bytes(k).id = 20;
r.extra_bytes(k).type = 'double';
r.extra_bytes(k).n_values = 2;
r.extra_bytes(k).byte_length = 8;
k = k + 1;

r.extra_bytes(k).id = 21;
r.extra_bytes(k).type = 'uchar';
r.extra_bytes(k).n_values = 3;
r.extra_bytes(k).byte_length = 1;
k = k + 1;

r.extra_bytes(k).id = 22;
r.extra_bytes(k).type = 'char';
r.extra_bytes(k).n_values = 3;
r.extra_bytes(k).byte_length = 1;
k = k + 1;

r.extra_bytes(k).id = 23;
r.extra_bytes(k).type = 'uint16';
r.extra_bytes(k).n_values = 3;
r.extra_bytes(k).byte_length = 2;
k = k + 1;

r.extra_bytes(k).id = 24;
r.extra_bytes(k).type = 'int16';
r.extra_bytes(k).n_values = 3;
r.extra_bytes(k).byte_length = 2;
k = k + 1;

r.extra_bytes(k).id = 25;
r.extra_bytes(k).type = 'uint32';
r.extra_bytes(k).n_values = 3;
r.extra_bytes(k).byte_length = 4;
k = k + 1;

r.extra_bytes(k).id = 26;
r.extra_bytes(k).type = 'int32';
r.extra_bytes(k).n_values = 3;
r.extra_bytes(k).byte_length = 4;
k = k + 1;

r.extra_bytes(k).id = 27;
r.extra_bytes(k).type = 'uint64';
r.extra_bytes(k).n_values = 3;
r.extra_bytes(k).byte_length = 8;
k = k + 1;

r.extra_bytes(k).id = 28;
r.extra_bytes(k).type = 'int64';
r.extra_bytes(k).n_values = 3;
r.extra_bytes(k).byte_length = 8;
k = k + 1;

r.extra_bytes(k).id = 29;
r.extra_bytes(k).type = 'single';
r.extra_bytes(k).n_values = 3;
r.extra_bytes(k).byte_length = 4;
k = k + 1;

r.extra_bytes(k).id = 30;
r.extra_bytes(k).type = 'double';
r.extra_bytes(k).n_values = 3;
r.extra_bytes(k).byte_length = 8;


%% adjust header format based on the LAS version

% filter version specific attributes
idx_header_format = false(length(r.header),1);

for j = 1:length(r.header)
    
    [flag_include_header, ind_header] = ismember(las_version, r.header(j).compatibility);
    
    if flag_include_header
        
        idx_header_format(j) = true;
        
        r.header(j).type = r.header(j).type{ind_header};
        r.header(j).byte_length = r.header(j).byte_length(ind_header);
        r.header(j).n_values = r.header(j).n_values(ind_header);
        r.header(j).flag_bit_field = r.header(j).flag_bit_field(ind_header);
        r.header(j).print_format = r.header(j).print_format{ind_header};
        
        if r.header(j).flag_bit_field
            
            r.header(j).bit_position = r.header(j).bit_position(ind_header);
            r.header(j).parent_type = r.header(j).parent_type{ind_header};
            
        end
        
        if isempty(r.header(j).value)
            
            try

                r.header(j).value = s.header.(r.header(j).short_name);
                
            catch
                
                r.header(j).value = r.header(j).default_value;
                
                if arg.Results.verbose
                    
                    fprintf('WARNING: Setting %s to default: %.3f\n', r.header(j).full_name, r.header(j).default_value);
                    
                end
                
            end
            
        end
        
    end
    
end

r.header(~idx_header_format) = [];


%% generate key-value pairs for public header block

for j = 1:length(r.header)
    
    phb_skeys.(r.header(j).short_name) = j;
    
end


%% adjust point record based on the Point Data Format ID (0-99 for spec)

idxl_extra_byte = ~ismember(fieldnames(s.record), {r.record.short_name});
idxl_valid_version = cellfun(@(x) ismember(arg.Results.version, x), {r.record.version_compatibility});
idxl_valid_format = cellfun(@(x) ismember(arg.Results.recordFormat, x), {r.record.format_compatibility});

r.record = r.record(idxl_valid_version & idxl_valid_format);

for j = 1:length(r.record)
    
    [~, ind_record] = ismember(arg.Results.recordFormat, r.record(j).format_compatibility);
    
    r.record(j).format_compatibility = r.record(j).format_compatibility(ind_record);
    r.record(j).bit_position = r.record(j).bit_position(ind_record);
    r.record(j).byte_position = r.record(j).byte_position(ind_record);
    r.record(j).type = r.record(j).type(ind_record);
    r.record(j).storage_type = r.record(j).storage_type(ind_record);
    r.record(j).bit_offset = r.record(j).bit_offset(ind_record);
    r.record(j).bit_length = r.record(j).bit_length(ind_record);
    r.record(j).byte_length = r.record(j).byte_length(ind_record);
    
end

[~, ind_record_sort] = sort([r.record.bit_position]);
ind_last_record = ind_record_sort(end);
ind_record = length(r.record); % index of last record


%% add extra bytes to variable length headers and records

record_abbreviations = fieldnames(s.record);

flag_extra_byte = any(idxl_extra_byte);

if flag_extra_byte
    
    if arg.Results.verbose
        
        fprintf('WARNING: The record contains extra bytes\n');
        
    end
    
    extra_record_abbreviations = record_abbreviations(idxl_extra_byte);
    ind_extra_records = find([s.variable_length_records.record_id] == 4);
    
    for j = 1:length(ind_extra_records)
  
        for k = 1:length(s.variable_length_records(ind_extra_records(j)).value)
            
            ind_record = ind_record + 1;
            ind_last_record = ind_last_record + 1;
            ind_data_type = s.variable_length_records(ind_extra_records(j)).value(k).data_type + 1;
            
            r.record(ind_record).full_name = s.variable_length_records(ind_extra_records(j)).value(k).name;
            r.record(ind_record).short_name = extra_record_abbreviations{k};
            r.record(ind_record).format_compatibility = r.record(ind_record - 1).format_compatibility;
            r.record(ind_record).bit_position = r.record(ind_last_record - 1).bit_position + (r.record(ind_last_record - 1).byte_length * 8);
            r.record(ind_record).byte_position =  r.record(ind_record).bit_position ./ 8;
            r.record(ind_record).bit_offset = 1;
            
            r.record(ind_record).type = {r.extra_bytes(ind_data_type).type};
            r.record(ind_record).storage_type = {r.extra_bytes(ind_data_type).type};
            r.record(ind_record).byte_length = r.extra_bytes(ind_data_type).n_values * r.extra_bytes(ind_data_type).byte_length;
            r.record(ind_record).bit_length = r.record(ind_record).byte_length * 8;

            r.record(ind_record).flag_bit_field = false;
            r.record(ind_record).flag_mandatory = false;
            r.record(ind_record).flag_transform = true;
            
            % scale and offset transformations
            if  s.variable_length_records(ind_extra_records(j)).value(k).scale(1) == 0
                
                r.record(ind_record).scale = 1;
                
            else
                
                r.record(ind_record).scale = s.variable_length_records(ind_extra_records(j)).value(k).scale(1);
                
            end
            
            r.record(ind_record).offset = s.variable_length_records(ind_extra_records(j)).value(k).offset(1);
            
            
            r.record(ind_record).value = [];
            
        end
        
    end
    
end


%% generate key-value pairs for point data record

for j = 1:length(r.record)
    
    pdr_skeys.(r.record(j).short_name) = j;
    
end

%pdr_keys = fieldnames(r.record);


%% populate point data records

for j = 1:length(r.record)
    
    try
        
        r.record(j).value = s.record.(r.record(j).short_name);

    catch
        
        % fill with zeros, if the record values are not specified
        r.record(j).value = cast(zeros(size(s.record.x)), r.record(j).storage_type{:});
        fprintf('WARNING: the point data record "%s" (%s) was not specified, setting to 0\n', r.record(j).short_name, r.record(j).full_name);
        
    end
    
end


%% check point data record

% check availability of mandatory record fields
mandatory_fields = {r.record([r.record.flag_mandatory]).short_name};
available_fields = ismember(mandatory_fields, {r.record.short_name}');

if ~all(available_fields)
    
    error('LASwrite:missingMandatoryField', 'The following mandatory fields are missing: %s', strjoin(mandatory_fields(~available_fields)',', '));
    
end

% update array length
n_records_per_attribute = zeros(length({r.record.full_name}),1);

for k = 1:length(n_records_per_attribute)
    
    n_records_per_attribute(k,1) = length(r.record(k).value);
    
end

if ~(range(n_records_per_attribute) == 0)
    
    error('LASwrite:unequalRecordLengths', 'The lengths of point data records are unequal');
    
end


%% update public header block values

for j = 1:length(r.header)
    
    switch r.header(j).short_name
        
        case 'n_variable_length_records'
            
            try
                
                r.header(j).value = length(s.variable_length_records); % Number of Variable Length Records, unsigned long, 4 bytes, *
                
            catch
                
                fprintf(r.header(j).warning_message);
                r.header(j).value = 0;
                
            end
            
        case 'offset_to_data'
            
            try
                n_variable_length_records = length(s.variable_length_records);
                
                if n_variable_length_records > 0
                    
                    variable_header_size = n_variable_length_records * 54 + sum([s.variable_length_records.record_length_after_header]);
                    
                else
                    
                    variable_header_size = 0;
                    
                end
                
                r.header(j).value = r.header(phb_skeys.header_size).value + variable_header_size;
                
            catch
                
                r.header(j).value = r.header(phb_skeys.header_size).value;
                
            end
            
        case 'point_data_record_length'
            
            point_data_format_size = sum([r.record.byte_length]); % ADD EXTRA BYTES
            r.header(j).value = point_data_format_size; % Point Data Record Length, unsigned short, 2 bytes, *
            
        case {'n_points_by_return', 'n_points_by_return_extended'}
            
            for k = 1:r.header(j).n_values
                
                r.header(j).value(k) = sum(r.record(pdr_skeys.return_number).value == k); % Number of points by return N
                
            end
            
    end
    
end

flag_evlr = isfield(s, 'extended_variable_length_records') & (las_version >= 14);

if flag_evlr
    
    for j = 1:length(r.header)
        
        switch r.header(j).short_name
            
            case 'offset_to_evlr'
                
                r.header(j).value = r.header(phb_skeys.offset_to_data).value + r.header(phb_skeys.point_data_record_length).value * r.header(phb_skeys.n_point_records).value;
                
            case 'number_of_evlr'
                
                r.header(j).value = length(s.extended_variable_length_records);
                
        end
        
    end
    
end

clear s.record
clear s.header


%% apply scale and offset when necessary

r.record(pdr_skeys.x).scale = r.header(phb_skeys.x_scale_factor).value;
r.record(pdr_skeys.x).offset = r.header(phb_skeys.x_offset).value;
r.record(pdr_skeys.y).scale = r.header(phb_skeys.y_scale_factor).value;
r.record(pdr_skeys.y).offset = r.header(phb_skeys.y_offset).value;
r.record(pdr_skeys.z).scale = r.header(phb_skeys.z_scale_factor).value;
r.record(pdr_skeys.z).offset = r.header(phb_skeys.z_offset).value;

for j = 1:length(r.record)
    
    if any([r.record(j).scale ~= 1,  r.record(j).offset ~= 0])
    %if r.record(j).flag_transform
        
        r.record(j).value = cast((r.record(j).value - r.record(j).offset) / r.record(j).scale, r.record(j).type{:});
        
    end
    
end


%% adjust scan angle format

switch r.record(pdr_skeys.scan_angle).type{:}
    
    case 'uint8'
        
        % Scan Angle Rank is the scan angle rounded to the nearest integer from from +90° to -90°,
        % with -ve being left of 0° (nadir) looking forward
        r.record(pdr_skeys.scan_angle).value = int16(r.record(pdr_skeys.scan_angle).value);
        idx_scan_angle_sign = r.record(pdr_skeys.scan_angle).value < 0;
        r.record(pdr_skeys.scan_angle).value(idx_scan_angle_sign) = r.record(pdr_skeys.scan_angle).value(idx_scan_angle_sign) + 256;
        
    case 'int16'
        
        % Scan Angle is a rotational position in 0.006° increments with values from +30,000 to -30,000
        % covering +180° to -180° (0° is at nadir/down, -ve is left/clockwise of nadir looking forward)
        r.record(pdr_skeys.scan_angle).value = r.record(pdr_skeys.scan_angle).value * 15000 / 90;
        
end


%% reformat UUID

% TEST
if ismember('uuid', fieldnames(pdr_skeys))
    
    [uuid_ref, ~, ic] = unique(s.record.uuid);
    uuid_ref_num = zeros(length(uuid_ref), 16, 'uint8');
    
    for j = 1:length(uuid_ref)
        
        uuid_ref_num(j,:) = uint8(hex2dec([uuid_ref{j}(1:2:end)', uuid_ref{j}(2:2:end)']))';
        
    end
    
    r.record(pdr_skeys.uuid).value = uuid_ref_num(ic,:);
    
end

%% write structure to LAS file

if ~isempty(arg.Results.filepath)
    
    % open file
    fid = fopen(filepath, 'w+', MACHINE_FORMAT);
    
    if fid == -1
        
        error('Could not open file at specified location')
        
    end
    
    %% write public header to file
    
    byte_offset = 0;
    integer_byte_offset = 0;
    fwrite(fid, zeros(r.header(phb_skeys.header_size).value, 1, 'uint8'), 'uint8', 0, MACHINE_FORMAT); % prefill with zeros
    
    % write public header content
    for j = 1:length(r.header)
        
        switch r.header(j).type
            
            case {'uchar', 'char'}
                
                fseek(fid, integer_byte_offset, 'bof');
                string = horzcat(r.header(j).value, char(zeros(1,r.header(j).n_values))); % add null terminator
                fwrite(fid, string(1:r.header(j).n_values), r.header(j).type, 0, MACHINE_FORMAT);
                
            case {'uint8', 'uint16', 'ulong', 'uint32', 'uint64', 'double'}
                
                fseek(fid, integer_byte_offset, 'bof');
                fwrite(fid, cast(r.header(j).value, r.header(j).type), r.header(j).type, 0, MACHINE_FORMAT);
                
            otherwise % pack bits
                
                fseek(fid, integer_byte_offset, 'bof');
                B = fread(fid, 1, r.header(j).parent_type, 0, MACHINE_FORMAT) + bitset(zeros(1, r.header(j).parent_type), r.header(j).bit_position+1, r.header(j).value);
                fseek(fid, integer_byte_offset, 'bof');
                fwrite(fid, B, r.header(j).parent_type, 0, MACHINE_FORMAT);
                
        end
        
        byte_offset = byte_offset + r.header(j).n_values * r.header(j).byte_length;
        
        if mod(byte_offset,1) == 0
            
            integer_byte_offset = byte_offset;
            
        end
        
        if arg.Results.verbose
            
            fprintf(sprintf('%s: %s\\n', r.header(j).full_name, r.header(j).print_format), r.header(j).value);
            
        end
        
    end
    
    
    %% write variable length records to file
    
    flag_vlr = (r.header(phb_skeys.n_variable_length_records).value > 0);
    byte_offset = r.header(phb_skeys.header_size).value;
    
    if flag_vlr
        
        for j = 1:r.header(phb_skeys.n_variable_length_records).value
            
            fseek(fid, byte_offset, 'bof');
            
            % write common header part
            fwrite(fid, s.variable_length_records(j).reserved, 'uint16', 0, MACHINE_FORMAT); % Reserved, unsigned short, 2 bytes;
            user_id = horzcat(s.variable_length_records(j).user_id, char(zeros(1, 16))); % zero padding % TODO
            fwrite(fid, user_id(1:16), 'char', 0, MACHINE_FORMAT); % User ID, char[16], 16 bytes, *
            fwrite(fid, s.variable_length_records(j).record_id , 'uint16', 0, MACHINE_FORMAT); % Record ID, unsigned short, 2 bytes, *
            fwrite(fid, s.variable_length_records(j).record_length_after_header , 'uint16', 0, MACHINE_FORMAT); % Record Length After Header, unsigned short, 2 bytes, *
            description = horzcat(s.variable_length_records(j).description, char(zeros(1, 32))); % zero padding % TODO
            fwrite(fid, description(1:32) , 'char', 0, MACHINE_FORMAT); % Description, char[32], 32 bytes
            
            % write defined variable length header entries
            switch s.variable_length_records(j).record_id
                
                case 2111 % OGC Math Transform WKT Record (optional)
                    
                    string = horzcat(s.variable_length_records(j).value, char(zeros(1, s.variable_length_records(j).record_length_after_header))); % zero padding
                    fwrite(fid, string(1:s.variable_length_records(j).record_length_after_header), 'char', 0, MACHINE_FORMAT);
                    
                case 2112 % OGC Coordinate System WKT Record (optional)
                    
                    string = horzcat(s.variable_length_records(j).value, char(zeros(1, s.variable_length_records(j).record_length_after_header))); % zero padding
                    fwrite(fid, string(1:s.variable_length_records(j).record_length_after_header), 'char', 0, MACHINE_FORMAT);
                    
                case 34735 % GeoKeyDirectoryTag Record (mandatory)
                    
                    s.variable_length_records(j).value.w_key_directory_version = 1; % wKeyDirectoryVersion, unsigned short, wKeyDirectoryVersion = 1 Always
                    s.variable_length_records(j).value.w_key_revision = 1; % wKeyRevision, unsigned short, wKeyRevision = 1 Always
                    s.variable_length_records(j).value.w_minor_revision = 0; % wMinorRevision, unsigned short, wMinorRevision = 0 Always
                    s.variable_length_records(j).value.w_number_of_keys = length(s.variable_length_records(j).value.s_key_entry); % wNumberOfKeys, unsigned short, Number of sets of 4 unsigned shorts to follow
                    
                    fwrite(fid, s.variable_length_records(j).value.w_key_directory_version, 'uint16', 0, MACHINE_FORMAT);
                    fwrite(fid, s.variable_length_records(j).value.w_key_revision, 'uint16', 0, MACHINE_FORMAT);
                    fwrite(fid, s.variable_length_records(j).value.w_minor_revision, 'uint16', 0, MACHINE_FORMAT);
                    fwrite(fid, s.variable_length_records(j).value.w_number_of_keys, 'uint16', 0, MACHINE_FORMAT);
                    
                    for i=1:s.variable_length_records(j).value.w_number_of_keys
                        
                        fwrite(fid, s.variable_length_records(j).value.s_key_entry{i}(1), 'uint16', 0, MACHINE_FORMAT);
                        fwrite(fid, s.variable_length_records(j).value.s_key_entry{i}(2), 'uint16', 0, MACHINE_FORMAT);
                        fwrite(fid, s.variable_length_records(j).value.s_key_entry{i}(3), 'uint16', 0, MACHINE_FORMAT);
                        fwrite(fid, s.variable_length_records(j).value.s_key_entry{i}(4), 'uint16', 0, MACHINE_FORMAT);
                        
                    end
                    
                case 34736 % GeoDoubleParamsTag Record (optional)
                    
                    fwrite(fid, s.variable_length_records(j).value, 'double', 0, MACHINE_FORMAT);
                    
                case 34737 % GeoAsciiParamsTag Record (optional)
                    
                    string = horzcat(s.variable_length_records(j).value, char(zeros(1, s.variable_length_records(j).record_length_after_header))); % zero padding
                    fwrite(fid, string(1:s.variable_length_records(j).record_length_after_header), 'char', 0, MACHINE_FORMAT);
                    
                case 0 % Classification lookup (optional)
                    
                case 1 % Reserved (optional)
                    
                case 2 % Histogram (optional)
                    
                case 3 % Text area description (optional)
                    
                    string = horzcat(s.variable_length_records(j).value, char(zeros(1, s.variable_length_records(j).record_length_after_header))); % zero padding
                    fwrite(fid, string(1:s.variable_length_records(j).record_length_after_header), 'char', 0, MACHINE_FORMAT);
                    
                case 4 % Extra bytes (optional)
                    
                    n_extra_records = s.variable_length_records(j).record_length_after_header / 192;
                    
                    for k = 1:n_extra_records
                        
                        % reserved, unsigned char, 2 bytes
                        fwrite(fid, 0, 'uint16', 0, MACHINE_FORMAT);
                        
                        % data type, unsigned char, 1 byte
                        fwrite(fid, s.variable_length_records(j).value(k).data_type, 'uint8', 0, MACHINE_FORMAT);
                        
                        % options, unsigned char, 1 byte
                        extra_byte_options = [s.variable_length_records(j).value(k).options.no_data_bit,...
                            s.variable_length_records(j).value(k).options.min_bit,...
                            s.variable_length_records(j).value(k).options.max_bit,...
                            s.variable_length_records(j).value(k).options.scale_bit,...
                            s.variable_length_records(j).value(k).options.offset_bit];
                        
                        if any(extra_byte_options)
                            
                            s.variable_length_records(j).value(k).options = sum(bitset(zeros(1, nnz(extra_byte_options), 'uint8'), find(extra_byte_options)+0));
                            
                        else
                            
                            s.variable_length_records(j).value(k).options = uint8(0);
                            
                        end
                        
                        fwrite(fid, s.variable_length_records(j).value(k).options, 'uint8', 0, MACHINE_FORMAT);
                        
                        % name, char, 32 bytes
                        name = horzcat(s.variable_length_records(j).value(k).name, char(zeros(1, 32))); % zero padding
                        fwrite(fid, name(1:32), 'char', 0, MACHINE_FORMAT);
                        
                        % unused, unsigned char, 4 bytes
                        fwrite(fid, 0, 'uint32', 0, MACHINE_FORMAT);
                        
                        % no data, anytype, 24 bytes
                        fwrite(fid, s.variable_length_records(j).value(k).no_data, 'double', 0, MACHINE_FORMAT);
                        
                        % min, anytype, 24 bytes
                        fwrite(fid, s.variable_length_records(j).value(k).min, 'double', 0, MACHINE_FORMAT);
                        
                        % max, anytype, 24 bytes
                        fwrite(fid, s.variable_length_records(j).value(k).max, 'double', 0, MACHINE_FORMAT);
                        
                        % scale, double, 24 bytes
                        fwrite(fid, s.variable_length_records(j).value(k).scale, 'double', 0, MACHINE_FORMAT);
                        
                        % offset, double, 24 bytes
                        fwrite(fid, s.variable_length_records(j).value(k).offset, 'double', 0, MACHINE_FORMAT);
                        
                        % description, char, 32 bytes
                        description = horzcat(s.variable_length_records(j).value(k).description, char(zeros(1, 32))); % zero padding
                        fwrite(fid, description(1:32), 'char', 0, MACHINE_FORMAT);
                        
                    end
                    
                case num2cell(100:354) % Waveform Packet descriptor (required when using point formats 4, 5, 9, 10)
                    
                    fwrite(fid, s.variable_length_records(j).value.bits_per_sample, 'uint8', 0, MACHINE_FORMAT); % Bits per sample, Unsigned char, 1 byte, *
                    fwrite(fid, s.variable_length_records(j).value.compression_type, 'uint8', 0, MACHINE_FORMAT); % Waveform compression type, Unsigned char, 1 byte, *
                    fwrite(fid, s.variable_length_records(j).value.number_samples, 'uint32', 0, MACHINE_FORMAT); % Number of samples, Unsigned long, 4 bytes, *
                    fwrite(fid, s.variable_length_records(j).value.temporal_sample_spacing, 'uint32', 0, MACHINE_FORMAT); % Temporal Sample Spacing, Unsigned long, 4 bytes, *
                    fwrite(fid, s.variable_length_records(j).value.digitizer_gain, 'double', 0, MACHINE_FORMAT); % Digitizer Gain, double, 8 bytes, *
                    fwrite(fid, s.variable_length_records(j).value.digitizer_offset, 'double', 0, MACHINE_FORMAT); % Digitizer Offset, double, 8 bytes, *
                    
                case 30001
                    
                    fwrite(fid, s.variable_length_records(j).value, 'uint32', 0, MACHINE_FORMAT);
                    
                case 30002
                    
                    fwrite(fid, s.variable_length_records(j).value, 'uint32', 0, MACHINE_FORMAT);
                    
                    % You may add custom records here
                    
                otherwise % Other
                    
                    fwrite(fid, typecast(s.variable_length_records(j).value, 'uint8'), 'uint8', 0, MACHINE_FORMAT);

            end
            
            byte_offset = byte_offset + 54 + s.variable_length_records(j).record_length_after_header;
            
        end
        
    end
    
    
    %% write point data record to file
    
    if arg.Results.verbose
        
        tic; % start timer
        fprintf('writing point data record...');
        
    end
    
    % create blob
    n_bytes = r.header(phb_skeys.n_point_records).value * r.header(phb_skeys.point_data_record_length).value;
    blob = zeros(n_bytes, 1, 'uint8');
    
    byte = zeros(r.header(phb_skeys.n_point_records).value, 1, 'uint8');
    bit_count = 0;
    
    for j = 1:size(r.record,2)

        if r.record(j).flag_bit_field % pack bit fields into byte

            byte = bitor(bitshift(uint8(r.record(j).value), bit_count), byte);
            bit_count = bit_count + r.record(j).bit_length;
            
            if mod(bit_count, 8) == 0
                
                % check if array is larger than max memory size
                ind_byte = reshape(cell2mat(arrayfun(@(x) x:r.header(phb_skeys.point_data_record_length).value:n_bytes, floor(r.record(j).byte_position) + 1, 'UniformOutput', false)), 1, []);
                
                blob(ind_byte) = typecast(cast(byte, 'uint8'), 'uint8'); % add values to blob part when byte is full
                bit_count = 0; % reinitialize bit count
                byte = zeros(r.header(phb_skeys.n_point_records).value, 1, 'uint8'); % reinitialize byte to zero
                
            end
            
        else
            
            seq = (r.record(j).byte_position + 1:r.header(phb_skeys.point_data_record_length).value:n_bytes)';
            ind_byte = (seq * ones(1, r.record(j).byte_length) + repmat(0:r.record(j).byte_length-1, r.header(phb_skeys.n_point_records).value, 1))';
            ind_byte = ind_byte(:);
            
            blob(ind_byte) = typecast(cast(reshape(r.record(j).value',[],1), r.record(j).type{:}), 'uint8');

            bit_count = 0;
            
        end
          
    end
    
    % write blob to file
    fseek(fid, r.header(phb_skeys.offset_to_data).value, 'bof');
    fwrite(fid, blob, 'uint8');
    fprintf('done!\n');
    
    %% write extended variable length records to file

        if flag_evlr
            
            byte_offset = ftell(fid); %r.header(phb_skeys.offset_to_evlr).value;
            
            for j = 1:r.header(phb_skeys.number_of_evlr).value
                
                fseek(fid, byte_offset, 'bof');
                
                % write common header part
                fwrite(fid, s.extended_variable_length_records(j).reserved, 'uint16', 0, MACHINE_FORMAT); % Reserved, unsigned short, 2 bytes;
                user_id = horzcat(s.extended_variable_length_records(j).user_id, char(zeros(1, 16))); % zero padding
                fwrite(fid, user_id(1:16), 'char', 0, MACHINE_FORMAT); % User ID, char[16], 16 bytes, *
                fwrite(fid, s.extended_variable_length_records(j).record_id , 'uint16', 0, MACHINE_FORMAT); % Record ID, unsigned short, 2 bytes, *
                fwrite(fid, s.extended_variable_length_records(j).record_length_after_header , 'uint64', 0, MACHINE_FORMAT); % Record Length After Header, unsigned short, 2 bytes, *
                description = horzcat(s.extended_variable_length_records(j).description, char(zeros(1, 32))); % zero padding
                fwrite(fid, description(1:32) , 'char', 0, MACHINE_FORMAT); % Description, char[32], 32 bytes
                
                % write defined variable length header entries
                switch s.extended_variable_length_records(j).record_id
                    
                    case 0 % Classification lookup (optional)
                        
                    case 1 % Reserved (optional)
                        
                    case 2 % Histogram (optional)
                        
                    case 3 % Text area description (optional)
                        
                        string = horzcat(s.extended_variable_length_records(j).value, char(zeros(1, s.extended_variable_length_records(j).record_length_after_header))); % zero padding % TODO
                        fwrite(fid, string(1:s.extended_variable_length_records(j).record_length_after_header), 'char', 0, MACHINE_FORMAT);
                        fprintf('WARNING: writing Text area description (optional) EVLR\n');
                        
                    case 4 % Extra bytes (optional)
                        
                        n_extra_records = s.extended_variable_length_records(j).record_length_after_header / 192;
                        
                        for k = 1:n_extra_records
                            
                            fwrite(fid, s.extended_variable_length_records(j).value(k).reserved, 'uint16', 0, MACHINE_FORMAT);
                            fwrite(fid, s.extended_variable_length_records(j).value(k).data_type, 'uint8', 0, MACHINE_FORMAT);
                            fwrite(fid, s.extended_variable_length_records(j).value(k).options, 'uint8', 0, MACHINE_FORMAT);
                            name = horzcat(s.extended_variable_length_records(j).value(k).name, char(zeros(1, 32))); % zero padding
                            fwrite(fid, name(1:32), 'char', 0, MACHINE_FORMAT);
                            
                            fwrite(fid, s.extended_variable_length_records(j).value(k).unused, 'uint32', 0, MACHINE_FORMAT);
                            fwrite(fid, s.extended_variable_length_records(j).value(k).no_data, 'double', 0, MACHINE_FORMAT);
                            fwrite(fid, s.extended_variable_length_records(j).value(k).min, 'double', 0, MACHINE_FORMAT);
                            fwrite(fid, s.extended_variable_length_records(j).value(k).max, 'double', 0, MACHINE_FORMAT);
                            fwrite(fid, s.extended_variable_length_records(j).value(k).scale, 'double', 0, MACHINE_FORMAT);
                            fwrite(fid, s.extended_variable_length_records(j).value(k).offset, 'double', 0, MACHINE_FORMAT);
                            description = horzcat(s.extended_variable_length_records(j).value(k).description, char(zeros(1, 32))); % zero padding
                            fwrite(fid, description(1:32), 'char', 0, MACHINE_FORMAT);
                            
                        end
                        
                        fprintf('WARNING: writing Extra bytes (optional) EVLR\n');
                        
                    case num2cell(100:354) % Waveform Packet descriptor (required when using point formats 4, 5, 9, 10)
                        
                        fwrite(fid, s.extended_variable_length_records(j).value.bits_per_sample, 'uint8', 0, MACHINE_FORMAT); % Bits per sample, Unsigned char, 1 byte, *
                        fwrite(fid, s.extended_variable_length_records(j).value.compression_type, 'uint8', 0, MACHINE_FORMAT); % Waveform compression type, Unsigned char, 1 byte, *
                        fwrite(fid, s.extended_variable_length_records(j).value.number_samples, 'uint32', 0, MACHINE_FORMAT); % Number of samples, Unsigned long, 4 bytes, *
                        fwrite(fid, s.extended_variable_length_records(j).value.temporal_sample_spacing, 'uint32', 0, MACHINE_FORMAT); % Temporal Sample Spacing, Unsigned long, 4 bytes, *
                        fwrite(fid, s.extended_variable_length_records(j).value.digitizer_gain, 'double', 0, MACHINE_FORMAT); % Digitizer Gain, double, 8 bytes, *
                        fwrite(fid, s.extended_variable_length_records(j).value.digitizer_offset, 'double', 0, MACHINE_FORMAT); % Digitizer Offset, double, 8 bytes, *
                        fprintf('WARNING: writing Waveform Packet descriptor EVLR\n');
                        
                    case 2111 % OGC Math Transform WKT Record (optional)
                        
                        string = horzcat(s.extended_variable_length_records(j).value, char(zeros(1, s.extended_variable_length_records(j).record_length_after_header))); % zero padding % TODO
                        fwrite(fid, string(1:s.extended_variable_length_records(j).record_length_after_header), 'char', 0, MACHINE_FORMAT);
                        fprintf('WARNING: writing OGC Math Transform WKT Record (optional) EVLR\n');
                        
                    case 2112 % OGC Coordinate System WKT Record (optional)
                        
                        string = horzcat(s.extended_variable_length_records(j).value, char(zeros(1, s.extended_variable_length_records(j).record_length_after_header))); % zero padding % TODO
                        fwrite(fid, string(1:s.extended_variable_length_records(j).record_length_after_header), 'char', 0, MACHINE_FORMAT);
                        fprintf('WARNING: writing OGC Coordinate System WKT Record (optional) EVLR\n');
        
                    case 34735 % GeoKeyDirectoryTag Record (mandatory)
                        
                        s.extended_variable_length_records(j).value.w_key_directory_version = 1; % wKeyDirectoryVersion, unsigned short, wKeyDirectoryVersion = 1 Always
                        s.extended_variable_length_records(j).value.w_key_revision = 1; % wKeyRevision, unsigned short, wKeyRevision = 1 Always
                        s.extended_variable_length_records(j).value.w_minor_revision = 0; % wMinorRevision, unsigned short, wMinorRevision = 0 Always
                        s.extended_variable_length_records(j).value.w_number_of_keys = length(s.extended_variable_length_records(j).value.s_key_entry); % wNumberOfKeys, unsigned short, Number of sets of 4 unsigned shorts to follow
                        
                        fwrite(fid, s.extended_variable_length_records(j).value.w_key_directory_version, 'uint16', 0, MACHINE_FORMAT);
                        fwrite(fid, s.extended_variable_length_records(j).value.w_key_revision, 'uint16', 0, MACHINE_FORMAT);
                        fwrite(fid, s.extended_variable_length_records(j).value.w_minor_revision, 'uint16', 0, MACHINE_FORMAT);
                        fwrite(fid, s.extended_variable_length_records(j).value.w_number_of_keys, 'uint16', 0, MACHINE_FORMAT);
                        
                        for i=1:s.extended_variable_length_records(j).value.w_number_of_keys
                            
                            fwrite(fid, s.extended_variable_length_records(j).value.s_key_entry{i}(1), 'uint16', 0, MACHINE_FORMAT);
                            fwrite(fid, s.extended_variable_length_records(j).value.s_key_entry{i}(2), 'uint16', 0, MACHINE_FORMAT);
                            fwrite(fid, s.extended_variable_length_records(j).value.s_key_entry{i}(3), 'uint16', 0, MACHINE_FORMAT);
                            fwrite(fid, s.extended_variable_length_records(j).value.s_key_entry{i}(4), 'uint16', 0, MACHINE_FORMAT);
                            
                        end
                        
                        fprintf('WARNING: writing GeoKeyDirectoryTag Record to EVLR\n');
                        
                    case 34736 % GeoDoubleParamsTag Record (optional)
                        
                        fwrite(fid, s.extended_variable_length_records(j).value, 'double', 0, MACHINE_FORMAT);
                        fprintf('WARNING: writing GeoDoubleParamsTag Record to EVLR\n');
                        
                    case 34737 % GeoAsciiParamsTag Record (optional)
                        
                        string = horzcat(s.extended_variable_length_records(j).value, char(zeros(1, s.extended_variable_length_records(j).record_length_after_header))); % zero padding % TODO
                        fwrite(fid, string(1:s.extended_variable_length_records(j).record_length_after_header), 'char', 0, MACHINE_FORMAT);
                        fprintf('WARNING: writing GeoAsciiParamsTag Record to EVLR\n');
                        
                    case {5013, 5017, 5018, 5043, 5044, 5061, 5063, 5071} % Custom uint8 field (optional)
                        
                        fwrite(fid, s.extended_variable_length_records(j).value, 'uint8', 0, MACHINE_FORMAT);
                        fprintf('WARNING: writing Custom uint8 field to EVLR\n');
                        
                    case {5002, 5014, 5015, 5016} % Custom uint16 field (optional)
                        
                        fwrite(fid, s.extended_variable_length_records(j).value, 'uint16', 0, MACHINE_FORMAT);
                        fprintf('WARNING: writing Custom uint16 field to EVLR\n');
                    
                    case 5001 % Custom char(32) field (optional) % TEST
                        
                        char32 = char(cellfun(@(x) horzcat(x, char(zeros(1, 32-length(x)))), s.extended_variable_length_records(j).value, 'UniformOutput', false)); % TODO
                        fwrite(fid, char32, 'char', 0, MACHINE_FORMAT);
                        fprintf('WARNING: writing Custom char(32) field to EVLR\n');
                        
                    case {5000, 5019, 5020, 5021} % Custom single field (optional)
                        
                        fwrite(fid, s.extended_variable_length_records(j).value, 'single', 0, MACHINE_FORMAT);
                        fprintf('WARNING: writing Custom single field to EVLR\n');
                        
                    case {5010, 5011, 5012, 5050, 5062} % Custom double field (optional)
                        
                        fwrite(fid, s.extended_variable_length_records(j).value, 'double', 0, MACHINE_FORMAT);
                        fprintf('WARNING: writing Custom double field to EVLR\n');
                        
                    case {5041} % Custom char(12) field (optional)
                
                        char12 = char(cellfun(@(x) horzcat(x, char(zeros(1, 12-length(x)))), s.extended_variable_length_records(j).value, 'UniformOutput', false)); % TODO
                        fwrite(fid, char12, 'char', 0, MACHINE_FORMAT);
                        fprintf('WARNING: writing Custom char(12) field to EVLR\n');
                        
                        % You may add custom records here
                        
                        
                    otherwise % Other
                        
                        fwrite(fid, typecast(s.extended_variable_length_records(j).value, 'uint8'), 'uint8', 0, MACHINE_FORMAT);
                        fprintf('WARNING: writing Custom uint8 field to EVLR\n');
                        
                end

                byte_offset = byte_offset + 60 + s.extended_variable_length_records(j).record_length_after_header;
                
            end
            
        end

    
    %% close file
    
    fclose(fid);
    
    if arg.Results.verbose
        
        tElapsed = toc; % stop timer
        
        fprintf('File successfully written to "%s"...\n', arg.Results.filepath);
        fprintf('%u point records written in %s\n', r.header(phb_skeys.n_point_records).value, datestr(tElapsed/(24*3600), 'HH:MM:SS.FFF'));
        
    end
    
end

varargout{1} = r;
