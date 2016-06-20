function s = LASread(filepath, varargin)
%LASREAD - Reads the ASPRS LAS point cloud format (LAS versions 1.0 to 1.4 are supported)
% S = LASREAD(FILEPATH, HEADERONLY, VERBOSE) reads the LAS file specified in FILEPATH
% into the structure S.
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
% http://www.asprs.org/Committee-General/LASer-LAS-File-Format-Exchange-Activitier.html
%
% Syntax:  [s] = LASread(filepath, headerOnly, verbose)
%
% Inputs:
%    filepath - The path to the input LAS file
%    headerOnly (optional, default: false) - If set to true, will read only the public header content
%    verbose (optional, default: true) - If set to true, will display the public header information
%              in the console
%
% Outputs:
%    s - A structure containing the fixed and variable length headers and point data records
%
% Example:
%    s = LASread('E:\data\points.las', false, true)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2015b, GNU Octave 4.0.0 (configured for "i686-w64-mingw32")
%
% See also: LASWRITE
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory
% Website: http://lasig.epfl.ch/
% Last revision: June 20, 2016
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund, WHFF (OFEV) - project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details

fclose('all'); % close any open files

%% setup constants

MACHINE_FORMAT = 'ieee-le'; % all data is in little-endian format
OCTAVE_FLAG = (exist('OCTAVE_VERSION', 'builtin') ~= 0); % determine if system is Matlab or GNU Octave
AUTH_RECORD_FORMAT_ID = 0:10; % authorized point record format IDs
AUTH_LAS_MAJOR_VERSIONS = [1]; % authorized LAS major versions
AUTH_LAS_MINOR_VERSIONS = [0, 1, 2, 3, 4]; % authorized LAS minor versions
AUTH_LAS_VERSIONS = [10, 11, 12, 13, 14]; % authorized LAS versions
PUBLIC_HEADER_SIZES = [227, 227, 227, 235, 375]; % public header sizes for LAS 1.0-1.4
STANDARD_RECORD_LENGTHS = [20, 28, 26, 34, 57, 63, 30, 36, 38, 59, 67]; % standard point data record lengths


%% check argument validity

arg = inputParser;

addRequired(arg, 'filepath', @ischar);
addOptional(arg, 'headerOnly', false, @(x) islogical(x) && (numel(x) == 1));
addOptional(arg, 'verbose', true, @(x) islogical(x) && (numel(x) == 1));

parse(arg, filepath, varargin{:});


%% check if file is accessible and in correct format

fid = fopen(filepath, 'r');
fseek(fid, 0, 'eof');
filesize = ftell(fid);
fseek(fid, 0, 'bof');

if fid == -1
    
    error('Could not open file')
    
end

if ~strcmp(sprintf('%c', fread(fid, 4, 'char', 0, MACHINE_FORMAT)), 'LASF')
    
    error('Input file is not a valid LAS file')
    
end


%% determine LAS version

% read LAS major version
fseek(fid, 24, 'bof');
las_version_major = fread(fid, 1, 'uint8', 0, MACHINE_FORMAT);

% read LAS minor version
fseek(fid, 25, 'bof');
las_version_minor = fread(fid, 1, 'uint8', 0, MACHINE_FORMAT);
las_version = las_version_major*10 + las_version_minor;


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
r.header(k).value = [];
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
r.header(k).value = [];
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
r.header(k).value = [];
r.header(k).validation = []; % @(x) (x == 0) && (numel(x) == 1);
r.header(k).error_id = [];
r.header(k).error_message = [];
r.header(k).warning_message = 'The Reserved field (LAS 1.1 only) was not set to 0';
k = k + 1;

% Global Encoding - GPS Time Type, bit 0, 1 bit, *
% The meaning of GPS Time in the Point Records
% 0 (not set) -> GPS time in the point record fields is GPS Week Time (the same as previous versions of LAS)
% 1 (set) -> GPS Time is standard GPS Time (satellite GPS Time) minus 1 x 10^9.
% The offset moves the time back to near zero to improve floating point resolution.
r.header(k).compatibility = [12, 13, 14];
r.header(k).full_name = 'Global Encoding - GPS Time Type';
r.header(k).short_name = 'global_encoding_gps_time_type';
r.header(k).type = {'ubit1', 'ubit1', 'ubit1'};
r.header(k).storage_type = {'logical', 'logical', 'logical'};
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
r.header(k).storage_type = {'double', 'double', 'double'};
r.header(k).byte_length = [15/8, 12/8, 11/8];
r.header(k).bit_position = [1, 4, 5];
r.header(k).n_values = [1, 1, 1];
r.header(k).flag_bit_field = [true, true, true];
r.header(k).print_format = {'%u', '%u', '%u'};
r.header(k).default_value = 0;
r.header(k).value = [];
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
r.header(k).default_value = [];
r.header(k).value = [];
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
r.header(k).default_value = [];
r.header(k).value = [];
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
r.header(k).default_value = [];
r.header(k).value = [];
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
r.header(k).default_value = [];
r.header(k).value = [];
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
r.header(k).default_value = [];
r.header(k).value = [];
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
r.header(k).default_value = [];
r.header(k).value = [];
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
r.header(k).default_value = [];
r.header(k).value = [];
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
r.header(k).default_value = [];
r.header(k).value = [];
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
r.header(k).default_value = [];
r.header(k).value = [];
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
r.header(k).default_value = [];
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
r.header(k).default_value = [];
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
r.header(k).default_value = [];
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
r.header(k).default_value = [];
r.header(k).value = [];
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
r.header(k).default_value = [];
r.header(k).value = [];
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
r.header(k).default_value = [];
r.header(k).value = [];
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
r.header(k).default_value = [];
r.header(k).value = [];
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
r.header(k).default_value = [];
r.header(k).value = [];
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
r.header(k).default_value = [];
r.header(k).value = [];
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
r.header(k).validation = @(x) isa(x, 'uint64') && (numel(x) == 1); % x <= intmax('uint64')
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
r.header(k).default_value = [];
r.header(k).value = [];
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
r.record(k).format_compatibility = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
r.record(k).bit_position = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32'};
r.record(k).storage_type = {'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32'};
r.record(k).byte_length = [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4];
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
r.record(k).format_compatibility = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
r.record(k).bit_position = [32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32'};
r.record(k).storage_type = {'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32'};
r.record(k).byte_length = [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4];
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
r.record(k).format_compatibility = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
r.record(k).bit_position = [64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32'};
r.record(k).storage_type = {'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32', 'int32'};
r.record(k).byte_length = [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4];
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
r.record(k).format_compatibility = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
r.record(k).bit_position = [96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16'};
r.record(k).storage_type = {'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16'};
r.record(k).byte_length = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = false;
r.record(k).flag_transform = false;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% Return Number, 3 bits (bits 0, 1, 2), 3 bits, *
r.record(k).full_name = 'Return Number';
r.record(k).short_name = 'return_number';
r.record(k).format_compatibility = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
r.record(k).bit_position = [112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'ubit3', 'ubit3', 'ubit3', 'ubit3', 'ubit3', 'ubit3', 'ubit4', 'ubit4', 'ubit4', 'ubit4', 'ubit4'};
r.record(k).storage_type = {'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8'};
r.record(k).byte_length = [3/8, 3/8, 3/8, 3/8, 3/8, 3/8, 4/8, 4/8, 4/8, 4/8, 4/8];
r.record(k).flag_bit_field = true;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% Number of Returns (given pulse), 3 bits (bits 3, 4, 5), 3 bits, *
r.record(k).full_name = 'Number of Returns (given pulse)';
r.record(k).short_name = 'number_of_returns';
r.record(k).format_compatibility = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
r.record(k).bit_position = [115, 115, 115, 115, 115, 115, 116, 116, 116, 116, 116];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'ubit3', 'ubit3', 'ubit3', 'ubit3', 'ubit3', 'ubit3', 'ubit4', 'ubit4', 'ubit4', 'ubit4', 'ubit4'};
r.record(k).storage_type = {'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8'};
r.record(k).byte_length = [3/8, 3/8, 3/8, 3/8, 3/8, 3/8, 4/8, 4/8, 4/8, 4/8, 4/8];
r.record(k).flag_bit_field = true;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% Classification flag 1 - Synthetic, 1 bit (bit 0), 1 bit, *
r.record(k).full_name = 'Classification flag 1 - Synthetic';
r.record(k).short_name = 'classification_flag_1';
r.record(k).format_compatibility = [6, 7, 8, 9, 10];
r.record(k).bit_position = [120, 120, 120, 120, 120];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1'};
r.record(k).storage_type = {'uint8', 'uint8', 'uint8', 'uint8', 'uint8'};
r.record(k).byte_length = [1/8, 1/8, 1/8, 1/8, 1/8];
r.record(k).flag_bit_field = true;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% Classification flag 2 - Key-point, 1 bit (bit 1), 1 bit, *
r.record(k).full_name = 'Classification flag 2 - Key-point';
r.record(k).short_name = 'classification_flag_2';
r.record(k).format_compatibility = [6, 7, 8, 9, 10];
r.record(k).bit_position = [121, 121, 121, 121, 121];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1'};
r.record(k).storage_type = {'uint8', 'uint8', 'uint8', 'uint8', 'uint8'};
r.record(k).byte_length = [1/8, 1/8, 1/8, 1/8, 1/8];
r.record(k).flag_bit_field = true;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% Classification flag 3 - Withheld, 1 bit (bit 2), 1 bit, *
r.record(k).full_name = 'Classification flag 3 - Withheld';
r.record(k).short_name = 'classification_flag_3';
r.record(k).format_compatibility = [6, 7, 8, 9, 10];
r.record(k).bit_position = [122, 122, 122, 122, 122];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1'};
r.record(k).storage_type = {'uint8', 'uint8', 'uint8', 'uint8', 'uint8'};
r.record(k).byte_length = [1/8, 1/8, 1/8, 1/8, 1/8];
r.record(k).flag_bit_field = true;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% Classification flag 4 - Overlap, 1 bit (bit 3), 1 bit, *
r.record(k).full_name = 'Classification flag 4 - Overlap';
r.record(k).short_name = 'classification_flag_4';
r.record(k).format_compatibility = [6, 7, 8, 9, 10];
r.record(k).bit_position = [123, 123, 123, 123, 123];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1'};
r.record(k).storage_type = {'uint8', 'uint8', 'uint8', 'uint8', 'uint8'};
r.record(k).byte_length = [1/8, 1/8, 1/8, 1/8, 1/8];
r.record(k).flag_bit_field = true;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% Scanner channel, 2 bits (bit 4-5), 1 bit, *
r.record(k).full_name = 'Scanner channel';
r.record(k).short_name = 'scanner_channel';
r.record(k).format_compatibility = [6, 7, 8, 9, 10];
r.record(k).bit_position = [124, 124, 124, 124, 124];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'ubit2', 'ubit2', 'ubit2', 'ubit2', 'ubit2'};
r.record(k).storage_type = {'uint8', 'uint8', 'uint8', 'uint8', 'uint8'};
r.record(k).byte_length = [2/8, 2/8, 2/8, 2/8, 2/8];
r.record(k).flag_bit_field = true;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% Scan Direction Flag, 1 bit (bit 6), 1 bit, *
r.record(k).full_name = 'Scan Direction Flag';
r.record(k).short_name = 'scan_direction_flag';
r.record(k).format_compatibility = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
r.record(k).bit_position = [118, 118, 118, 118, 118, 118, 126, 126, 126, 126, 126];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1'};
r.record(k).storage_type = {'logical', 'logical', 'logical', 'logical', 'logical', 'logical', 'logical', 'logical', 'logical', 'logical', 'logical'};
r.record(k).byte_length = [1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8];
r.record(k).flag_bit_field = true;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% Edge of Flight Line, 1 bit (bit 7), 1 bit, *
r.record(k).full_name = 'Edge of Flight Line';
r.record(k).short_name = 'flightline_edge_flag';
r.record(k).format_compatibility = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
r.record(k).bit_position = [119, 119, 119, 119, 119, 127, 127, 127, 127, 127];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1'};
r.record(k).storage_type = {'logical', 'logical', 'logical', 'logical', 'logical', 'logical', 'logical', 'logical', 'logical', 'logical', 'logical'};
r.record(k).byte_length = [1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8];
r.record(k).flag_bit_field = true;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% ASPRS classification number (0-31), 5 bits (bits 0, 1, 2, 3, 4), *
r.record(k).full_name = 'Classification (number)';
r.record(k).short_name = 'classification'; % 'classification_number';
r.record(k).format_compatibility = [0, 1, 2, 3, 4, 5];
r.record(k).bit_position = [120, 120, 120, 120, 120, 120];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'ubit5', 'ubit5', 'ubit5', 'ubit5', 'ubit5', 'ubit5'};
r.record(k).storage_type = {'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8'};
r.record(k).byte_length = [5/8, 5/8, 5/8, 5/8, 5/8, 5/8];
r.record(k).flag_bit_field = true;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% Synthetic - If set, then this point was created by a technique other than LIDAR collection such as digitized from a photogrammetric stereo model, 1 bit, *
r.record(k).full_name = 'Classification (synthetic flag)';
r.record(k).short_name = 'classification_synthetic';
r.record(k).format_compatibility = [0, 1, 2, 3, 4, 5];
r.record(k).bit_position = [125, 125, 125, 125, 125, 125];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1'};
r.record(k).storage_type = {'logical', 'logical', 'logical', 'logical', 'logical', 'logical'};
r.record(k).byte_length = [1/8, 1/8, 1/8, 1/8, 1/8, 1/8];
r.record(k).flag_bit_field = true;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% Keypoint - If set, this point is considered to be a model key-point and thus generally should not be withheld in a thining algorithm, 1 bit, *
r.record(k).full_name = 'Classification (key-point flag)';
r.record(k).short_name = 'classification_keypoint';
r.record(k).format_compatibility = [0, 1, 2, 3, 4, 5];
r.record(k).bit_position = [126, 126, 126, 126, 126, 126];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1'};
r.record(k).storage_type = {'logical', 'logical', 'logical', 'logical', 'logical', 'logical'};
r.record(k).byte_length = [1/8, 1/8, 1/8, 1/8, 1/8, 1/8];
r.record(k).flag_bit_field = true;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% Withheld - If set, this point should not be included in processing (synonymous with Deleted), 1 bit, *
r.record(k).full_name = 'Classification (withheld flag)';
r.record(k).short_name = 'classification_withheld';
r.record(k).format_compatibility = [0, 1, 2, 3, 4, 5];
r.record(k).bit_position = [127, 127, 127, 127, 127, 127];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1', 'ubit1'};
r.record(k).storage_type = {'logical', 'logical', 'logical', 'logical', 'logical', 'logical'};
r.record(k).byte_length = [1/8, 1/8, 1/8, 1/8, 1/8, 1/8];
r.record(k).flag_bit_field = true;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% Scan Angle, unsigned char, 1 byte, *
r.record(k).full_name = 'Scan Angle';
r.record(k).short_name = 'scan_angle';
r.record(k).format_compatibility = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
r.record(k).bit_position = [128, 128, 128, 128, 128, 128, 144, 144, 144, 144, 144];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'int16', 'int16', 'int16', 'int16', 'int16'};
r.record(k).storage_type = {'uchar', 'uchar', 'uchar', 'uchar', 'uchar', 'uchar', 'int16', 'int16', 'int16', 'int16', 'int16'};
r.record(k).byte_length = [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2];
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% Classification, unsigned char, 1 byte, *
r.record(k).full_name = 'Classification';
r.record(k).short_name = 'classification';
r.record(k).format_compatibility = [6, 7, 8, 9, 10];
r.record(k).bit_position = [128, 128, 128, 128, 128];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'uint8', 'uint8', 'uint8', 'uint8', 'uint8'};
r.record(k).storage_type = {'uint8', 'uint8', 'uint8', 'uint8', 'uint8'};
r.record(k).byte_length = [1, 1, 1, 1, 1];
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% User Data, unsigned char, 1 byte
r.record(k).full_name = 'User Data';
r.record(k).short_name = 'user_data';
r.record(k).format_compatibility = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
r.record(k).bit_position = [136, 136, 136, 136, 136, 136, 136, 136, 136, 136, 136];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8'};
r.record(k).storage_type = {'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8', 'uint8'};
r.record(k).byte_length = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = false;
r.record(k).flag_transform = false;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% Point Source ID, unsigned short, 2 bytes, *
r.record(k).full_name = 'Point Source ID';
r.record(k).short_name = 'point_source_id';
r.record(k).format_compatibility = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
r.record(k).bit_position = [144, 144, 144, 144, 144, 144, 160, 160, 160, 160, 160];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16'};
r.record(k).storage_type = {'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16'};
r.record(k).byte_length = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% GPS Time, double, 8 bytes, *
r.record(k).full_name = 'GPS Time';
r.record(k).short_name = 'gps_time';
r.record(k).format_compatibility = [1, 3, 4, 5, 6, 7, 8, 9, 10];
r.record(k).bit_position = [160, 160, 160, 160, 176, 176, 176, 176, 176];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'uint64', 'uint64', 'uint64', 'uint64', 'uint64', 'uint64', 'uint64', 'uint64', 'uint64'};
r.record(k).storage_type = {'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'};
r.record(k).byte_length = [8, 8, 8, 8, 8, 8, 8, 8, 8];
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% Red, unsigned short, 2 bytes, *
r.record(k).full_name = 'Red';
r.record(k).short_name = 'red';
r.record(k).format_compatibility = [2, 3, 5, 7, 8, 10];
r.record(k).bit_position = [160, 224, 224, 240, 240, 240];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16'};
r.record(k).storage_type = {'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16'};
r.record(k).byte_length = [2, 2, 2, 2, 2, 2];
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% Green, unsigned short, 2 bytes, *
r.record(k).full_name = 'Green';
r.record(k).short_name = 'green';
r.record(k).format_compatibility = [2, 3, 5, 7, 8, 10];
r.record(k).bit_position = [176, 240, 240, 256, 256, 256];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16'};
r.record(k).storage_type = {'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16'};
r.record(k).byte_length = [2, 2, 2, 2, 2, 2];
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% Blue, unsigned short, 2 bytes, *
r.record(k).full_name = 'Blue';
r.record(k).short_name = 'blue';
r.record(k).format_compatibility = [2, 3, 5, 7, 8, 10];
r.record(k).bit_position = [192, 256, 256, 272, 272, 272];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16'};
r.record(k).storage_type = {'uint16', 'uint16', 'uint16', 'uint16', 'uint16', 'uint16'};
r.record(k).byte_length = [2, 2, 2, 2, 2, 2];
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% NIR, unsigned short, 2 bytes, *
r.record(k).full_name = 'Near infrared';
r.record(k).short_name = 'nir';
r.record(k).format_compatibility = [8, 10];
r.record(k).bit_position = [288, 288];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'uint16', 'uint16'};
r.record(k).storage_type = {'uint16', 'uint16'};
r.record(k).byte_length = [2, 2];
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% Wave Packet Descriptor Index, Unsigned char, 1 byte, *
r.record(k).full_name = 'Wave Packet Descriptor Index';
r.record(k).short_name = 'wave_packet_descriptor_index';
r.record(k).format_compatibility = [4, 5, 9, 10];
r.record(k).bit_position = [224, 272, 240, 304];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'uint8', 'uint8', 'uint8', 'uint8'};
r.record(k).storage_type = {'uint8', 'uint8', 'uint8', 'uint8'};
r.record(k).byte_length = [1, 1, 1, 1];
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% Byte offset to waveform data, Unsigned long long, 8 bytes, *
r.record(k).full_name = 'Byte offset to waveform data';
r.record(k).short_name = 'offset_to_waveform_data';
r.record(k).format_compatibility = [4, 5, 9, 10];
r.record(k).bit_position = [232, 280, 248, 312];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'uint64', 'uint64', 'uint64', 'uint64'};
r.record(k).storage_type = {'uint64', 'uint64', 'uint64', 'uint64'};
r.record(k).byte_length = [8, 8, 8, 8];
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% Waveform packet size in bytes, Unsigned long, 4 bytes, *
r.record(k).full_name = 'Waveform packet size in bytes';
r.record(k).short_name = 'waveform_packet_size';
r.record(k).format_compatibility = [4, 5, 9, 10];
r.record(k).bit_position = [296, 344, 312, 376];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'uint32', 'uint32', 'uint32', 'uint32'};
r.record(k).storage_type = {'uint32', 'uint32', 'uint32', 'uint32'};
r.record(k).byte_length = [4, 4, 4, 4];
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% Return Point Waveform Location, float, 4 bytes, *
r.record(k).full_name = 'Return Point Waveform Location';
r.record(k).short_name = 'return_point_waveform_location';
r.record(k).format_compatibility = [4, 5, 9, 10];
r.record(k).bit_position = [328, 376, 344, 408];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'uint32', 'uint32', 'uint32', 'uint32'};
r.record(k).storage_type = {'uint32', 'uint32', 'uint32', 'uint32'};
r.record(k).byte_length = [4, 4, 4, 4];
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% X(t), float, 4 bytes, *
r.record(k).full_name = 'X(t)';
r.record(k).short_name = 'x_t';
r.record(k).format_compatibility = [4, 5, 9, 10];
r.record(k).bit_position = [360, 408, 376, 440];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'uint32', 'uint32', 'uint32', 'uint32'};
r.record(k).storage_type = {'single', 'single', 'single', 'single'};
r.record(k).byte_length = [4, 4, 4, 4];
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% Y(t), float, 4 bytes, *
r.record(k).full_name = 'Y(t)';
r.record(k).short_name = 'y_t';
r.record(k).format_compatibility = [4, 5, 9, 10];
r.record(k).bit_position = [392, 440, 408, 472];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'uint32', 'uint32', 'uint32', 'uint32'};
r.record(k).storage_type = {'single', 'single', 'single', 'single'};
r.record(k).byte_length = [4, 4, 4, 4];
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];
k = k + 1;

% Z(t), float, 4 bytes, *
r.record(k).full_name = 'Z(t)';
r.record(k).short_name = 'z_t';
r.record(k).format_compatibility = [4, 5, 9, 10];
r.record(k).bit_position = [424, 472, 440, 504];
r.record(k).byte_position = r.record(k).bit_position ./ 8;
r.record(k).type = {'uint32', 'uint32', 'uint32', 'uint32'};
r.record(k).storage_type = {'single', 'single', 'single', 'single'};
r.record(k).byte_length = [4, 4, 4, 4];
r.record(k).flag_bit_field = false;
r.record(k).flag_mandatory = true;
r.record(k).flag_transform = false;
r.record(k).scale = [];
r.record(k).offset = [];
r.record(k).value = [];


%% extra bytes format definition

k = 1;
r.extra_bytes(k).id = 0;
r.extra_bytes(k).type = 'NaN';
r.extra_bytes(k).n_values = 1;
r.extra_bytes(k).byte_length = 'NaN';
k = k + 1;

r.extra_bytes(k).id = 1;
r.extra_bytes(k).type = 'uchar';
r.extra_bytes(k).n_values = 1;
r.extra_bytes(k).byte_length = 1;
k = k + 1;

r.extra_bytes(k).id = 2;
r.extra_bytes(k).type = 'char';
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
r.extra_bytes(k).type = 'uint64';
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


%% adjust public header based on the LAS version

% filter version specific attributes
idx_header_format = false(length(r.header),1);

if ismember(las_version, AUTH_LAS_VERSIONS)
    
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
            
            end
            
        end
        
    end
    
else
    
    error('LASread:wrongPointDataFormatID', 'The point data format is not in the valid range (0-10)');
    
end

r.header(~idx_header_format) = [];


%% generate key-value pairs for public header block

for j = 1:length(r.header)
    
    phb_skeys.(r.header(j).short_name) = j;
    
end


%% read public header block

% read public header size
fseek(fid, 94, 'bof');
r.header(phb_skeys.header_size).value = fread(fid, 1, 'uint16', 0, MACHINE_FORMAT);

byte_offset = 0;
integer_byte_offset = 0;

% read public header content
for j = 1:length(r.header)
    
    fseek(fid, integer_byte_offset, 'bof');
    
    switch r.header(j).type
        
        case {'uchar','char'}
            
            r.header(j).value = regexp(fread(fid, r.header(j).n_values, '*char', 0, MACHINE_FORMAT)', '[^\0]*', 'match', 'once'); % read until null terminator
            
        case {'uint8', 'uint16', 'uint32', 'uint64', 'double' }
            
            r.header(j).value = fread(fid, r.header(j).n_values, r.header(j).type, 0, MACHINE_FORMAT);
            
        otherwise
            
            % move to previous byte
            fseek(fid, integer_byte_offset-1, 'bof');
            
            % skip offset bits
            fread(fid, 1, 'ubit8', r.header(j).bit_position, MACHINE_FORMAT);
            
            % read bit(s)
            r.header(j).value = fread(fid, r.header(j).n_values, r.header(j).type, 0, MACHINE_FORMAT);
            
    end
    
    byte_offset = byte_offset + r.header(j).n_values * r.header(j).byte_length;
    
    if mod(byte_offset,1) == 0
        
        integer_byte_offset = floor(byte_offset);
        
    end
    
    if arg.Results.verbose
        
        fprintf(sprintf('%s: %s\\n', r.header(j).full_name, r.header(j).print_format), r.header(j).value);
        
    end
    
end

% For LAS >= 1.4 update "number of point records" and "number of points by
% return" so they match the values in the extended fields
if (las_version >= 14) && ~isempty(r.header(phb_skeys.n_point_records_extended).value)
    
    r.header(phb_skeys.n_point_records).value = r.header(phb_skeys.n_point_records_extended).value;
    r.header(phb_skeys.n_points_by_return).value = r.header(phb_skeys.n_points_by_return_extended).value;
    
end

% convert GUID data entries to hexadecimal
% A GUID is a 128-bit value consisting of one group of 8 hexadecimal digits,
% followed by three groups of 4 hexadecimal digits each, followed by one group of 12 hexadecimal digits
% 4 bytes | 2 bytes  | 2 bytes  | 2 bytes |  6 bytes
% r.header(phb_skeys.project_id_1).value = dec2hex(r.header(phb_skeys.project_id_1).value);
% r.header(phb_skeys.project_id_2).value = dec2hex(r.header(phb_skeys.project_id_2).value);
% r.header(phb_skeys.project_id_3).value = dec2hex(r.header(phb_skeys.project_id_3).value);
% r.header(phb_skeys.project_id_4).value = dec2hex(r.header(phb_skeys.project_id_4).value);


%% adjust point record based on the Point Data Format ID (0-99 for spec)

idx_record_format = false(length(r.record),1);

[flag_record_format, ~] = ismember(r.header(phb_skeys.point_data_format_id).value, AUTH_RECORD_FORMAT_ID);

if flag_record_format
    
    for j = 1:length(r.record)
        
        [flag_include_record, ind_record] = ismember(r.header(phb_skeys.point_data_format_id).value, r.record(j).format_compatibility);
        
        if flag_include_record
            
            idx_record_format(j) = true;
            
            r.record(j).format_compatibility = r.record(j).format_compatibility(ind_record);
            r.record(j).bit_position = r.record(j).bit_position(ind_record);
            r.record(j).byte_position = r.record(j).byte_position(ind_record);
            r.record(j).type = r.record(j).type(ind_record);
            r.record(j).byte_length = r.record(j).byte_length(ind_record);
            r.record(j).storage_type = r.record(j).storage_type(ind_record);
            
        end
        
    end
    
else
    
    error('LASread:wrongPointDataFormatID', 'The point data format is not in the valid range (0-10)');
    
end

r.record(~idx_record_format) = [];

[~, ind_record_sort] = sort([r.record.bit_position]);
ind_last_record = ind_record_sort(end);
ind_record = length(r.record); % index of last record


%% read variable length records

flag_vlr = (r.header(phb_skeys.n_variable_length_records).value > 0);

if flag_vlr
    
    byte_offset = r.header(phb_skeys.header_size).value;
    
    for j = 1:r.header(phb_skeys.n_variable_length_records).value
        
        fseek(fid, byte_offset, 'bof');
        
        r.variable_length_records(j).reserved = fread(fid, 1, 'uint16', 0, MACHINE_FORMAT); % Reserved, unsigned short, 2 bytes;
        r.variable_length_records(j).user_id = regexp(fread(fid, 16, '*char', 0, MACHINE_FORMAT)', '[^\0]*', 'match', 'once'); % read until null terminator
        r.variable_length_records(j).record_id = fread(fid, 1, 'uint16', 0, MACHINE_FORMAT); % Record ID, unsigned short, 2 bytes, *
        r.variable_length_records(j).record_length_after_header = fread(fid, 1, 'uint16', 0, MACHINE_FORMAT); % Record Length After Header, unsigned short, 2 bytes, *
        r.variable_length_records(j).description = regexp(fread(fid, 32, '*char', 0, MACHINE_FORMAT)', '[^\0]*', 'match', 'once'); % Description, char[32], 32 bytes
        
        % read any defined variable length records
        switch r.variable_length_records(j).record_id
            
            case 0 % Classification lookup (optional)
                
                r.variable_length_records(j).value.class_number = fread(fid, r.variable_length_records(j).record_length_after_header, '*uint8', 0, MACHINE_FORMAT);
                r.variable_length_records(j).value.class_description = sprintf('%c', fread(fid, r.variable_length_records(j).record_length_after_header, 'char', 0, MACHINE_FORMAT));
                
            case 1 % Reserved (optional)
                
                
            case 2 % Histogram (optional)
                
                
            case 3 % Text area description (optional)
                
                r.variable_length_records(j).value = regexp(fread(fid, r.variable_length_records(j).record_length_after_header, '*char', 0, MACHINE_FORMAT)', '[^\0]*', 'match', 'once');
                
            case 4 % Extra bytes (optional)
                
                %flag_extra_records = true;
                n_extra_records = r.variable_length_records(j).record_length_after_header / 192;
                
                for k = 1:n_extra_records
                    
                    r.variable_length_records(j).value(k).reserved = fread(fid, 1, 'uint16', 0, MACHINE_FORMAT);
                    r.variable_length_records(j).value(k).data_type = fread(fid, 1, 'uint8', 0, MACHINE_FORMAT);
                    extra_byte_options = fread(fid, 1, 'uint8', 0, MACHINE_FORMAT);
                    r.variable_length_records(j).value(k).options.no_data_bit = bitget(extra_byte_options, 1, 'uint8');
                    r.variable_length_records(j).value(k).options.min_bit = bitget(extra_byte_options, 2, 'uint8');
                    r.variable_length_records(j).value(k).options.max_bit = bitget(extra_byte_options, 3, 'uint8');
                    r.variable_length_records(j).value(k).options.scale_bit = bitget(extra_byte_options, 4, 'uint8');
                    r.variable_length_records(j).value(k).options.offset_bit = bitget(extra_byte_options, 5, 'uint8');
                    r.variable_length_records(j).value(k).name = regexp(fread(fid, 32, '*char', 0, MACHINE_FORMAT)', '[^\0]*', 'match', 'once');
                    r.variable_length_records(j).value(k).unused = fread(fid, 1, 'uint32', 0, MACHINE_FORMAT);
                    r.variable_length_records(j).value(k).no_data = fread(fid, 3, 'double', 0, MACHINE_FORMAT);
                    r.variable_length_records(j).value(k).min = fread(fid, 3, 'double', 0, MACHINE_FORMAT);
                    r.variable_length_records(j).value(k).max = fread(fid, 3, 'double', 0, MACHINE_FORMAT);
                    r.variable_length_records(j).value(k).scale = fread(fid, 3, 'double', 0, MACHINE_FORMAT);
                    r.variable_length_records(j).value(k).offset = fread(fid, 3, 'double', 0, MACHINE_FORMAT);
                    r.variable_length_records(j).value(k).description = regexp(fread(fid, 32, '*char', 0, MACHINE_FORMAT)', '[^\0]*', 'match', 'once');
                    
                    % append extra bytes to the point data record
                    ind_data_type = r.variable_length_records(j).value(k).data_type + 1;
                    
                    ind_record = ind_record + 1;
                    ind_last_record = ind_last_record + 1;
                    
                    r.record(ind_record).full_name = r.variable_length_records(j).value(k).name;
                    r.record(ind_record).short_name = lower(strrep(strtrim(regexprep(r.variable_length_records(j).value(k).name,'[^a-zA-Z0-9_\s]','')), ' ', '_'));
                    r.record(ind_record).format_compatibility = r.record(ind_record - 1).format_compatibility;
                    r.record(ind_record).bit_position = r.record(ind_last_record - 1).bit_position + (r.record(ind_last_record - 1).byte_length * 8);
                    r.record(ind_record).byte_position =  r.record(ind_record).bit_position ./ 8;
                    r.record(ind_record).type = {r.extra_bytes(ind_data_type).type};
                    r.record(ind_record).storage_type = {r.extra_bytes(ind_data_type).type}; %{'double'};
                    r.record(ind_record).byte_length = r.extra_bytes(ind_data_type).n_values * r.extra_bytes(ind_data_type).byte_length;
                    r.record(ind_record).flag_bit_field = false;
                    r.record(ind_record).flag_mandatory = false;
                    r.record(ind_record).flag_transform = true;
                    r.record(ind_record).value = [];
                    
                    if r.variable_length_records(j).value(k).options.scale_bit
                        
                        r.record(ind_record).scale = r.variable_length_records(j).value(k).scale(1);
                        
                    else
                        
                        r.record(ind_record).scale = 1;
                        
                    end
                    
                    if r.variable_length_records(j).value(k).options.offset_bit
                        
                        r.record(ind_record).offset = r.variable_length_records(j).value(k).offset(1);
                        
                    else
                        
                        r.record(ind_record).offset = 0;
                        
                    end
                    
                    if arg.Results.verbose
                        
                        fprintf('WARNING: The point data record contains an extra %i byte field: "%s"\n', r.record(ind_record).byte_length, r.record(ind_record).short_name)
                        
                    end
                    
                end
                
            case 7 % Superseded (optional)
                
                
            case num2cell(100:354) % Waveform Packet descriptor (required when using point formats 4, 5, 9, 10)
                
                r.variable_length_records(j).value.bits_per_sample = fread(fid, 1, 'uint8', 0, MACHINE_FORMAT); % Bits per sample, Unsigned char, 1 byte, *
                r.variable_length_records(j).value.compression_type = fread(fid, 1, 'uint8', 0, MACHINE_FORMAT); % Waveform compression type, Unsigned char, 1 byte, *
                r.variable_length_records(j).value.number_samples = fread(fid, 1, 'uint32', 0, MACHINE_FORMAT); % Number of samples, Unsigned long, 4 bytes, *
                r.variable_length_records(j).value.temporal_sample_spacing = fread(fid, 1, 'uint32', 0, MACHINE_FORMAT); % Temporal Sample Spacing, Unsigned long, 4 bytes, *
                r.variable_length_records(j).value.digitizer_gain = fread(fid, 1, 'double', 0, MACHINE_FORMAT); % Digitizer Gain, double, 8 bytes, *
                r.variable_length_records(j).value.digitizer_offset = fread(fid, 1, 'double', 0, MACHINE_FORMAT); % Digitizer Offset, double, 8 bytes, *
                
            case 2111 % OGC Math Transform WKT Record (optional)
                
                r.variable_length_records(j).value = regexp(fread(fid, r.variable_length_records(j).record_length_after_header, '*char', 0, MACHINE_FORMAT)', '[^\0]*', 'match', 'once');
                
            case 2112 % OGC Coordinate System WKT Record (optional)
                
                r.variable_length_records(j).value = regexp(fread(fid, r.variable_length_records(j).record_length_after_header, '*char', 0, MACHINE_FORMAT)', '[^\0]*', 'match', 'once');
                
            case 34735 % GeoKeyDirectoryTag Record (mandatory)
                
                r.variable_length_records(j).value.w_key_directory_version = fread(fid, 1, 'uint16', 0, MACHINE_FORMAT); % wKeyDirectoryVersion, unsigned short, wKeyDirectoryVersion = 1 Always
                r.variable_length_records(j).value.w_key_revision = fread(fid, 1, 'uint16', 0, MACHINE_FORMAT); % wKeyRevision, unsigned short, wKeyRevision = 1 Always
                r.variable_length_records(j).value.w_minor_revision = fread(fid, 1, 'uint16', 0, MACHINE_FORMAT); % wMinorRevision, unsigned short, wMinorRevision = 0 Always
                r.variable_length_records(j).value.w_number_of_keys = fread(fid, 1, 'uint16', 0, MACHINE_FORMAT); % wNumberOfKeys, unsigned short, Number of sets of 4 unsigned shorts to follow
                
                for k = 1:r.variable_length_records(j).value.w_number_of_keys
                    
                    r.variable_length_records(j).value.s_key_entry{k} = fread(fid, 4, 'uint16', 0, MACHINE_FORMAT);
                    
                end
                
            case 34736 % GeoDoubleParamsTag Record (optional)
                
                r.variable_length_records(j).value = fread(fid, r.variable_length_records(j).record_length_after_header / 8, 'double', 0, MACHINE_FORMAT);
                
            case 34737 % GeoAsciiParamsTag Record (optional)
                
                r.variable_length_records(j).value = sprintf('%c', fread(fid, r.variable_length_records(j).record_length_after_header, 'char', 0, MACHINE_FORMAT));
                
                % You may add custom records here
                
            otherwise % Other
                
                r.variable_length_records(j).value = fread(fid, r.variable_length_records(j).record_length_after_header, 'uint8', 0, MACHINE_FORMAT);
                
        end
        
        byte_offset = byte_offset + 54 + r.variable_length_records(j).record_length_after_header;
        
    end
    
    % read remaining bytes in header
    user_defined_bytes = r.header(phb_skeys.header_size).value - byte_offset;
    
    if user_defined_bytes > 0
        
        r.header.user_defined_bytes = fread(fid, user_defined_bytes, '*uint8', 0, MACHINE_FORMAT); % user defined bytes
        
        if arg.Results.verbose
            
            %warning(fprintf('The public header block is followed by %u user-defined bytes\n', user_defined_bytes));
            fprintf('WARNING: The public header block is followed by %u user-defined bytes\n', user_defined_bytes);
            
        end
        
    end
    
else
    
    % check number of variable length records
    if arg.Results.verbose
        
        %warning('The mandatory GeoKeyDirectoryTag was not specified in the LAS variable length records');
        fprintf('WARNING: The mandatory GeoKeyDirectoryTag was not specified in the LAS variable length records\n');
        
    end
    
end


%% read point data record

if ~arg.Results.headerOnly && r.header(phb_skeys.n_point_records).value > 0
    
    if arg.Results.verbose
        
        tic; % start timer
        fprintf('reading point data record...\n');
        
    end
    
    % generate key-value pairs for point data record
    for j = 1:length(r.record)
        
        pdr_skeys.(r.record(j).short_name) = j;
        
    end
    
    % read blob
    fseek(fid, r.header(phb_skeys.offset_to_data).value, 'bof');
    blob = fread(fid, '*uint8');
    
    % decode blob
    step = r.header(phb_skeys.point_data_record_length).value; % step is data record length
    ind_end = double(r.header(phb_skeys.point_data_record_length).value) * double(r.header(phb_skeys.n_point_records).value); % end
    
    % read byte chunks from blob
    for j = 1:size(r.record,2)
        
        ind_start = floor(r.record(j).byte_position) + 1;
        
        nbitshift = 0;
        
        if r.record(j).flag_bit_field
            
            bytes = blob(ind_start : step : ind_end); % cast to uint8
            shift_left = floor(r.record(j).byte_position + 1) * 8 - (r.record(j).bit_position + r.record(j).byte_length * 8);
            shift_right =  r.record(j).byte_length * 8 - 8;
            r.record(j).value = cast(bitshift(bitshift(bytes, shift_left), shift_right), r.record(j).storage_type{:});
            
        else
            
            % initialize byte array
            bytes = zeros(r.header(phb_skeys.n_point_records).value, r.record(j).byte_length, r.record(j).type{:});
            
            % add columns to byte array (one column per byte)
            for k = 1:r.record(j).byte_length

                bytes(:,k) = bitshift(cast(blob(ind_start: step : ind_end), r.record(j).type{:}), nbitshift);
                ind_start = ind_start + 1; % increment byte address
                nbitshift = nbitshift + 8;  % increment bit shift
                
            end
            
            % merge byte array columns
            byte_a = bytes(:,1);
            if r.record(j).byte_length > 1
                
                for k = 1:r.record(j).byte_length-1
                    
                    byte_a = bitor(byte_a, bytes(:,k+1));
                    
                end
                
            end
     
            % adjust format
            switch r.record(j).storage_type{:}
                
                case {'single','double'}
                    
                    r.record(j).value = typecast(byte_a, r.record(j).storage_type{:});
                    
                case {'uchar'}
                    
                    r.record(j).value = cast(byte_a, 'int8');
                    idx_sign = (byte_a >= 128);
                    r.record(j).value(idx_sign) = int8(int16(byte_a(idx_sign)) - 256);
                    
                otherwise
                    
                    r.record(j).value = cast(byte_a, r.record(j).storage_type{:});
                    
            end
            
        end
        
    end
    
    clear byte byte_a
    
    % apply scale and offset when necessary
    r.record(pdr_skeys.x).scale = r.header(phb_skeys.x_scale_factor).value;
    r.record(pdr_skeys.x).offset = r.header(phb_skeys.x_offset).value;
    r.record(pdr_skeys.y).scale = r.header(phb_skeys.y_scale_factor).value;
    r.record(pdr_skeys.y).offset = r.header(phb_skeys.y_offset).value;
    r.record(pdr_skeys.z).scale = r.header(phb_skeys.z_scale_factor).value;
    r.record(pdr_skeys.z).offset = r.header(phb_skeys.z_offset).value;
    
    for j = 1:length(r.record)
        
        if r.record(j).flag_transform

            r.record(j).value = double(r.record(j).value) * r.record(j).scale + r.record(j).offset;
            
        end
        
    end
    
    % adjust scan angle format
    switch r.record(pdr_skeys.scan_angle).type{:}
        
        case 'uint8' %'uint8'
            
            % Scan Angle Rank is the scan angle rounded to the nearest integer from from +90° to -90°,
            % with -ve being left of 0° (nadir) looking forward
            idx_scan_angle_sign = r.record(pdr_skeys.scan_angle).value >= 128;
            r.record(pdr_skeys.scan_angle).value(idx_scan_angle_sign) = r.record(pdr_skeys.scan_angle).value(idx_scan_angle_sign) - 256;
            
        case 'int16'
            
            % Scan Angle is a rotational position in 0.006° increments with values from +30,000 to -30,000
            % covering +180° to -180° (0° is at nadir/down, -ve is left/clockwise of nadir looking forward)
            r.record(pdr_skeys.scan_angle).value = double(r.record(pdr_skeys.scan_angle).value) / 15000 * 90;
            
    end
    
end


%% read extended variable length record (waveform data packets and custom)

flag_evlr = false;

if ~arg.Results.headerOnly && (las_version >= 14)
    
    if r.header(phb_skeys.number_of_evlr).value > 0
        
        byte_offset = r.header(phb_skeys.offset_to_evlr).value;
        
        for j = 1:r.header(phb_skeys.number_of_evlr).value
            
            fseek(fid, byte_offset, 'bof');
            
            % read common header part
            r.extended_variable_length_records(j).reserved = fread(fid, 1, 'uint16', 0, MACHINE_FORMAT); % Reserved, unsigned short, 2 bytes;
            r.extended_variable_length_records(j).user_id = regexp(fread(fid, 16, '*char', 0, MACHINE_FORMAT)', '[^\0]*', 'match', 'once'); % read until null terminator
            r.extended_variable_length_records(j).record_id = fread(fid, 1, 'uint16', 0, MACHINE_FORMAT); % Record ID, unsigned short, 2 bytes, *
            r.extended_variable_length_records(j).record_length_after_header = fread(fid, 1, 'uint64', 0, MACHINE_FORMAT); % Record Length After Header, unsigned long long, 8 bytes, *
            r.extended_variable_length_records(j).description = regexp(fread(fid, 32, '*char', 0, MACHINE_FORMAT)', '[^\0]*', 'match', 'once'); % Description, char[32], 32 bytes
            
            % read any defined extended variable length records
            switch r.extended_variable_length_records(j).record_id
                
                case 0 % Classification lookup (optional)
                    
                    r.extended_variable_length_records(j).value.class_number = fread(fid, r.extended_variable_length_records(j).record_length_after_header, '*uint8', 0, MACHINE_FORMAT);
                    r.extended_variable_length_records(j).value.class_description = sprintf('%c', fread(fid, r.extended_variable_length_records(j).record_length_after_header, 'char', 0, MACHINE_FORMAT));
                    
                case 1 % Reserved (optional)
                    
                    
                case 2 % Histogram (optional)
                    
                    
                case 3 % Text area description (optional)
                    
                    r.extended_variable_length_records(j).value = regexp(fread(fid, r.extended_variable_length_records(j).record_length_after_header, '*char', 0, MACHINE_FORMAT)', '[^\0]*', 'match', 'once');
                    
                case 4 % Extra bytes (optional)
                    
                    n_extra_records = r.extended_variable_length_records(j).record_length_after_header / 192;
                    
                    for k = 1:n_extra_records
                        
                        r.extended_variable_length_records(j).value(k).reserved = fread(fid, 1, 'uint16', 0, MACHINE_FORMAT);
                        r.extended_variable_length_records(j).value(k).data_type = fread(fid, 1, 'uint8', 0, MACHINE_FORMAT);
                        
                        extra_byte_options = fread(fid, 1, 'uint8', 0, MACHINE_FORMAT);
                        r.extended_variable_length_records(j).value(k).options.no_data_bit = bitget(extra_byte_options, 1, 'uint8');
                        r.extended_variable_length_records(j).value(k).options.min_bit = bitget(extra_byte_options, 2, 'uint8');
                        r.extended_variable_length_records(j).value(k).options.max_bit = bitget(extra_byte_options, 3, 'uint8');
                        r.extended_variable_length_records(j).value(k).options.scale_bit = bitget(extra_byte_options, 4, 'uint8');
                        r.extended_variable_length_records(j).value(k).options.offset_bit = bitget(extra_byte_options, 5, 'uint8');
                        
                        r.extended_variable_length_records(j).value(k).name = regexp(fread(fid, 32, '*char', 0, MACHINE_FORMAT)', '[^\0]*', 'match', 'once');
                        r.extended_variable_length_records(j).value(k).unused = fread(fid, 1, 'uint32', 0, MACHINE_FORMAT);
                        r.extended_variable_length_records(j).value(k).no_data = fread(fid, 3, 'double', 0, MACHINE_FORMAT);
                        r.extended_variable_length_records(j).value(k).min = fread(fid, 3, 'double', 0, MACHINE_FORMAT);
                        r.extended_variable_length_records(j).value(k).max = fread(fid, 3, 'double', 0, MACHINE_FORMAT);
                        r.extended_variable_length_records(j).value(k).scale = fread(fid, 3, 'double', 0, MACHINE_FORMAT);
                        r.extended_variable_length_records(j).value(k).offset = fread(fid, 3, 'double', 0, MACHINE_FORMAT);
                        r.extended_variable_length_records(j).value(k).description = regexp(fread(fid, 32, '*char', 0, MACHINE_FORMAT)', '[^\0]*', 'match', 'once');
                        
                    end
                    
                case 7 % Superseded (optional)
                    
                    
                case num2cell(100:354) % Waveform Packet descriptor (required when using point formats 4, 5, 9, 10)
                    
                    r.extended_variable_length_records(j).value.bits_per_sample = fread(fid, 1, 'uint8', 0, MACHINE_FORMAT); % Bits per sample, Unsigned char, 1 byte, *
                    r.extended_variable_length_records(j).value.compression_type = fread(fid, 1, 'uint8', 0, MACHINE_FORMAT); % Waveform compression type, Unsigned char, 1 byte, *
                    r.extended_variable_length_records(j).value.number_samples = fread(fid, 1, 'uint32', 0, MACHINE_FORMAT); % Number of samples, Unsigned long, 4 bytes, *
                    r.extended_variable_length_records(j).value.temporal_sample_spacing = fread(fid, 1, 'uint32', 0, MACHINE_FORMAT); % Temporal Sample Spacing, Unsigned long, 4 bytes, *
                    r.extended_variable_length_records(j).value.digitizer_gain = fread(fid, 1, 'double', 0, MACHINE_FORMAT); % Digitizer Gain, double, 8 bytes, *
                    r.extended_variable_length_records(j).value.digitizer_offset = fread(fid, 1, 'double', 0, MACHINE_FORMAT); % Digitizer Offset, double, 8 bytes, *
                    
                case 2111 % OGC Math Transform WKT Record (optional)
                    
                    r.extended_variable_length_records(j).value = regexp(fread(fid, r.extended_variable_length_records(j).record_length_after_header, '*char', 0, MACHINE_FORMAT)', '[^\0]*', 'match', 'once');
                    
                case 2112 % OGC Coordinate System WKT Record (optional)
                    
                    r.extended_variable_length_records(j).record_length_after_header
                    r.extended_variable_length_records(j).value = regexp(fread(fid, r.extended_variable_length_records(j).record_length_after_header, '*char', 0, MACHINE_FORMAT)', '[^\0]*', 'match', 'once');
                    
                case 34735 % GeoKeyDirectoryTag Record (mandatory)
                    
                    r.extended_variable_length_records(j).value.w_key_directory_version = fread(fid, 1, 'uint16', 0, MACHINE_FORMAT); % wKeyDirectoryVersion, unsigned short, wKeyDirectoryVersion = 1 Always
                    r.extended_variable_length_records(j).value.w_key_revision = fread(fid, 1, 'uint16', 0, MACHINE_FORMAT); % wKeyRevision, unsigned short, wKeyRevision = 1 Always
                    r.extended_variable_length_records(j).value.w_minor_revision = fread(fid, 1, 'uint16', 0, MACHINE_FORMAT); % wMinorRevision, unsigned short, wMinorRevision = 0 Always
                    r.extended_variable_length_records(j).value.w_number_of_keys = fread(fid, 1, 'uint16', 0, MACHINE_FORMAT); % wNumberOfKeys, unsigned short, Number of sets of 4 unsigned shorts to follow
                    
                    for k = 1:r.extended_variable_length_records(j).value.w_number_of_keys
                        
                        r.extended_variable_length_records(j).value.s_key_entry{k} = fread(fid, 4, 'uint16', 0, MACHINE_FORMAT);
                        
                    end
                    
                case 34736 % GeoDoubleParamsTag Record (optional)
                    
                    r.extended_variable_length_records(j).value = fread(fid, r.extended_variable_length_records(j).record_length_after_header / 8, 'double', 0, MACHINE_FORMAT);
                    
                case 34737 % GeoAsciiParamsTag Record (optional)
                    
                    r.extended_variable_length_records(j).value = sprintf('%c', fread(fid, r.extended_variable_length_records(j).record_length_after_header, 'char', 0, MACHINE_FORMAT));
                    
                case 65535 % Waveform data packets (internal/external storage)
                    
                    if r.header(phb_skeys.global_encoding_waveform_data_packets_internal).value
                        
                        % internal waveform storage
                        %r.extended_variable_length_records(j).value = fread(fid, r.extended_variable_length_records(j).record_length_after_header / 8, 'double', 0, MACHINE_FORMAT);
                        
                    else
                        
                        % external waveform storage
                        
                        
                    end
                    
                case {5005, 5006} % Custom uint8 field (optional)
                    
                    r.extended_variable_length_records(j).value = fread(fid, r.extended_variable_length_records(j).record_length_after_header, 'uint8', 0, MACHINE_FORMAT);
                    
                case {5001, 5002, 5003, 5004, 5007} % Custom uint16 field (optional)
                    
                    r.extended_variable_length_records(j).value = fread(fid, r.extended_variable_length_records(j).record_length_after_header / 2, 'uint16', 0, MACHINE_FORMAT);
                    
                case {5000, 5008, 5009} % Custom uint32 field (optional)
                    
                    r.extended_variable_length_records(j).value = fread(fid, r.extended_variable_length_records(j).record_length_after_header / 4, 'uint32', 0, MACHINE_FORMAT);
                    
                    % You may add custom records here
                    
                otherwise % Other
                    
                    r.extended_variable_length_records(j).value = fread(fid, r.extended_variable_length_records(j).record_length_after_header, 'uint8', 0, MACHINE_FORMAT);
                    
            end
            
            byte_offset = byte_offset + 60 + r.extended_variable_length_records(j).record_length_after_header;
            
        end
        
        flag_evlr = true;
        
    end
    
end


%% format data for export

% public header block
for j = 1:length(r.header)
    
    s.header.(r.header(j).short_name) = r.header(j).value;
    
end

% variable length records
if flag_vlr
    
    s.variable_length_records = r.variable_length_records;
    
end

% point data record
if ~arg.Results.headerOnly
    
    for j = 1:length(r.record)
        
        s.record.(r.record(j).short_name) = r.record(j).value;
        
    end
    
end

% extended variable length records
if flag_evlr
    
    s.extended_variable_length_records = r.extended_variable_length_records;
    
end


%% close file

if arg.Results.verbose && ~arg.Results.headerOnly
    
    tElapsed = toc; % stop timer
    fprintf('%u point records read in %s\n', r.header(phb_skeys.n_point_records).value, datestr(tElapsed/(24*3600), 'HH:MM:SS.FFF'));
    
end

fclose(fid);

end
