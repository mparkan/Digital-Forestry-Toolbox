function trajectory = TRJread(filepath, varargin)
%TRJREAD - Reads a trajectory file in the Terrascan TRJ format.
% TRAJECTORY = TRJREAD(FILEPATH, HEADERONLY, VERBOSE) reads the TRJ file specified in FILEPATH
% into the structure TRAJECTORY.
%
% The Terrascan TRJ format stores trajectories as binary files with TRJ extension. These files
% contain a header followed by a number of trajectory position records.
%
% For more information see: Terrasolid, TerraScan User's Guide,
% https://www.terrasolid.com/download/tscan.pdf
%
% Syntax:  trajectory = TRJread(filepath, headerOnly, verbose)
%
% Inputs:
%    filepath - The path to the input TRJ file
%    headerOnly - If set to true, will read only the header information
%    verbose - If set to true, will display the header information
%              in the console
%
% Outputs:
%    trajectory - A structure containing the header information and trajectory records
%
% Example:
%    trajectory = TRJread('E:\trajectories\1701.trj', ...
%        'headerOnly', false, ...
%        'verbose', true);
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

machineformat = 'ieee-le'; % little-endian format

arg = inputParser;

addRequired(arg, 'filepath', @ischar);
addParameter(arg, 'headerOnly', false, @(x) islogical(x) && (numel(x) == 1));
addParameter(arg, 'verbose', true, @(x) islogical(x) && (numel(x) == 1));

parse(arg, filepath, varargin{:});

% check file format
fclose('all');
fid = fopen(filepath, 'r');

if fid == -1
    
    error('Could not open file');
    
end

if ~strcmp(sprintf('%c', fread(fid, 8, 'char', 0, machineformat)), 'TSCANTRJ')
    
    error('Input file is not a valid TRJ file');
    
end

fclose(fid);


%% header block format definition

% File Signature (TSCANTRJ), 8 bytes
k = 1;
r.TrajHdr(k).full_name = 'File Signature (TSCANTRJ)';
r.TrajHdr(k).short_name = 'Recog';
r.TrajHdr(k).type = 'char';
r.TrajHdr(k).byte_length = 1;
r.TrajHdr(k).n_values = 8;
r.TrajHdr(k).print_format = '%s';
r.TrajHdr(k).value = [];
k = k + 1;

% File version, 4 bytes
r.TrajHdr(k).full_name = 'Version';
r.TrajHdr(k).short_name = 'Version';
r.TrajHdr(k).type = 'int32';
r.TrajHdr(k).byte_length = 4;
r.TrajHdr(k).n_values = 1;
r.TrajHdr(k).print_format = '%u';
r.TrajHdr(k).value = [];
k = k + 1;

% TrajHdr size, sizeof(TrajHdr), 4 bytes
r.TrajHdr(k).full_name = 'Header size';
r.TrajHdr(k).short_name = 'HdrSize';
r.TrajHdr(k).type = 'uint32';
r.TrajHdr(k).byte_length = 4;
r.TrajHdr(k).n_values = 1;
r.TrajHdr(k).print_format = '%u';
r.TrajHdr(k).value = [];
k = k + 1;

% Number of position records, 4 bytes
r.TrajHdr(k).full_name = 'Number of position records';
r.TrajHdr(k).short_name = 'PosCnt';
r.TrajHdr(k).type = 'int32';
r.TrajHdr(k).byte_length = 4;
r.TrajHdr(k).n_values = 1;
r.TrajHdr(k).print_format = '%u';
r.TrajHdr(k).value = [];
k = k + 1;

% Size of position records, 4 bytes
r.TrajHdr(k).full_name = 'Size of position records';
r.TrajHdr(k).short_name = 'PosSize';
r.TrajHdr(k).type = 'int32';
r.TrajHdr(k).byte_length = 4;
r.TrajHdr(k).n_values = 1;
r.TrajHdr(k).print_format = '%u';
r.TrajHdr(k).value = [];
k = k + 1;

% Description, 78 bytes
r.TrajHdr(k).full_name = 'Description';
r.TrajHdr(k).short_name = 'Desc';
r.TrajHdr(k).type = 'char';
r.TrajHdr(k).byte_length = 1;
r.TrajHdr(k).n_values = 78;
r.TrajHdr(k).print_format = '%s';
r.TrajHdr(k).value = [];
k = k + 1;

% System identifier (for lever arms), 1 byte
r.TrajHdr(k).full_name = 'System identifier (for lever arms)';
r.TrajHdr(k).short_name = 'SysIdv';
r.TrajHdr(k).type = 'int8';
r.TrajHdr(k).byte_length = 1;
r.TrajHdr(k).n_values = 1;
r.TrajHdr(k).print_format = '%u';
r.TrajHdr(k).value = [];
k = k + 1;

% Quality for whole trajectory (1-5), 1 byte
r.TrajHdr(k).full_name = 'Quality for whole trajectory (1-5)';
r.TrajHdr(k).short_name = 'Quality';
r.TrajHdr(k).type = 'int8';
r.TrajHdr(k).byte_length = 1;
r.TrajHdr(k).n_values = 1;
r.TrajHdr(k).print_format = '%u';
r.TrajHdr(k).value = [];
k = k + 1;

% First time stamp, 8 bytes
r.TrajHdr(k).full_name = 'First time stamp';
r.TrajHdr(k).short_name = 'BegTime';
r.TrajHdr(k).type = 'double';
r.TrajHdr(k).byte_length = 8;
r.TrajHdr(k).n_values = 1;
r.TrajHdr(k).print_format = '%g';
r.TrajHdr(k).value = [];
k = k + 1;

% Last time stamp, 8 bytes
r.TrajHdr(k).full_name = 'Last time stamp';
r.TrajHdr(k).short_name = 'EndTime';
r.TrajHdr(k).type = 'double';
r.TrajHdr(k).byte_length = 8;
r.TrajHdr(k).n_values = 1;
r.TrajHdr(k).print_format = '%g';
r.TrajHdr(k).value = [];
k = k + 1;

% Original number (before any splitting), 4 bytes
r.TrajHdr(k).full_name = 'Original number (before any splitting)';
r.TrajHdr(k).short_name = 'OrigNbr';
r.TrajHdr(k).type = 'int32';
r.TrajHdr(k).byte_length = 4;
r.TrajHdr(k).n_values = 1;
r.TrajHdr(k).print_format = '%u';
r.TrajHdr(k).value = [];
k = k + 1;

% Flightline number (in laser points), 4 bytes
r.TrajHdr(k).full_name = 'Flightline number (in laser points)';
r.TrajHdr(k).short_name = 'Number';
r.TrajHdr(k).type = 'int32';
r.TrajHdr(k).byte_length = 4;
r.TrajHdr(k).n_values = 1;
r.TrajHdr(k).print_format = '%u';
r.TrajHdr(k).value = [];
k = k + 1;

% Vertical facing video, 400 bytes
r.TrajHdr(k).full_name = 'Vertical facing video';
r.TrajHdr(k).short_name = 'VrtVideo';
r.TrajHdr(k).type = 'char';
r.TrajHdr(k).byte_length = 1;
r.TrajHdr(k).n_values = 400;
r.TrajHdr(k).print_format = '%s';
r.TrajHdr(k).value = [];
k = k + 1;

% Start time of VrtVideo[], 8 bytes
r.TrajHdr(k).full_name = 'Start time of VrtVideo[]';
r.TrajHdr(k).short_name = 'VrtBeg';
r.TrajHdr(k).type = 'double';
r.TrajHdr(k).byte_length = 8;
r.TrajHdr(k).n_values = 1;
r.TrajHdr(k).print_format = '%g';
r.TrajHdr(k).value = [];
k = k + 1;

% End time of VrtVideo[], 8 bytes
r.TrajHdr(k).full_name = 'End time of VrtVideo[]';
r.TrajHdr(k).short_name = 'VrtEnd';
r.TrajHdr(k).type = 'double';
r.TrajHdr(k).byte_length = 8;
r.TrajHdr(k).n_values = 1;
r.TrajHdr(k).print_format = '%g';
r.TrajHdr(k).value = [];
k = k + 1;

% Forward facing video, 400 bytes
r.TrajHdr(k).full_name = 'Forward facing video';
r.TrajHdr(k).short_name = 'FwdVideo';
r.TrajHdr(k).type = 'char';
r.TrajHdr(k).byte_length = 1;
r.TrajHdr(k).n_values = 400;
r.TrajHdr(k).print_format = '%s';
r.TrajHdr(k).value = [];
k = k + 1;

% Start time of FwdVideo[], 8 bytes
r.TrajHdr(k).full_name = 'Start time of FwdVideo[]';
r.TrajHdr(k).short_name = 'FwdBeg';
r.TrajHdr(k).type = 'double';
r.TrajHdr(k).byte_length = 8;
r.TrajHdr(k).n_values = 1;
r.TrajHdr(k).print_format = '%g';
r.TrajHdr(k).value = [];
k = k + 1;

% End time of FwdVideo[], 8 bytes
r.TrajHdr(k).full_name = 'End time of FwdVideo[]';
r.TrajHdr(k).short_name = 'FwdEnd';
r.TrajHdr(k).type = 'double';
r.TrajHdr(k).byte_length = 8;
r.TrajHdr(k).n_values = 1;
r.TrajHdr(k).print_format = '%g';
r.TrajHdr(k).value = [];
k = k + 1;

% Waveform data file, 400 bytes
r.TrajHdr(k).full_name = 'Waveform data file';
r.TrajHdr(k).short_name = 'WaveFile';
r.TrajHdr(k).type = 'char';
r.TrajHdr(k).byte_length = 1;
r.TrajHdr(k).n_values = 400;
r.TrajHdr(k).print_format = '%s';
r.TrajHdr(k).value = [];
k = k + 1;

% Group (session description), 16 bytes
r.TrajHdr(k).full_name = 'Group (session description)';
r.TrajHdr(k).short_name = 'Group';
r.TrajHdr(k).type = 'char';
r.TrajHdr(k).byte_length = 1;
r.TrajHdr(k).n_values = 16;
r.TrajHdr(k).print_format = '%s';
r.TrajHdr(k).value = [];

%% generate key-value pairs for header block

for j = 1:length(r.TrajHdr)
    
    TrajHdr_skeys.(r.TrajHdr(j).short_name) = j;
    
end

%% open file

fid = fopen(filepath);

if fid == -1
    error('Could not open file')
end

%% read header block

fseek(fid, 12, 'bof');
r.TrajHdr(TrajHdr_skeys.HdrSize).value = fread(fid, 1, 'uint32', 0, machineformat); % read header size

byte_offset = 0;
fseek(fid, byte_offset, 'bof');
j = 1;

while byte_offset < r.TrajHdr(TrajHdr_skeys.HdrSize).value
    
    if ismember(r.TrajHdr(j).type, {'uchar','char'})
        
        val = fread(fid, r.TrajHdr(j).n_values, r.TrajHdr(j).type, 0, machineformat);
        str = strtrim(cellstr(sprintf('%c', val)));
        r.TrajHdr(j).value = str{:};

    else
        
        r.TrajHdr(j).value = fread(fid, r.TrajHdr(j).n_values, r.TrajHdr(j).type, 0, machineformat);
        
    end
    
    if arg.Results.verbose
        
        fprintf(sprintf('%s: %s\\n', r.TrajHdr(j).full_name, r.TrajHdr(j).print_format), r.TrajHdr(j).value);
        
    end
    
    byte_offset = byte_offset + r.TrajHdr(j).n_values * r.TrajHdr(j).byte_length;
    j = j + 1;
    
end

%% trajectory position record block format definition

% Time stamp (seconds in some system), 8 bytes
k = 1;
r.TrajPos(k).full_name = 'Time stamp (seconds in some system)';
r.TrajPos(k).short_name = 'Time';
r.TrajPos(k).type = 'int64';
r.TrajPos(k).byte_length = 8;
r.TrajPos(k).mandatory = true;
r.TrajPos(k).storage_type = 'double';
r.TrajPos(k).value = [];
k = k + 1;

% X, 8 bytes
r.TrajPos(k).full_name = 'X';
r.TrajPos(k).short_name = 'X';
r.TrajPos(k).type = 'int64';
r.TrajPos(k).byte_length = 8;
r.TrajPos(k).mandatory = true;
r.TrajPos(k).storage_type = 'double';
r.TrajPos(k).value = [];
k = k + 1;

% Y, 8 bytes
r.TrajPos(k).full_name = 'Y';
r.TrajPos(k).short_name = 'Y';
r.TrajPos(k).type = 'int64';
r.TrajPos(k).byte_length = 8;
r.TrajPos(k).mandatory = true;
r.TrajPos(k).storage_type = 'double';
r.TrajPos(k).value = [];
k = k + 1;

% Z, 8 bytes
r.TrajPos(k).full_name = 'Z';
r.TrajPos(k).short_name = 'Z';
r.TrajPos(k).type = 'int64';
r.TrajPos(k).byte_length = 8;
r.TrajPos(k).mandatory = true;
r.TrajPos(k).storage_type = 'double';
r.TrajPos(k).value = [];
k = k + 1;

% Heading (degrees), 8 bytes
r.TrajPos(k).full_name = 'Heading (degrees)';
r.TrajPos(k).short_name = 'Head';
r.TrajPos(k).type = 'int64';
r.TrajPos(k).byte_length = 8;
r.TrajPos(k).mandatory = true;
r.TrajPos(k).storage_type = 'double';
r.TrajPos(k).value = [];
k = k + 1;

% Roll (degrees), 8 bytes
r.TrajPos(k).full_name = 'Roll (degrees)';
r.TrajPos(k).short_name = 'Roll';
r.TrajPos(k).type = 'int64';
r.TrajPos(k).byte_length = 8;
r.TrajPos(k).mandatory = true;
r.TrajPos(k).storage_type = 'double';
r.TrajPos(k).value = [];
k = k + 1;

% Pitch (degrees), 8 bytes
r.TrajPos(k).full_name = 'Pitch (degrees)';
r.TrajPos(k).short_name = 'Pitch';
r.TrajPos(k).type = 'int64';
r.TrajPos(k).byte_length = 8;
r.TrajPos(k).mandatory = true;
r.TrajPos(k).storage_type = 'double';
r.TrajPos(k).value = [];
k = k + 1;

% Quality for xy, 0=not set, 1 byte
r.TrajPos(k).full_name = 'Quality for xy, 0=not set';
r.TrajPos(k).short_name = 'QtyXy';
r.TrajPos(k).type = 'int8';
r.TrajPos(k).byte_length = 1;
r.TrajPos(k).mandatory = true;
r.TrajPos(k).storage_type = 'int8';
r.TrajPos(k).value = [];
k = k + 1;

% Quality for z, 0=not set, 1 byte
r.TrajPos(k).full_name = 'Quality for z, 0=not set';
r.TrajPos(k).short_name = 'QtyZ';
r.TrajPos(k).type = 'int8';
r.TrajPos(k).byte_length = 1;
r.TrajPos(k).mandatory = true;
r.TrajPos(k).storage_type = 'int8';
r.TrajPos(k).value = [];
k = k + 1;

% Quality for heading, 0=not set, 1 byte
r.TrajPos(k).full_name = 'Quality for headingy, 0=not set';
r.TrajPos(k).short_name = 'QtyH';
r.TrajPos(k).type = 'int8';
r.TrajPos(k).byte_length = 1;
r.TrajPos(k).mandatory = true;
r.TrajPos(k).storage_type = 'int8';
r.TrajPos(k).value = [];
k = k + 1;

% Quality for roll/pitch, 0=not set, 1 byte
r.TrajPos(k).full_name = 'Quality for roll/pitch, 0=not set';
r.TrajPos(k).short_name = 'QtyRp';
r.TrajPos(k).type = 'int8';
r.TrajPos(k).byte_length = 1;
r.TrajPos(k).mandatory = true;
r.TrajPos(k).storage_type = 'int8';
r.TrajPos(k).value = [];
k = k + 1;

% Run time flag, 2 bytes
r.TrajPos(k).full_name = 'Run time flag 1';
r.TrajPos(k).short_name = 'Mark';
r.TrajPos(k).type = 'int16';
r.TrajPos(k).byte_length = 2;
r.TrajPos(k).mandatory = true;
r.TrajPos(k).storage_type = 'int16';
r.TrajPos(k).value = [];
k = k + 1;

% Run time flag, 2 bytes
r.TrajPos(k).full_name = 'Run time flag 2';
r.TrajPos(k).short_name = 'Flag';
r.TrajPos(k).type = 'int16';
r.TrajPos(k).byte_length = 2;
r.TrajPos(k).mandatory = true;
r.TrajPos(k).storage_type = 'int16';
r.TrajPos(k).value = [];

%% read trajectory position record

if ~arg.Results.headerOnly
    
    if arg.Results.verbose
        
        fprintf('reading trajectory position record...');
        
    end
    
    % generate key-value pairs for point data record
    for j = 1:length(r.TrajPos)
        
        TrajPos_skeys.(r.TrajPos(j).short_name) = j;
        
    end
    
    % read blob
    byte_offset = r.TrajHdr(TrajHdr_skeys.HdrSize).value;
    fseek(fid, byte_offset, 'bof');
    blob = fread(fid, '*uint8');
    fclose(fid); % close file
    
    % parse blob
    idx_step = r.TrajHdr(TrajHdr_skeys.PosSize).value; % step
    idx_end = r.TrajHdr(TrajHdr_skeys.PosCnt).value * r.TrajHdr(TrajHdr_skeys.PosSize).value; % end
    offset = 1;
    
    for j = 1:size(r.TrajPos,2)
        
        byte_length = r.TrajPos(j).byte_length;
        nbitshift = 0;
        
        byte = zeros(r.TrajHdr(TrajHdr_skeys.PosCnt).value, byte_length, r.TrajPos(j).type);
        
        % fill byte array
        for k = 1:byte_length
            
            byte(:,k) = bitshift(cast(blob(offset : idx_step : idx_end), r.TrajPos(j).type), nbitshift);
            offset = offset + 1; % shift to next byte
            nbitshift = nbitshift + 8; % increase bitshift by 8 bits
            
        end
        
        byte_a = byte(:,1);
        
        if byte_length > 1
            
            for k = 1:byte_length-1
                
                byte_a = bitor(byte_a, byte(:,k+1));
                
            end
        end
        
        if ismember(r.TrajPos(j).storage_type, {'single','double'})
            
            r.TrajPos(j).value = typecast(byte_a, r.TrajPos(j).storage_type); % cast uint32/uint64 types to single/double precision
            
        else
            
            r.TrajPos(j).value = cast(byte_a, r.TrajPos(j).storage_type);
            
        end
        
    end
    
    clear byte byte_a
    
    % apply scale to quality values
    r.TrajPos(TrajPos_skeys.QtyXy).value = double(r.TrajPos(TrajPos_skeys.QtyXy).value).^1.5 * 0.001; % quality for xy in meters
    r.TrajPos(TrajPos_skeys.QtyZ).value = double(r.TrajPos(TrajPos_skeys.QtyZ).value).^1.5 * 0.001; % quality for z in meters
    r.TrajPos(TrajPos_skeys.QtyH).value = double(r.TrajPos(TrajPos_skeys.QtyH).value).^1.5 * 0.0001; % quality for heading in degrees
    r.TrajPos(TrajPos_skeys.QtyRp).value = double(r.TrajPos(TrajPos_skeys.QtyRp).value).^1.5 * 0.0001; % quality for roll/pitch in degrees
    
end

%% format data for export

% header block
for j = 1:length(r.TrajHdr)
    
    trajectory.header.(r.TrajHdr(j).short_name) = r.TrajHdr(j).value;
    
end

% record block
if ~arg.Results.headerOnly
    
    for j = 1:length(r.TrajPos)
        
        trajectory.record.(r.TrajPos(j).short_name) = r.TrajPos(j).value;
        
    end
    
end

if arg.Results.verbose
    
    fprintf('done!\n');
    
end

end
