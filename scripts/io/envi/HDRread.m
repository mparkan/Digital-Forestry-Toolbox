function [info] = HDRread(filepath, varargin)
%HDRREAD - Reads a ENVI header (.hdr) file
% [INFO] = HDRREAD(FILEPATH, ...) reads the ENVI header file (.hdr)
% specified in FILEPATH into the structure INFO.
%
% For more information, see: Exelis, ENVI Header Files,
% http://www.exelisvis.com/docs/ENVIHeaderFiles.html
%
% Syntax:  [info] = HDRread(filepath, ...)
%
% Inputs:
%    filepath - The path to the input ENVI header (.hdr) file
%    verbose (optional parameter) - If set to true, will display header attributes
%              in the console
%
% Outputs:
%    info - A structure containing the header attributes
%
% Example: 
%    info = HDRread('E:\images\image.hdr', 'verbose', true)
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

addRequired(arg, 'filepath', @(x) ischar(x));
addParameter(arg, 'verbose', true, @(x) islogical(x) && (numel(x) == 1));

parse(arg, filepath, varargin{:});


%% read header file

fclose('all');
fid = fopen(filepath, 'r');
if fid == -1
    
    error('Could not open header file')
    
end

str = fscanf(fid, '%c');
fclose(fid);


%% remove carriage returns (windows formatted eol -> Unix formatted eol)

str = regexprep(str,'\r\n','\n');


%% define ENVI header format

k = 1;
hdr(k).full_name = 'Acquisition time';
hdr(k).short_name = 'acquisition_time';
hdr(k).pattern = 'acquisition time';
hdr(k).type = 'str';
hdr(k).print_format = '%s';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Band names';
hdr(k).short_name = 'band_names';
hdr(k).pattern = 'band names';
hdr(k).type = 'str_arr';
hdr(k).print_format = '%s';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Bands';
hdr(k).short_name = 'bands';
hdr(k).pattern = 'bands';
hdr(k).type = 'num';
hdr(k).print_format = '%u';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Bad multiplier values';
hdr(k).short_name = 'bbl';
hdr(k).pattern = 'bbl';
hdr(k).type = 'num_arr';
hdr(k).print_format = '%u';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Byte order';
hdr(k).short_name = 'byte_order';
hdr(k).pattern = 'byte order';
hdr(k).type = 'num';
hdr(k).print_format = '%u';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Class lookup';
hdr(k).short_name = 'class_lookup';
hdr(k).pattern = 'class lookup';
hdr(k).type = 'str_arr';
hdr(k).print_format = '%s';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Class names';
hdr(k).short_name = 'class_names';
hdr(k).pattern = 'class names';
hdr(k).type = 'str_arr';
hdr(k).print_format = '%s';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Classes';
hdr(k).short_name = 'classes';
hdr(k).pattern = 'classes';
hdr(k).type = 'num';
hdr(k).print_format = '%f';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Cloud cover';
hdr(k).short_name = 'cloud_cover';
hdr(k).pattern = 'cloud cover';
hdr(k).type = 'num';
hdr(k).print_format = '%f';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Complex function';
hdr(k).short_name = 'complex_function';
hdr(k).pattern = 'complex function';
hdr(k).type = 'str_arr';
hdr(k).print_format = '%s';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Coordinate system string';
hdr(k).short_name = 'coordinate_system_string';
hdr(k).pattern = 'coordinate system string';
hdr(k).type = 'str';
hdr(k).print_format = '%s';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Data gain values';
hdr(k).short_name = 'data_gain_values';
hdr(k).pattern = 'data gain values';
hdr(k).type = 'num_arr';
hdr(k).print_format = '%f';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Data ignore value';
hdr(k).short_name = 'data_ignore_value';
hdr(k).pattern = 'data ignore value';
hdr(k).type = 'num_arr';
hdr(k).print_format = '%f';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Data offset values';
hdr(k).short_name = 'data_offset_values';
hdr(k).pattern = 'data offset values';
hdr(k).type = 'num_arr';
hdr(k).print_format = '%f';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Data reflectance gain values';
hdr(k).short_name = 'data_reflectance_gain_values';
hdr(k).pattern = 'data reflectance gain values';
hdr(k).type = 'num_arr';
hdr(k).print_format = '%f';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Data reflectance offset values';
hdr(k).short_name = 'data_reflectance_offset_values';
hdr(k).pattern = 'data reflectance offset values';
hdr(k).type = 'num_arr';
hdr(k).print_format = '%f';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Data Type';
hdr(k).short_name = 'data_type';
hdr(k).pattern = 'data type';
hdr(k).type = 'num';
hdr(k).print_format = '%s';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Default bands';
hdr(k).short_name = 'default_bands';
hdr(k).pattern = 'default bands';
hdr(k).type = 'num_arr';
hdr(k).print_format = '%f';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Default stretch';
hdr(k).short_name = 'default_stretch';
hdr(k).pattern = 'default stretch';
hdr(k).type = 'str';
hdr(k).print_format = '%s';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'DEM band';
hdr(k).short_name = 'dem_band';
hdr(k).pattern = 'dem band';
hdr(k).type = 'num';
hdr(k).print_format = '%u';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'DEM file';
hdr(k).short_name = 'dem_file';
hdr(k).pattern = 'dem file';
hdr(k).type = 'str';
hdr(k).print_format = '%s';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Description';
hdr(k).short_name = 'description';
hdr(k).pattern = 'description';
hdr(k).type = 'str_arr';
hdr(k).print_format = '%s';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'File Type';
hdr(k).short_name = 'file_type';
hdr(k).pattern = 'file type';
hdr(k).type = 'str';
hdr(k).print_format = '%s';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Full-Width-Half-Maximum (FWHM) Values';
hdr(k).short_name = 'fwhm';
hdr(k).pattern = 'fwhm';
hdr(k).type = 'num_arr';
hdr(k).print_format = '%f';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Geographic corners for non-georeferenced files';
hdr(k).short_name = 'geo_points';
hdr(k).pattern = 'geo points';
hdr(k).type = 'num_arr';
hdr(k).print_format = '%s';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Header offset';
hdr(k).short_name = 'header_offset';
hdr(k).pattern = 'header offset';
hdr(k).type = 'num';
hdr(k).print_format = '%u';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Interleave';
hdr(k).short_name = 'interleave';
hdr(k).pattern = 'interleave';
hdr(k).type = 'str';
hdr(k).print_format = '%s';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Lines';
hdr(k).short_name = 'lines';
hdr(k).pattern = 'lines';
hdr(k).type = 'num';
hdr(k).print_format = '%u';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Map information';
hdr(k).short_name = 'map_info';
hdr(k).pattern = 'map info';
hdr(k).type = 'str_arr';
hdr(k).print_format = '%s';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Major frame offsets';
hdr(k).short_name = 'major_frame_offsets';
hdr(k).pattern = 'major frame offsets';
hdr(k).type = 'num';
hdr(k).print_format = '%u';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Minor frame offsets';
hdr(k).short_name = 'minor_frame_offsets';
hdr(k).pattern = 'minor frame offsets';
hdr(k).type = 'num';
hdr(k).print_format = '%u';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Pixel size';
hdr(k).short_name = 'pixel_size';
hdr(k).pattern = 'pixel size';
hdr(k).type = 'num_arr';
hdr(k).print_format = '%f';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Product Type';
hdr(k).short_name = 'product_type';
hdr(k).pattern = 'product type';
hdr(k).type = 'str';
hdr(k).print_format = '%s';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Projection info';
hdr(k).short_name = 'projection_info';
hdr(k).pattern = 'projection info';
hdr(k).type = 'str_arr';
hdr(k).print_format = '%s';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Read procedures';
hdr(k).short_name = 'read_procedures';
hdr(k).pattern = 'read procedures';
hdr(k).type = 'str';
hdr(k).print_format = '%s';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Reflectance scale factor';
hdr(k).short_name = 'reflectance_scale_factor';
hdr(k).pattern = 'reflectance scale factor';
hdr(k).type = 'num';
hdr(k).print_format = '%f';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Rational polynomial coefficient (RPC) geolocation information';
hdr(k).short_name = 'rpc_info';
hdr(k).pattern = 'rpc info';
hdr(k).type = 'num_arr';
hdr(k).print_format = '%f';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Samples';
hdr(k).short_name = 'samples';
hdr(k).pattern = 'samples';
hdr(k).type = 'num';
hdr(k).print_format = '%u';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Security tag';
hdr(k).short_name = 'security_tag';
hdr(k).pattern = 'security tag';
hdr(k).type = 'str';
hdr(k).print_format = '%s';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Sensor type';
hdr(k).short_name = 'sensor_type';
hdr(k).pattern = 'sensor type';
hdr(k).type = 'str';
hdr(k).print_format = '%s';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Solar irradiance';
hdr(k).short_name = 'solar_irradiance';
hdr(k).pattern = 'solar irradiance';
hdr(k).type = 'num';
hdr(k).print_format = '%f';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Spectra names';
hdr(k).short_name = 'spectra_names';
hdr(k).pattern = 'spectra names';
hdr(k).type = 'str_arr';
hdr(k).print_format = '%s';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Sun azimuth';
hdr(k).short_name = 'sun_azimuth';
hdr(k).pattern = 'sun azimuth';
hdr(k).type = 'num';
hdr(k).print_format = '%f';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Sun elevation';
hdr(k).short_name = 'sun_elevation';
hdr(k).pattern = 'sun elevation';
hdr(k).type = 'num';
hdr(k).print_format = '%f';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Wavelength';
hdr(k).short_name = 'wavelength';
hdr(k).pattern = 'wavelength';
hdr(k).type = 'num_arr';
hdr(k).print_format = '%f';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Wavelength units';
hdr(k).short_name = 'wavelength_units';
hdr(k).pattern = 'wavelength units';
hdr(k).type = 'str';
hdr(k).print_format = '%s';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'X start';
hdr(k).short_name = 'x_start';
hdr(k).pattern = 'x start';
hdr(k).type = 'num';
hdr(k).print_format = '%f';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Y start';
hdr(k).short_name = 'y_start';
hdr(k).pattern = 'y start';
hdr(k).type = 'num';
hdr(k).print_format = '%f';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'Z plot average';
hdr(k).short_name = 'z_plot_average';
hdr(k).pattern = 'z plot average';
hdr(k).type = 'num_arr';
hdr(k).print_format = '%f';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'z plot range';
hdr(k).short_name = 'z_plot_range';
hdr(k).pattern = 'z plot range';
hdr(k).type = 'num_arr';
hdr(k).print_format = '%f';
hdr(k).value = [];
k = k + 1;

hdr(k).full_name = 'z plot titles';
hdr(k).short_name = 'z_plot_titles';
hdr(k).pattern = 'z plot titles';
hdr(k).type = 'str_arr';
hdr(k).print_format = '%s';
hdr(k).value = [];

%% parse header

for j = 1:size(hdr,2)

    switch hdr(j).type
        
        case 'num'
        
            pattern = [hdr(j).pattern '[\s]{0,100}[=][\s]{0,100}([0-9.-]{0,100})']; % '[\s]{0,100}[=][\s]{0,100}(\d*)'
            hdr(j).value = cellfun(@str2num, regexp(str, pattern, 'tokens', 'once'), 'UniformOutput', true);
            
        case 'num_arr'
            
            pattern = [hdr(j).pattern '[\s]{0,100}[=][\s]{0,100}[\{]([^}]+)'];
            substr = regexp(str, pattern, 'tokens', 'once');
            
            if ~isempty(substr)
                
                hdr(j).value = cellfun(@str2num, strsplit(substr{:}, ','), 'UniformOutput', true);
                
            end
            
        case 'str'
            
            pattern = [hdr(j).pattern '[\s]{0,100}[=][\s]{0,100}([^\n]+)'];
            %pattern = [hdr(j).pattern '[\s]{0,100}[=][\s]{0,100}([\w*$]+)'];
            substr = regexp(str, pattern, 'tokens', 'once');
            
            if ~isempty(substr);
                
                hdr(j).value = substr{:};
                
            end
            
        case 'str_arr'
            
            pattern = [hdr(j).pattern '[\s]{0,100}[=][\s]{0,100}[\{]([^}]+)'];
            substr = regexp(str, pattern, 'tokens', 'once');
            
            if ~isempty(substr)
                
                hdr(j).value = strsplit(substr{:}, ',');
                
                for k = 1:length(hdr(j).value)
                    
                    if ~isempty(regexp(hdr(j).value{k}, '^[0-9.\s]*$','ONCE'))
                       
                        hdr(j).value{k} = str2double(hdr(j).value(k));
                        
                    else
                        
                        hdr(j).value{k} = strtrim(hdr(j).value{k});
                        
                    end
                    
                end
                
            end
            
    end

end

%% generate key-value pairs

for j = 1:length(hdr)
    
    header_skeys.(hdr(j).short_name) = j;
    
end


%% reformat data type attribute

% The type of data representation:
% 1 = Byte: 8-bit unsigned integer
% 2 = Integer: 16-bit signed integer
% 3 = Long: 32-bit signed integer
% 4 = Floating-point: 32-bit single-precision
% 5 = Double-precision: 64-bit double-precision floating-point
% 6 = Complex: Real-imaginary pair of single-precision floating-point
% 9 = Double-precision complex: Real-imaginary pair of double precision floating-point
% 12 = Unsigned integer: 16-bit
% 13 = Unsigned long integer: 32-bit
% 14 = 64-bit long integer (signed)
% 15 = 64-bit unsigned long integer (unsigned)
 
try 
    
    data_types = {'uint8','int16','int32','float32','double','uint32','','','uint64','','','uint16','uint32','int64','uint64'};
    hdr(header_skeys.data_type).value = data_types{hdr(header_skeys.data_type).value};
    
catch
    
    warning('The mandatory Data Type attribute was not specified in the HDR file');
    
end


%% reformat for export

% remove empty fields
empty_idx = cellfun(@isempty, {hdr.value}, 'UniformOutput', 1);
hdr(empty_idx) = [];

for j = 1:size(hdr,2)
    
    info.(hdr(j).short_name) = hdr(j).value;
    
    % print info to console
    if arg.Results.verbose
        
        % fprintf(':\n', hdr(j).value, );
        if ismember(hdr(j).type,{'num','str'})
            
            fprintf(sprintf('%s: %s\\n', hdr(j).full_name, hdr(j).print_format), hdr(j).value);
            
        else
            
            fprintf('%s: too large to print\n', hdr(j).full_name);
            
        end
        
    end
    
end

