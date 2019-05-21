function ASCwrite(filepath, A, refmat, varargin)
% ASCWRITE - Writes matrix A to the ESRI ARC/INFO ASCII grid file specified with
% FILEPATH. The georeference parameters are specified with REFMAT. 
%
% Syntax:  ASCwrite(filepath, A, refmat, ...)
%
% Inputs:
%    filepath - string, path to the output ASC file
%
%    A - MxN numeric matrix, grid values
%
%    refmat - 3x2 numeric matrix, spatial referencing matrix, such that xy_map = [row, col, ones(nrows,1)] * refmat
%
%    precision (optional, default: 2) - number of decimals for the grid values
%
%    noData (optional, default: -9999) - no data value, any NaN value in A will be replaced by this value.
%
%    verbose (optional, default: true) - boolean value, verbosiy switch
%
% Example:
%    ASCwrite(mygrid.asc, ...
%        B, ...
%        refmat, ...
%        'precision', 2, ...
%        'noData', -9999);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2018b, GNU Octave 5.1.0 (configured for "x86_64-w64-mingw32")
%
% See also:
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory (LASIG)
% Website: http://mparkan.github.io/Digital-Forestry-Toolbox/
% Last revision: May 21, 2019
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund (WHFF, OFEV), project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;
addRequired(arg, 'filepath', @(x) ischar(x));
addRequired(arg, 'A', @(x) isnumeric(x));
addRequired(arg, 'refmat', @(x) all(size(x) == [3, 2]) && isnumeric(x));
addParameter(arg, 'precision', 2, @(x) isnumeric(x) && (floor(x) == x) && (numel(x) == 1));
addParameter(arg, 'noData', -9999, @(x) isnumeric(x) && (floor(x) == x) && (numel(x) == 1));
addParameter(arg, 'verbose', true, @(x) islogical(x) && (numel(x) == 1));

parse(arg, filepath, A, refmat, varargin{:});

% check filepath extension
[fpath, fname, fext] = fileparts(filepath);

if isempty(fname)
   
    error('filename is invalid');
    
end

if ~strcmp(fext, '.asc')
    
    fext = '.asc';
    
end

if isempty(fpath)
  
  filepath = sprintf('%s%s', fname, fext);
  
else
  
  filepath = sprintf('%s\\%s%s', fpath, fname, fext);
  
end


%% open file

fid = fopen(filepath, 'w+');

if fid == -1
    
    error('Could not open file at specified location')
    
end


%% write header

[nrows, ncols] = size(A);
xy_llc = [nrows+0.5, 0.5, 1] * refmat; % lower left corner coordinates

fprintf(fid, 'ncols         %.0f\n', ncols);
fprintf(fid, 'nrows         %.0f\n', nrows);
fprintf(fid, 'xllcorner     %.12f\n', xy_llc(1));
fprintf(fid, 'yllcorner     %.12f\n', xy_llc(2));
fprintf(fid, 'cellsize      %.12f\n', abs(refmat(2,1)));
fprintf(fid, 'NODATA_value  %.0f\n', arg.Results.noData);


%% write grid values

if arg.Results.verbose
    
    fprintf('writing grid to "%s"...', filepath);
    
end

% substitute no data values 
A((A == arg.Results.noData) | isnan(A)) = arg.Results.noData; 

% replicate line format
fmt = [strtrim(repmat(sprintf('%%.%uf ', arg.Results.precision), 1, size(A,2))), '\n'];

% print grid values to file
A = A';
fprintf(fid, fmt, A(:));

if arg.Results.verbose
    
    fprintf('done!\n');
    
end

%% close file

fclose(fid);
