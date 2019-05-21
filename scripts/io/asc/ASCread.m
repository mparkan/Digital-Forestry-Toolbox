function [A, refmat] = ASCread(filepath, varargin)
% ASCREAD - Reads an ESRI ARC/INFO ASCII grid file specified with
% FILEPATH. The georeference parameters are returned in REFMAT. 
%
% Syntax:  [A, refmat] = ASCread(filepath, ...)
%
% Inputs:
%    filepath - string, path to the input ASC file
%
%    verbose (optional, default: true) - boolean value, verbosiy switch
%
% Outputs:
%    A - MxN numeric matrix, grid values
%
%    refmat - 3x2 numeric matrix, spatial referencing matrix, such that xy_map = [row, col, ones(nrows,1)] * refmat
%
% Example:
%    [A, refmat] = ASCread(mygrid.asc);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2017b, GNU Octave 4.4.1 (configured for "x86_64-w64-mingw32")
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
addParameter(arg, 'verbose', true, @(x) islogical(x) && (numel(x) == 1));

parse(arg, filepath, varargin{:});

% check filepath extension
[~, fname, fext] = fileparts(filepath);

if isempty(fname)
   
    error('filename is invalid');
    
end

if ~ismember(fext, {'.asc', '.grd'})
    
    error('extension is invalid');
    
end


%% open file

fid = fopen(filepath, 'r');

if fid == -1
    
    error('Could not open file at specified location')
    
end


%% read header

ncols = sscanf(lower(fgetl(fid)), 'ncols %u');
nrows = sscanf(lower(fgetl(fid)), 'nrows %u');
xllcorner = sscanf(lower(fgetl(fid)), 'xllcorner %f');
yllcorner = sscanf(lower(fgetl(fid)), 'yllcorner %f');
cellsize = sscanf(lower(fgetl(fid)), 'cellsize %f');
nodata = sscanf(lower(fgetl(fid)), 'nodata_value %f');

% create spatial reference matrix
dx = cellsize;
dy = -cellsize;
refmat = [0, dy; dx, 0; xllcorner-dx/2, yllcorner+(nrows*cellsize)-dy/2];
% xy_map = [1, 1, 1] * refmat; % check

if arg.Results.verbose
    
    fprintf('ncols: %u\n', ncols);
    fprintf('nrows: %u\n', nrows);
    fprintf('xllcorner: %f\n', xllcorner);
    fprintf('yllcorner: %f\n', yllcorner);
    fprintf('cellsize: %f\n', cellsize);
    fprintf('nodata: %f\n', nodata);
    
end


%% read grid values

if arg.Results.verbose
    
    fprintf('reading grid...');
    
end

A = fscanf(fid, '%f', [nrows, ncols])';

% substitute no data values with NAN
A(A == nodata) = nan;

if arg.Results.verbose
    
    fprintf('done!\n');
    
end

%% close file

fclose(fid);
