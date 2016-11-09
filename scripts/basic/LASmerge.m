function s = LASmerge(inputFilePaths, varargin)
% LASMERGE - Merge multiple LAS files or LAS structures
%
% LASMERGE(INPUTFILEPATHS, ...) Merge multiple LAS files or LAS structures. Optionaly writes the
% clipped point cloud to the LAS file specified by OUTPUTFILEPATH
%
% Syntax:  LASmerge(inputPath, ...)
%
% Inputs:
%    inputPath - A cell array of LAS structures or a cell array of paths to LAS files
%
%    outputFilePath (optional: default: []) - The path to the clipped output LAS file
%
%    verbose (optional, default: false) - boolean value, verbosiy switch
%
% Example:
%    inputPath = {'E:\data\1143424.las', 'E:\data\1143425.las'} ;
%    outputFilepath = 'E:\data\merged.las';
%    s = LASmerge(inputPath, outputFilepath, 'verbose', true);
%
% Other m-files required: LASread.m, LASwrite.m
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2016a
%
% See also:
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory
% Website: http://lasig.epfl.ch/
% Last revision: November 11, 2016
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund (WHFF, OFEV), project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'inputFilePaths', @(x) iscell(x) && length(x) >= 1); % && length(x) >= 2 && length(unique(cellfun(@class, x, 'UniformOutput', false))) == 1
addOptional(arg, 'outputFilePath', [], @(x) ischar(x) || isempty(x));
addParameter(arg, 'verbose', false, @(x) islogical(x) && (numel(x) == 1));

parse(arg, inputFilePaths, varargin{:});


%% check if input is structure or file path

n_parts = length(inputFilePaths);

switch class(inputFilePaths{1})
    
    case 'struct'
        
        parts = inputFilePaths;
        
        if isrow(parts)
            
            parts = parts';
            
        end
        
    case 'char'
        
        parts = cell(n_parts,1);
        for j = 1:n_parts
            
            parts{j,1} = LASread(inputFilePaths{j}, false, arg.Results.verbose);
            
        end
        
end


%% merge point parts

if arg.Results.verbose
    
    fprintf('merging point clouds...');
    
end

pdr_skeys = fieldnames(parts{1,1}.record);
    
for i = 1:length(pdr_skeys)
    
    for k = 1:n_parts-1
        
        % concatenate parts 2-N to part 1
        parts{1,1}.record.(pdr_skeys{i}) = [parts{1,1}.record.(pdr_skeys{i}); parts{k+1,1}.record.(pdr_skeys{i})];
        
    end
    
end

if arg.Results.verbose
    
    fprintf('done!\n');
    
end

s = parts{1,1};


%% write LAS file

if ~isempty(arg.Results.outputFilePath)

        las_version = s.header.version_major * 10 + s.header.version_minor;
        
        LASwrite(s, arg.Results.outputFilePath, ...
        'version', las_version, ...
        'systemID', 'MERGE', ...
        'recordFormat', s.header.point_data_format_id, ...
        'verbose', arg.Results.verbose);

end