function composite = LASmerge(inputFilePaths, varargin)
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
%    outputFilePath (optional: default: []) - The path to the merged output LAS file
%
%    updatePointSourceID (optional, default: false) - boolean value, if set
%    to true the point source ID record will be overwritten with a file
%    source index
%
%    UUID (optional, default: auto) - 32 character Univerally Unique
%    Identifier (e.g. '401db5f076cc277abc4aa217f593c48b').
%    If not specified a new UUID is automatically assigned.
%
%    verbose (optional, default: false) - boolean value, verbosiy switch
%
% Example:
%    inputPath = {'E:\data\1143424.las', 'E:\data\1143425.las'} ;
%    outputFilepath = 'E:\data\merged.las';
%    s = LASmerge(inputPath, ...
%         outputFilepath, ...
%         'updatePointSourceID', true, ...
%         'verbose', true);

% Other m-files required: LASread.m, LASwrite.m
% Subfunctions: none
% MAT-files required: none
% Compatibility: tested on Matlab R2017b, GNU Octave 4.2.1 (configured for "x86_64-w64-mingw32")
%
% See also:
%
% This code is part of the Matlab Digital Forestry Toolbox
%
% Author: Matthew Parkan, EPFL - GIS Research Laboratory
% Website: http://mparkan.github.io/Digital-Forestry-Toolbox/
% Last revision: March 28, 2018
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund (WHFF, OFEV), project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'inputFilePaths', @(x) iscell(x) && length(x) >= 1);
addOptional(arg, 'outputFilePath', [], @(x) ischar(x) || isempty(x));
addParameter(arg, 'updatePointSourceID', false, @(x) islogical(x) && (numel(x) == 1));
addParameter(arg, 'UUID', lower(strcat(dec2hex(randi(16,32,1)-1)')), @(x) ischar(x) && length(x) == 32);
addParameter(arg, 'verbose', false, @(x) islogical(x) && (numel(x) == 1));

parse(arg, inputFilePaths, varargin{:});

fclose('all');

%% check if input is structure or file path

n_parts = length(inputFilePaths);

switch class(inputFilePaths{1})
    
    case 'struct'
        
        component = inputFilePaths;
        
        if isrow(component)
            
            component = component';
            
        end
        
    case 'char'
        
        component = cell(n_parts,1);
        for j = 1:n_parts
            
            component{j,1} = LASread(inputFilePaths{j}, false, false);
            
        end
        
end


%% merge extended variable length records (EVLR)

evlr_id = cell(n_parts,1);
composite = component{1,1};

if isfield(composite, 'extended_variable_length_records')
    
    if arg.Results.verbose
        
        fprintf('merging common EVLR...');
        
    end
    
    % determine common EVLR
    for j = 1:n_parts
        
        evlr_id{j,1} = [component{j,1}.extended_variable_length_records.record_id]';
        
    end
    
    [G, ~, grp] = unique(cell2mat(evlr_id));
    idxl_common = (accumarray(grp, grp, [], @numel, nan) == n_parts);
    record_id_com = G(idxl_common);
    
    % merge common EVLR
    if any(idxl_common)
        
        n_evlr = nnz(idxl_common);
        
        for j = 1:n_parts
            
            [Lia, Locb] = ismember(record_id_com, [component{j,1}.extended_variable_length_records.record_id]);
            component{j,1}.extended_variable_length_records = component{j,1}.extended_variable_length_records(Locb(Lia));
   
        end
        
        % sort by record id
        [~, idxn_sort] = sort([composite.extended_variable_length_records.record_id]);
        composite.extended_variable_length_records = composite.extended_variable_length_records(idxn_sort);
        
        % update "value" field
        for j = 1:n_parts-1
            
            for k = 1:n_evlr
                
                composite.extended_variable_length_records(k).value = [composite.extended_variable_length_records(k).value; component{j+1,1}.extended_variable_length_records(k).value];
                composite.extended_variable_length_records(k).record_length_after_header = composite.extended_variable_length_records(k).record_length_after_header + component{j+1,1}.extended_variable_length_records(k).record_length_after_header;
                
            end
            
        end
        
    end
    
    if arg.Results.verbose
        
        fprintf('done!\n');
        
    end
    
end


%% merge point records

if arg.Results.verbose
    
    fprintf('merging point records...');
    
end

% merge records
pdr_skeys = fieldnames(composite.record);

for j = 1:length(pdr_skeys)
    
    for k = 1:n_parts-1
        
        % concatenate component 2-N to part 1
        if isfield(component{k+1,1}.record, pdr_skeys{j})
            
            
            composite.record.(pdr_skeys{j}) = [composite.record.(pdr_skeys{j}); component{k+1,1}.record.(pdr_skeys{j})];
            
        else
            
            n_entries = length(component{k+1,1}.record.x);
            composite.record.(pdr_skeys{j}) = [composite.record.(pdr_skeys{j}); zeros(n_entries,1)];
            
        end
        
    end
    
end

if arg.Results.verbose
    
    fprintf('done!\n');
    
end


%% update point source id record

if isfield(component{1,1}, 'extended_variable_length_records') && arg.Results.updatePointSourceID
    
    if arg.Results.verbose
        
        fprintf('Update point source id record...');
        
    end
    
    % update point source id record
    point_count = cellfun(@(rr) length(rr.record.x), component);
    [~, fid_record] = histc((1:sum(point_count))', cumsum([1; point_count]));
    composite.record.point_source_id = fid_record;
    
    if arg.Results.verbose
        
        fprintf('done!\n');
        
    end
    
end

%% merge variable length records (VLR)

% to be completed


%% write LAS file

if ~isempty(arg.Results.outputFilePath)
    
    las_version = composite.header.version_major * 10 + composite.header.version_minor;
    
    LASwrite(composite, ...
        arg.Results.outputFilePath, ...
        'version', las_version, ...
        'systemID', 'MERGE', ...
        'guid', arg.Results.UUID, ...
        'recordFormat', composite.header.point_data_format_id, ...
        'verbose', arg.Results.verbose);
    
end
