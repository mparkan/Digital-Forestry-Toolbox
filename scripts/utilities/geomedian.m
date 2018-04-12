function [M, I, SSD] = geomedian(X)
% GEOMEDIAN - computes the geometric median M of points X. Optionally
% outputs the index I and Sum of Squarred Distances SSD of the geometric
% median point
% 
% Syntax:  [M, I] = geomedian(xyz, ...)
%
% Inputs:
%    X - UxV numeric matrix, point coordinates
%
%    fig (optional, default: false) - boolean value, switch to plot figures
%
% Outputs:
%    M - 1xV numeric vector, geometric median coordinates
%
%    I - scalar, index of the geometric median point
%
%    SSD - scalar, sum of squarred distances from the the gemoetric median
%    point to all other points
%
% Example:
%
%    [M, I] = geomedian(X);
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
% Last revision: April 12, 2018
% Acknowledgments: This work was supported by the Swiss Forestry and Wood Research Fund, WHFF (OFEV) - project 2013.18
% Licence: GNU General Public Licence (GPL), see https://www.gnu.org/licenses/gpl.html for details


%% check argument validity

arg = inputParser;

addRequired(arg, 'X', @(x) isnumeric(x));

parse(arg, X);


%% compute pairwise squared distances between all points

if size(X,1) >= 2
    
    % find unique combinations of all point pairs
    pairs = nchoosek(1:size(X, 1), 2);
    
    % compute pairwise squared distances between all unique points
    dist = sum((X(pairs(:, 1),:) - X(pairs(:, 2),:)) .^ 2, 2);
    PAIRS = [pairs; pairs(:,2:-1:1)];
    DIST = [dist; dist];
    
    % compute the sum of squarred distances for each point
    SSD = accumarray(PAIRS(:,1), DIST, [], @(x) sum(x));
    
    % find the point that minimizes the sum of squarred distances to all others
    [SSD, I] = min(SSD);
    M = X(I,:);
    
else
    
    I = 1;
    M = X;
    
end

