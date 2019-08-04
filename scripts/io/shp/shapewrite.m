## Copyright (C) 2014-2019 Philip Nienhuis
## 
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {Function File} {@var{status} =} shapewrite (@var{shpstr}, @var{fname})
## @deftypefnx {Function File} {@var{status} =} shapewrite (@var{shpstr}, @var{fname}, @var{atts})
## Write contents of map- or geostruct to a GIS shape file.
##
## @var{shpstr} must be a valid mapstruct or geostruct, a struct array with an
## entry for each shape feature, with fields Geometry, BoundingBox, and X and Y
## (mapstruct) or Lat and Lon (geostruct).  For geostructs, Lat and Lon field
## data will be written as X and Y data.  Field Geometry can have data values
## of "Point", "MultiPoint", "Line" or "Polygon", all case-insensitive.
## For each shape feature, field BoundingBox should contain the minimum and
## maximum (X,Y) coordinates in a 2x2 array [minX, minY; maxX, maxY].
## The X and Y fields should contain X (or Latitude) and Y (or Longitude)
## coordinates for each point or vertex as row vectors; for
## poly(lines) and polygons vertices of each subfeature (if present) should be
## separated by NaN entries.
##
## <Geometry>M or <Geometry>Z types (e.g., PointM, PolygonZ) can also be
## written; shapewrite.m just checks if fields "M" and/or "Z" are present in
## input mapstruct.
##
## @var{fname} should be a valid shape file name, optionally with a '.shp'
## suffix.
##
## Optional third input argument @var{atts} is one attribute name or a cellstr
## array containing a list of attribute names; only those attributes will be
## written to the .dbf file.  Alternatively a struct can be supplied with
## attibute names contained in field "FieldName" (preferrably camelcased as
## shown, but Octave treats this field's name as case-insensitive).  If
## argument @var{atts} is omitted all attributes will be written to the
## shapefile.
##
## shapewrite produces 3 files, i.e. a .shp file (the actual shape file),
## a .shx file (index file), and a .dbf file (dBase type 3) with the contents
## of additional attribute fields other than Geometry, X/Y or Lat/Lon, M, Z,
## and BoundingBox.  If no additional attributes are present, a .dbf file is
## created with the shape feature numbers as contents of field "ID".
##
## Return argument @var{status} is set to 1 if the complete shape file set was
## written successfully, to 0 otherwise.
##
## Ref: http://www.esri.com/library/whitepapers/pdfs/shapefile.pdf
##
## @seealso{shapedraw, shapeinfo, shaperead}
## @end deftypefn

## Author: Philip Nienhuis <prnienhuis@users.sf.net>
## Created: 2014-12-30

function [status] = shapewrite (shp, fname="", atts=[])

  persistent n_dbfspec = 0;
  status = 0;

  ## Input validation
  if (nargin < 1)
    print_usage;
  endif

  ## Assess shape variable type (oct or ml/geo ml/map)
  if (! isstruct (shp))
    error ("shapewrite: [map-, geo-] struct expected for argument #1");
  else
    ## Yep. Find out what type
    fldn = fieldnames (shp);
    if (ismember ("vals", fldn) && ismember ("shpbox", fldn))
      ## Assume it is an Octave-style struct read by shaperead
      otype = 0;
      warning ("shapewrite: only Matlab-type map/geostructs can be written\n");
      return;
    elseif (ismember ("Geometry", fldn) && all (ismember ({"X", "Y"}, fldn)))
      ## Assume it is a Matlab-style mapstruct
      otype = 1;
    elseif (ismember ("Geometry", fldn) && all (ismember ({"Lat", "Lon"}, fldn)))
      ## Assume it is a Matlab-style geostruct
      otype = 2;
    else
      ## Not a supported struct type
      error ("shapewrite: unsupported struct type.\n")
    endif
  endif

  ## FIXME struct field type validation

  if (strcmpi (shp(1).Geometry, "MultiPatch"))
    error ("shapewrite.m: MultiPatch type not supported");
  endif

  ## Check for dbfwrite function
  if (isempty (which ("dbfwrite")))
    error ("shapewrite.m: dbfwrite function not found. (io package installed \
and loaded?)");
    return;
  endif

  ## Check file name
  if (isempty (fname))
    error ("shapewrite: filename expected for input argument #2");
  else
    [pth, fnm, ext] = fileparts (fname);
    if (isempty (ext))
      bname = fname;
      fname = [fname ".shp"];
    ## Later on bname.shx and bname.dbf will be created
    elseif (isempty (pth))
      bname = fnm;
    else
      bname = [pth filesep fnm];
    endif
  endif

  ## Check optional 3rd argument
  if (nargin > 2)
    if (isstruct (atts))
      if (! n_dbfspec)
        warning ("shapewrite.m: DbfSpec not implemented; including requested \
attributes\n");
        n_dbfspec = 1;
      endif
      ## Get attribute names from field "FieldName"; allow lowercase and camelcase
      ## Index of case-insensitive matches of "FieldName"
      fnidx = find (ismember ("fieldname", lower (fieldnames (atts))));
      if (! isempty (fnidx))
        ## Get fieldnames out of first match
        atts = {atts.(fieldnames (atts (fnidx)){1})};
      else
        warning ("shapewrite.m: no field 'fieldname' (case-insensitive) found \
in struct\n=> input arg. #3 ignored");
      endif
    elseif (! iscellstr (atts) && ! ischar (atts))
      error ("shapewrite.m: arg.#3: attribute name or cellstr array of attribute names expected");
    endif
    ## Check if requested attributes exist at all in shapestruct
    atts = unique (atts);
    mtch = ! ismember (atts, fieldnames (shp));
    if (any (mtch))
      warning ("shapewrite.m: requested attribute(s) '%s' not in shapestruct\n", ...
               strjoin (atts(find (mtch)), "', '"));
      atts(mtch) = [];
    endif
  endif

  ## Prepare a few things
  numfeat = numel (shp);
  if (abs (otype) >= 1)
    [~, stype] = ismember(lower(shp(1).Geometry), {"point", "multipoint", "line", "polygon"});

    if (isempty (stype))
      ## Not a supported struct type
      error ("shapewrite: unsupported struct type.\n")
    else
      stype = [1, 8, 3, 5](stype);
    endif

    ## Preprocess geostructs
    if (abs (otype) == 2)
      ## Change Lat/Lon fields into X/Y
      [shp.X] = deal (shp.Lon);
      [shp.Y] = deal (shp.Lat);
    endif
    ## "Point" need not have a BoundingBox field => add a dummy if not found
    if (stype == 1 && ! ismember ("BoundingBox", fldn))
      [shp.BoundingBox] = deal ([0, 0; 0, 0]);
    endif
    if (any (ismember ({"M", "Z"}, fldn)))
      if (ismember ({"Z"}, fldn))
        stype += 10;
      else
        stype += 20;
      endif
      otype = -otype;
    endif
  endif

  ## Only now (after input checks) open .shp and .shx files & rewind just to be sure
  fids = fopen (fname, "w");
  if (fids < 0)
    error ("shapewrite: shapefile %s can't be opened for writing\n", fname);
  endif
  fseek (fids, 0, "bof");
  fidx = fopen ([ bname ".shx" ], "w");
  if (fidx < 0)
    error ("shapewrite: index file %s can't be opened for writing\n", fname);
  endif
  fseek (fidx, 0, "bof");

  ## Write headers in .shp & .shx (identical). First magic number 9994 + six
  ## zeros, the last zero in .shp is a placeholder for the yet unknown .shp
  ## file length.
  fwrite (fids, [9994 0 0 0 0 0 0], "int32", 0, "ieee-be");
  fwrite (fidx, [9994 0 0 0 0 0],   "int32", 0, "ieee-be");
  ## For .shx the file length in 16-bit words (single) is known:
  fwrite (fidx, ((numfeat * 4) + 50), "int32", 0, "ieee-be");
  ## Next, shp file version
  fwrite (fids, 1000, "int32");
  fwrite (fidx, 1000, "int32");
  ## Shape feature type
  fwrite (fids, stype, "int32");
  fwrite (fidx, stype, "int32");
  ## Bounding box. Can be run later for ML type shape structs. Fill with zeros
  fwrite (fids, [0 0 0 0 0 0 0 0], "double");
  fwrite (fidx, [0 0 0 0 0 0 0 0], "double");
  ## Prepare BoundingBox limits
  xMin = yMin = zMin = mMin =  Inf; 
  xMax = yMax = zMax = mMax = -Inf;

  ## Skip to start of first record position
%  fseek (fids, 100, "bof");
%  fseek (fidx, 100, "bof");

  ## Write shape features one by one
  if (abs (otype) >= 1)
    for ishp=1:numfeat
      ## Write record start pos to .shx file
      fwrite (fidx, ftell (fids) / 2, "int32", 0, "ieee-be");

      ## Prepare multipart polygons/lines.
      ## Find pointers to separators
      idx = [ 0 find(isnan (shp(ishp).X)) ];
      ## Eliminate trailing NaN rows
      if (isnan (shp(ishp).X))
        idx(end) = [];
      endif
      ## Augment idx for later on
      idx = unique ([ idx (numel (shp(ishp).X)+1) ]);
      ## Remove NaN separators
      idn = find (! isfinite (shp(ishp).X));
      shp(ishp).X(idn) = [];
      shp(ishp).Y(idn) = [];
      if (stype >= 10)
        shp(ishp).Z(idn) = [];
      endif
      if (stype >= 10)
        shp(ishp).M(idn) = [];
        ## M as need not be present (we assume a NaN value then).
        ## Set M to a value less than -1e39
        idm = find (! isfinite (shp(ishp).M));
        shp(ishp).M(idm) = -1.1e39;
      endif
      ## Nr. of vertices
      npt = numel (shp(ishp).X);
      ## Pointers to parts
      nptr = idx(1:end-1) .- [0:numel(idx)-2];

      ## Write record contents
      switch (stype)
        case 1                                                ## Point
          reclen = 10;
          ## Write record index number & content length (fixed)
          fwrite (fids, [ishp reclen], "int32", 0, "ieee-be");
          fwrite (fidx, reclen, "int32", 0, "ieee-be");
          ## Write shape type
          fwrite (fids, stype, "int32");
          ## Simply write XY cordinates
          fwrite (fids, [shp(ishp).X shp(ishp).Y], "double");
          ## Update overall BoundingBox
          xMin = min (xMin, shp(ishp).X);
          xMax = max (xMax, shp(ishp).X);
          yMin = min (yMin, shp(ishp).Y);
          yMax = max (yMax, shp(ishp).Y);

        case 11                                               ## PointZ
          reclen = 18;
          ## Write record index number & content length (fixed)
          fwrite (fids, [ishp reclen], "int32", 0, "ieee-be");
          fwrite (fidx, reclen, "int32", 0, "ieee-be");
          ## Write shape type
          fwrite (fids, stype, "int32");
          ## Simply write XY coordinates & Z & M value
          fwrite (fids, [shp(ishp).X shp(ishp).Y shp(ishp).Z shp(ishp).M], "double");
          ## Update overall BoundingBox
          xMin = min (xMin, shp(ishp).X);
          xMax = max (xMax, shp(ishp).X);
          yMin = min (yMin, shp(ishp).Y);
          yMax = max (yMax, shp(ishp).Y);
          zMin = min (zMin, shp(ishp).Z);
          zMax = max (zMax, shp(ishp).Z);
          mMin = min (mMin, shp(ishp).M);
          mMax = max (mMax, shp(ishp).M);

        case 21                                               ## PointM
          reclen = 14;
          ## Write record index number & content length (fixed)
          fwrite (fids, [ishp reclen], "int32", 0, "ieee-be");
          fwrite (fidx, reclen, "int32", 0, "ieee-be");
          ## Write shape type
          fwrite (fids, stype, "int32");
          ## Simply write XY coordinates & M-value
          fwrite (fids, [shp(ishp).X shp(ishp).Y shp(ishp).M], "double");
          ## Update overall BoundingBox
          xMin = min (xMin, shp(ishp).X);
          xMax = max (xMax, shp(ishp).X);
          yMin = min (yMin, shp(ishp).Y);
          yMax = max (yMax, shp(ishp).Y);
          mMin = min (mMin, shp(ishp).M);
          mMax = max (mMax, shp(ishp).M);

        case 8                                                ## MultiPoint
          reclen = (40 + 16 * numel (shp(ishp).X)) / 2;
          ## Write record index number & content length (fixed)
          fwrite (fids, [ishp reclen], "int32", 0, "ieee-be");
          fwrite (fidx, reclen, "int32", 0, "ieee-be");
          ## Write shape type (+4)
          fwrite (fids, stype, "int32");
          ## Write bounding box (+32 -> 36)
          fwrite (fids, [shp(ishp).BoundingBox'(:)]', "double");
          ## Update overall BoundingBox
          xMin = min (xMin, min (shp(ishp).X));
          xMax = max (xMax, max (shp(ishp).X));
          yMin = min (yMin, min (shp(ishp).Y));
          yMax = max (yMax, max (shp(ishp).Y));
          ## Write nr of points and XY data (+4 -> 40 + Nx16)
          fwrite (fids, numel (shp(ishp).X), "int32");
          fwrite (fids, [shp(ishp).X' shp(ishp).Y']', "double");

        case 18                                               ## MultiPointZ
          reclen = (72 + 32 * numel (shp(ishp).X)) / 2;
          ## Write record index number & content length (fixed)
          fwrite (fids, [ishp reclen], "int32", 0, "ieee-be");
          fwrite (fidx, reclen, "int32", 0, "ieee-be");
          ## Write shape type
          fwrite (fids, stype, "int32");
          ## Write bounding box
          fwrite (fids, [shp(ishp).BoundingBox'(:)]', "double");
          ## Update overall BoundingBox
          xMin = min (xMin, min (shp(ishp).X));
          xMax = max (xMax, max (shp(ishp).X));
          yMin = min (yMin, min (shp(ishp).Y));
          yMax = max (yMax, max (shp(ishp).Y));
          zMin = min (zMin, min (shp(ishp).Z));
          zMax = max (zMax, max (shp(ishp).Z));
          mMin = min (mMin, min (shp(ishp).M));
          mMax = max (mMax, max (shp(ishp).M));
          ## Write Nr of points and XY data
          fwrite (fids, numel (shp(ishp).X), "int32");
          fwrite (fids, [shp(ishp).X' shp(ishp).Y']', "double");
          ## Write Z/M range and -data in turn
          fwrite (fids, [min(shp(ishp).Z) max(shp(ishp).Z) shp(ishp).Z], "double");
          fwrite (fids, [min(shp(ishp).M) max(shp(ishp).M) shp(ishp).M], "double");

        case 28                                               ## MultiPointM
          reclen = (56 + 24 * numel (shp(ishp).X)) / 2;
          ## Write record index number & content length (fixed)
          fwrite (fids, [ishp reclen], "int32", 0, "ieee-be");
          fwrite (fidx, reclen, "int32", 0, "ieee-be");
          ## Write shape type
          fwrite (fids, stype, "int32");
          ## Write bounding box
          fwrite (fids, [shp(ishp).BoundingBox'(:)]', "double");
          ## Update overall BoundingBox
          xMin = min (xMin, min (shp(ishp).X));
          xMax = max (xMax, max (shp(ishp).X));
          yMin = min (yMin, min (shp(ishp).Y));
          yMax = max (yMax, max (shp(ishp).Y));
          mMin = min (mMin, min (shp(ishp).M));
          mMax = max (mMax, max (shp(ishp).M));
          ## Write Nr of points and XY data
          fwrite (fids, numel (shp(ishp).X), "int32");
          fwrite (fids, [shp(ishp).X' shp(ishp).Y']', "double");
          ## Write M range and M data
          fwrite (fids, [min(shp(ishp).M) max(shp(ishp).M) shp(ishp).M], "double");

        case {3, 5}                                           ## Polyline/-gon
          ## Content length
          reclen = (44 + (numel(idx)-1)*4 + 2*8*npt) / 2;
          ## Write record index number & content length (fixed)
          fwrite (fids, [ishp reclen], "int32", 0, "ieee-be");
          fwrite (fidx, reclen, "int32", 0, "ieee-be");
          ## Write shape type
          fwrite (fids, stype, "int32");
          ## Write bounding box
          fwrite (fids, [shp(ishp).BoundingBox'(:)]', "double");
          ## Update overall BoundingBox
          xMin = min (xMin, min (shp(ishp).X));
          xMax = max (xMax, max (shp(ishp).X));
          yMin = min (yMin, min (shp(ishp).Y));
          yMax = max (yMax, max (shp(ishp).Y));
          ## Write number of parts, number of points, part pointers
          fwrite (fids, [(numel(idx)-1) npt nptr ], "int32");
          fwrite (fids, [shp(ishp).X' shp(ishp).Y']'(:), "double");

        case {13, 15}                                         ## Polyline/-gonZ
          ## Content length
          reclen = (76 + (numel(idx)-1)*4 + 4*8*npt) / 2;
          ## Write record index number & content length (fixed)
          fwrite (fids, [ishp reclen], "int32", 0, "ieee-be");
          fwrite (fidx, reclen, "int32", 0, "ieee-be");
          ## Shape type
          fwrite (fids, stype, "int32");
          ## Write bounding box
          fwrite (fids, [shp(ishp).BoundingBox'(:)]', "double");
          ## Update overall BoundingBox
          xMin = min (xMin, min (shp(ishp).X));
          xMax = max (xMax, max (shp(ishp).X));
          yMin = min (yMin, min (shp(ishp).Y));
          yMax = max (yMax, max (shp(ishp).Y));
          zMin = min (zMin, min (shp(ishp).Z));
          zMax = max (zMax, max (shp(ishp).Z));
          mMin = min (mMin, min (shp(ishp).M));
          mMax = max (mMax, max (shp(ishp).M));
          ## Write number of parts, number of points, part pointers
          fwrite (fids, [(numel(idx)-1) npt nptr ], "int32");
          ## Write XY data
          fwrite (fids, [shp(ishp).X' shp(ishp).Y']'(:), "double");
          fwrite (fids, [min(shp(ishp).Z) max(shp(ishp).Z) ...
                         shp(ishp).Z], "double");
          fwrite (fids, [min(shp(ishp).M) max(shp(ishp).M) ...
                         shp(ishp).M], "double");

        case {23, 25}                                         ## Polyline/-gonM
          ## Content length
          reclen = (60 + (numel(idx)-1)*4 + 3*8*npt) / 2;
          ## Write record index number & content length (fixed)
          fwrite (fids, [ishp reclen], "int32", 0, "ieee-be");
          fwrite (fidx, reclen, "int32", 0, "ieee-be");
          ## Write shape type
          fwrite (fids, stype, "int32");
          ## Write bounding box
          fwrite (fids, [shp(ishp).BoundingBox'(:)]', "double");
          ## Update overall BoundingBox
          xMin = min (xMin, min (shp(ishp).X));
          xMax = max (xMax, max (shp(ishp).X));
          yMin = min (yMin, min (shp(ishp).Y));
          yMax = max (yMax, max (shp(ishp).Y));
          mMin = min (mMin, min (shp(ishp).M));
          mMax = max (mMax, max (shp(ishp).M));
          ## Write number of parts, number of points, part pointers
          fwrite (fids, [(numel(idx)-1) npt nptr ], "int32");
          ## Write XY data
          fwrite (fids, [shp(ishp).X' shp(ishp).Y']'(:), "double");
          ## Write M range and M data
          fwrite (fids, [min(shp(ishp).M) max(shp(ishp).M) ...
                         shp(ishp).M], "double");

        otherwise
          ## Future shape types or types unsupported yet (MultiPatch)

      endswitch
    endfor
  endif

  ## Write file length and overall BoundingBox into .shp header
  flen = ftell (fids);
  fseek (fids, 24, "bof");
  fwrite (fids, flen/2, "int32", 0, "ieee-be");
  fseek (fids, 36, "bof");
  fwrite (fids, [xMin yMin xMax yMax], "double");
  ## Same for .shx header
  xlen = ftell (fidx);
  fseek (fidx, 24, "bof");
  fwrite (fidx, xlen/2, "int32", 0, "ieee-be");
  fseek (fidx, 36, "bof");
  fwrite (fidx, [xMin yMin xMax yMax], "double");
  if (stype > 10)
    ## +-Inf & NaN not allowed in shapefiles
    zm = [zMin zMax mMin mMax];
    zm (! isfinite (zm)) = -1e-39;
    fwrite (fids, zm, "double");
    fwrite (fidx, zm, "double");
  endif

  ## Close files
  fclose (fids);
  fclose (fidx);

  ## Write .dbf file.
  ## Remove basic attributes
  if (abs (otype) == 1)
    ## Attributes + shp data in mapstruct
    shp = rmfield (shp, {"Geometry", "BoundingBox", "X", "Y"});
  elseif (abs (otype) == 2)
    ## Attributes + shp data in geostruct
    shp = rmfield (shp, {"Geometry", "BoundingBox", "Lat", "Lon", "X", "Y"});
  endif
  if (otype < 1)
    shp = rmfield (shp, {"M"});
    if (isfield (shp, "Z"))
      shp = rmfield (shp, {"Z"});
    endif
  endif

  ## Write rest of attributes
  if (nargin == 3)
    ## Only write user-specified attribute selection
    fldn = fieldnames (shp);
    ## First remove regular attributes (with a value for each vertex)
    fldn (ismember (fldn, {"Geometry", "BoundingBox", "Lat", "Lon", "X", "Y", "M", "Z"})) = [];
    ## Next, only retain user-specified attributes
    shp = rmfield (shp, setdiff (fldn, atts));
  endif
  attribs = cell (numfeat + 1, numel (fieldnames (shp)));
  if (! isempty (attribs))
    attribs(1, :) = fieldnames (shp);
    attribs(2:end, :) = (squeeze (struct2cell (shp)))';
  else
    ## Substitute ID attribute
    attribs{1, 1} = "ID";
    [attribs{2:end}] = deal (num2cell ([1:size(shp, 2)]){:});
  endif
  try
    status = dbfwrite ([ bname ".dbf"], attribs);
    status = 1;
  catch
    warning ("shapewrite: writing attributes to file %s failed\n", [bname ".dbf"]);
  end_try_catch

endfunction


## Test various geometries: (1) Point
%!test
%! shp.Geometry = "Point";
%! shp.X = 10;
%! shp.Y = 20;
%! shp.Z = 30;
%! shp.M = -1;
%! shp.attr_1 = "Attribute1";
%! shp.attr_Z = "AttributeA";
%! shp(2).Geometry = "Point";
%! shp(2).X = 11;
%! shp(2).Y = 25;
%! shp(2).Z = 32;
%! shp(2).M = -2;
%! shp(2).attr_1 = "Attribute2";
%! shp(2).attr_Z = "AttributeB";
%! fname = tempname ();
%! sts = shapewrite (shp, fname);
%! assert (sts, 1, eps);
%! ## Check index file
%%! fx = fopen ([fname ".shx"], "r");
%%! fseek (fx, 100, "bof");
%%! shxinfo = fread (fx, "*int32", 0, "ieee-be");
%%! assert (shxinfo, int32 ([50; 36; 72; 36]));
%%! fclose (fx);
%! ## Check on filesizes, based on Esri shapewrite doc
%! assert (stat ([fname ".shp"]).size, 188, 1e-10);
%! assert (stat ([fname ".shx"]).size, 116, 1e-10);
%! shp2 = shaperead ([fname ".shp"]);
%! assert (size (shp2), [2 1]);
%! flds = fieldnames (shp2);
%! fields = {"Geometry", "X", "Y", "attr_1", "attr_Z"};
%! ism = ismember (fields, flds);
%! ## Do we have only those fields?
%! assert (numel (ism), numel (fields));
%! ## Do we have only those fields?
%! assert (sum (ism), numel (ism));
%! unlink ([fname ".shp"]);
%! unlink ([fname ".shx"]);
%! assert ([shp2.X shp2.Y], [10 11 20 25], 1e-10);
%! assert ({shp2.Geometry}, {"Point", "Point"});
%! assert ({shp2.attr_1}, {"Attribute1", "Attribute2"});

## Test various geometries: (2) Line & Polygon
%!test
%! shp.Geometry = "Line";
%! shp.BoundingBox = [9 110; 19 120];
%! shp.X = [10 110 NaN   9 109];
%! shp.Y = [20 120 NaN  19 119];
%! shp.Z = [30 130 NaN  29 129];
%! shp.M = [-1   1 NaN -11  11];
%! shp.attr_1 = "Attribute1";
%! shp.attr_Z = "AttributeA";
%! shp(2).Geometry = "Line";
%! shp(2).BoundingBox = [11 211; 24 225];
%! shp(2).X = [11 111 211 NaN  11 110 NaN  12  200];
%! shp(2).Y = [25 125 225 NaN  24 124 NaN  26  200];
%! shp(2).Z = [32 132 232 NaN  31 131 NaN  33  200];
%! shp(2).M = [-2 NaN  -3 NaN -22  22 NaN -33   33];
%! shp(2).attr_1 = "Attribute2";
%! shp(2).attr_Z = "AttributeB";
%! fname = tempname ();
%! sts = shapewrite (shp, fname);
%! assert (sts, 1, eps);
%! ## Check index file
%%! fx = fopen ([fname ".shx"], "r");
%%! fseek (fx, 100, "bof");
%%! shxinfo = fread (fx, "*int32", 0, "ieee-be");
%%! assert (shxinfo, int32 ([50; 106; 160; 156]));
%%! fclose (fx);
%! ## Check on filesizes, based on Esri shapewrite doc
%! assert (stat ([fname ".shp"]).size, 640, 1e-10);
%! assert (stat ([fname ".shx"]).size, 116, 1e-10);
%!
%! ## Check Matlab-style reading
%! shp2 = shaperead ([fname ".shp"]);
%! assert (size (shp2), [2 1]);
%! flds = fieldnames (shp2);
%! fields = {"Geometry", "BoundingBox", "X", "Y", "attr_1", "attr_Z"};
%! ism = ismember (fields, flds);
%! ## Do we have only those fields?
%! assert (numel (ism), numel (fields));
%! ## Do we have only those fields?
%! assert (sum (ism), numel (ism));
%! assert ([shp2.BoundingBox], [9 110 11 211; 19 120 24 225], 1e-10);
%! assert ({shp2.attr_Z}, {"AttributeA", "AttributeB"});
%!
%! ## Re-read using M & Z values, read only 2nd feature
%! shp2 = shaperead ([fname ".shp"], "e", "rec", 2);
%! assert (size (shp2), [1 1]);
%! flds = fieldnames (shp2);
%! fields = {"Geometry", "BoundingBox", "X", "Y", "Z", "M", "attr_1", "attr_Z"};
%! ism = ismember (fields, flds);
%! ## Do we have only those fields?
%! assert (numel (ism), numel (fields));
%! ## Do we have only those fields?
%! assert (sum (ism), numel (ism));
%! ## Check regular numerical values
%! assert ([shp2.X; shp2.Y; shp2.Z], ...
%!         [11, 111, 211, NaN, 11, 110, NaN, 12, 200; ...
%!          25, 125, 225, NaN, 24, 124, NaN, 26, 200; ...
%!          32, 132, 232, NaN, 31, 131, NaN, 33, 200], 1e-10);
%! ## Check missing M value
%! assert (shp2.M, [-2, NaN, -3, NaN, -22, 22, NaN, -33, 33], 1e-10);
%! assert ({shp2.attr_1, shp2.attr_Z}, {"Attribute2", "AttributeB"});
%!
%! unlink ([fname ".shp"]);
%! unlink ([fname ".shx"]);
