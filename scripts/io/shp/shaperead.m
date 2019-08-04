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
## @deftypefn {Function File}  [@var{outstruct} ] = shaperead (@var{shp_filename})
## @deftypefnx {Function File} [@var{outstruct} ] = shaperead (@var{shp_filename}, @var{outstyle})
## @deftypefnx {Function File} [@var{outstruct} ] = shaperead (@var{shp_filename}, @var{outstyle}, @var{opts})
## @deftypefnx {Function File} [@var{outstruct}, @var{atts} ] = shaperead (@var{shp_filename}, ...)
## Read an ArcGis shapefile set (shp, shx and dbf).
##
## Depending on the value of @var{outstyle} some different output formats
## will be returned:
##
## @table @code
## @item  0 (numeric)
## @itemx  ml (case-insensitive)
## @itemx  m
## Return a Matlab-compatible M X 1 struct with a separate entry for each shape
## feature in the shape file.  Each struct element contains fields "Geometry"
## (shape type), "BoundingBox" ([minX minY ; maxX maxY]), X, Y (coordinates of
## points in the shape item as row vectors).  For multi-part items, the
## coordinates of each part are separated by NaNs.  This output format supports
## neither M and Z type nor MultiPatch shape features.  For M and Z type shape
## features the M and Z values will simply be ignored.
## The struct is augmented with attributes found in the accompanying .dbf file,
## if found.
##
## For ML-style output, if only one output argument is requested the attributes
## in the .dbf file will be augmented to that struct.  If two output arguments
## are requested, the attributes will be returned separately in output struct
## @var{atts}.  
##
## @item  1 (numeric)
## @itemx  ext (case-insensitive)
## @itemx  e
## Same as 1 but M and Z type and MultiPatch shape features are accepted.  The
## resulting output struct is no more ML-compatible.  If the shapefile contains
## M and/or Z type shape features the mapstruct or geostruct has extra fields M
## and -optionally- Z.  Note that MultiPatch shape features may not have
## M-values even if Z-values are present.  For MultiPatch shapes another field
## Parts is added, a Px2 array with zero-based indices to the first vertex of
## each subfeature in the XYZ fields in column 1 and the type of each
## subfeature in column 2; P is the number of shape feature parts.
##
## @item  2 (numeric)
## @itemx  oct (case-insensitive)
## @itemx  o
## Return a struct containing a N X 6 double array "vals" containing the X, Y,
## and Z coordinates, M-values, record nr. and type of each point in the shape
## file.  If no M or Z values were present the relevant columns contain
## NaNs.  Individual shape features and shape parts are separated by a row of
## NaN values.  The field "idx" contains 1-based pointers into field vals to
## the first vertex of each shape feature.
## Field "bbox" contains an 8 X M double array of XYZ coordinates of the
## bounding boxes and min/max M-values corresponding to the M items found in
## the .shp file; for point shapes these contain NaNs.  
## Field "npt" contains a 1 X M array of the number of points for each item.
## Field "npr" contains a 1 X M cell array containing a row of P part indices
## (zero-based) for each Polyline, Polygon or MultiPatch part in the shape
## file; for multipatch each cell contains another row with the part types;
## for other item types (point etc.) the cell array contains empty rows.
## A separate field "shpbox" contains the overall bounding box  X, Y and Z
## coordinates and min/max M-values in a 4 X 2 double array.  If the shape file
## contains no Z or M values the corresponding columns are filled with NaNs.
##
## The struct field "fields" contains a cellstr array with names of the columns.
## If a corresponding .dbf file was read, the struct array also contains
## a field for each attribute found in the .dbf file with the corresponding
## field name, each containing a 1 X M array of attribute values matching the
## M items in the .shp file.  These arrays can be double, char or logical,
## depending on the type found in the .dbf file.
##
## @item  3 (numeric)
## @itemx  dat (case-insensitive)
## @itemx  d
## Same as OCT or 0 but without a row of NaN values between each shape
## file item in the VALS array.
## @end table
##
## If a character option is given, just one character will suffice. The default
## for @var{outstyle} is "ml".
##
## The output of 'shaperead' can be influenced by property-value pairs. The
## following properties are recognized (of which only the first three
## characters are significant, case doesn't matter):
##
## @table @code
## @item Attributes
## Normally all attributes associated with the shape features will be read
## and returned in the output struct(s).  To limit this to just some
## attributes, enter a value consisting of a cell array of attribute names to
## be read.  To have no attributes read at all, specify @{@}, an empty cell
## array.
##
## @item BoundingBox
## Select only those shape items (features) whose bounding box lies within, or
## intersets in at least one point with the limits of the BoundingBox value (a
## 2 X 2 double array [Minx, MinY; MaxX, MaxY]).
## No intersection or clipping with the BoundingBox value will be done by
## default!
##
## @item Clip
## (only useful in conjuction with the BoundingBox property) If a value of 1
## or true is supplied, clip all shapes to the bounding box limits.  This
## option may take quite a bit of processing time.  If a value of "0" or false
## is given, do not perform clipping.  The default value is 0.
## Clipping is merely meant to be performed in the XY plane.  Clipping 3D
## shapes is supported but may lead to unexpected results.
## For Z and M type polylines and polygons including MultiPatch and
## Longitude/Latitude/Height types, Z (Height) and M values for each vertex in
## the clipped shape feature are simply copied over from the nearest vertex in
## the original shape feature.  This implies that Z and M values of new
## vertices created on the bounding box edges may be less optimal.
##
## For clipping polylines and polygons the Octave-Forge geometry package needs
## to be installed and loaded.
## 
## @item Debug
## If a value of 'true' or 1 is given, shaperead echoes the current record
## number while reading.  Can be useful for very big shapefiles.  The default
## value is 0 (no feedback).  If a Matlab-compatible output structarray is
## requested and the Bounding Box property is specified, the extracted shape
## feature indices are added to the field "___Shape_feature_nr___".
##
## @item RecordNumbers
## Select only those records whose numbers are listed as integer values in
## an array following RecordNumbers property. Neither the size nor the class of
## the array matters as long as it is a numeric array.
##
## @item UseGeoCoords
## (Only applicable if a Matlab-style output struct is requested). If a value
## of 'true' (or 1) is supplied, return a geostruct rather than a mapstruct.
## If a value of 0 or false is given, return a mapstruct.
## The mere difference is that in a geostruct the fields 'X' and 'Y' (and
## optionally 'Z') are replaced by 'Long' and 'Lat' (and 'Hght').  The
## default value is 'false' (return a mapstruct').
## @end table
##
## Ref: http://www.esri.com/library/whitepapers/pdfs/shapefile.pdf
##
## @seealso{geoshow, mapshow, shapedraw, shapeinfo}
## @end deftypefn

## Author: Philip Nienhuis <prnienhuis@users.sf.net>
## Created: 2014-11-07

function [ outs, oatt ] = shaperead (fname, varargin);

## FIXME Implementation needed of these ML input arguments:
##       - Selector  (supposedly difficult. I'd prefer myself to leave this
##                   outside of shaperead, it needlessly complicates and slows
##                   this function /PRN)

  ## Check input
  if (nargin < 1)
    print_usage ();
  endif
   
  ## Check file name
  [pth, fnm, ext] = fileparts (fname);
  if (isempty (ext))
    bname = fname;
    fname = [fname ".shp"];
  elseif (isempty (pth))
    ## Later on bname.shx and bname.dbf will be read
    bname = fnm;
  else
    ## Later on bname.shx and bname.dbf will be read
    bname = [pth filesep fnm];
  endif

  ## Find out what args have been supplied
  if (nargin == 1)
    ## Only filename supplied. Set "ml" (Matlab) type as default
    outopts = 0;
  elseif (nargin == 2)
    ## Assume filename + outopts was supplied
    if (isempty (varargin{1}))
      ## Assume ML-style output
      outopts = 0;
    else
      outopts = varargin{1};
    endif
    varargin = {};
  elseif (rem (nargin, 2) == 0)
    ## Even number of input args => Outopts & at least one pair of varargin
    outopts = varargin{1};
    varargin(1) = [];
  else
    ## Odd nr; maybe only filename and prop/val(s) supplied, outstyle skipped
    if (ischar (varargin{1}))
      ## Check arg#2, must be a property name then
      if (! ismember (lower (varargin{1}(1:min(3, numel (varargin{1})))), ...
                          {"att", "bou", "cli", "deb", "rec", "sel", "use"}))
        error ("shaperead: property name expected for arg. #2");
      endif
    else
      ## no outstyle or property => wrong input
      print_usage ();
    endif
    outopts = 0;
  endif

  ## Check output type arg
  if (isnumeric (outopts))
    if (outopts < 0 || outopts > 3)
      error ("shaperead: arg. #2 integer value out of range 0-3\n");
    endif
  elseif (ischar (outopts))
    outopts = lower (outopts);
    if (! any (strncmp (outopts, {"ml", "ext", "oct", "dat"}, 1)))
      error ("shaperead: arg. #2 char value should be one of 'ml', 'ext', \
'oct' or 'dat'\n");
    endif
    switch outopts
      case {"ml", "m"}
        outopts = 0;
      case {"ext", "e"}
        outopts = 1;
      case {"oct", "o"}
        outopts = 2;
      case {"dat", "d"}
        outopts = 3;
      otherwise
        error ("shaperead: illegal value for arg. #2: '%s' - expected 'ml', \
'ext', 'oct' of 'dat'", ...
               outopts);
    endswitch
  else
    error ("shaperead: numeric or character type expected for arg. #2\n");
  endif

  ## Check .shp file existence
  fidp = fopen (fname, "r");
  if (fidp < 0)
    error ("shaperead: can't open file %s\n", fname);
  endif
  ## Temporarily close to avoid file handle leaks during further input checks
  fclose (fidp);
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Open .shx file to help speed up seeks to next records. We need this info
  ## for s_recs check below
  have_shx = 0;
  fidx = fopen ([bname ".shx"], "r");
  if (fidx < 0)
    warning ("shaperead: index file %s not found\n", [fnm ".shx"]);
    s_recs = [];
  else
    fseek (fidx, 24, "bof");
    fxlng = fread (fidx, 1, "int32", 72, "ieee-be");
    nrec = (fxlng - 50) / 4;
    fseek (fidx, 100, "bof");
    ## Get record start positions & -lengths in 16-bit words
    ridx = reshape (fread (fidx, nrec*2, "int32", 0, "ieee-be"), 2, [])';
    fclose (fidx);
    ## Get indices & lengths in bytes
    ridx *= 2;
    have_shx = 1;
    s_recs = [1 : size(ridx, 1)];
  endif

  ## Parse options; first set defaults
  clip = 0;
  dbug   = 0;
  s_atts = [];
  s_bbox = [];
  s_geo  = 0;
  ## Init collection of records that meet BB criteria
  bb_union = [];
  ## Process input args
  for ii = 1:2:length (varargin)
    if (! ischar (varargin{ii}))
      error ("shaperead: property %d: property name expected but got a %s value", ...
              (ii+1)/2, class (varargin{ii}));
    elseif (numel (varargin{ii}) < 3)
      warning ("shaperead: unknown option '%s' - ignored\n", varargin{ii});
    else
      switch (lower (varargin{ii})(1:3))
        case "att"
          ## Select records based on attribute values
          s_atts = varargin{ii+1};
        case "bou"
          ## Select whether record/shape features partly lie inside or outside limits
          try
            s_bbox = double (varargin{ii+1});
            if (numel (s_bbox) != 4)
              error ("shaperead: 2 X 2 numeric array expected for BoundingBox\n");
            endif
          catch
            error ("shaperead: numeric 2 X 2 array expected for BoundingBox\n");
          end_try_catch
          ## Initialize supplementary polygon array here
          sbox = [s_bbox(1) s_bbox(3); s_bbox(2) s_bbox(3); s_bbox(2) s_bbox(4); ...
                  s_bbox(1) s_bbox(4); s_bbox(1) s_bbox(3)];
        case "cli"
          ## Clip polygons to requested BoundingBox
          try
            clip = logical (varargin{ii+1});
          catch
            error ("numeric or logical value expected for 'Clip'\n");
          end_try_catch
        case "deb"
          ## Set verbose output (count records)
          try
            dbug = logical (varargin{ii+1});
          catch
            error ("numeric or logical value expected for 'Debug'\n");
          end_try_catch
        case "rec"
          ## Select record nrs directly. Check for proper type & clean up
          try
            s_recs = sort (unique (double (varargin{ii+1}(:))));
            if (have_shx && any (s_recs > nrec))
              printf ("shaperead: requ. record nos. > nr. of records (%d) ignored\n", ...
                      nrec);
              s_recs (find (srecs > nrec)) = [];
            endif
          catch
            error ("shaperead: numeric value or array expected for RecordNumbers\n");
          end_try_catch
        case "sel"
          ## A hard one, to be implemented later?
          printf ("shaperead: 'Selector' option not implemented, option ignored\n");
        case "use"
          ## Return a geostruct or a mapstruct (default). Only for ML-structs
          if (outopts != 0)
            error ("shaperead: UseGeoCoords only valid for Matlab-style output\n");
          endif
          try
            s_geo = logical (varargin{ii+1});
          catch
            error ("shaperead: logical value type expected for 'UseGeoCoords'\n");
          end_try_catch
        otherwise
          warning ("shaperead: unknown option '%s' - ignored\n", varargin{ii});
      endswitch
    endif
  endfor
  ## Post-processing
  if (clip)
    if (isempty (s_bbox))
      warning ("shaperead: no BoundingBox supplied => Clip option ignored.\n");
      clip = 0;
    endif
    if (isempty (which ("clipPolygon_clipper")))
      ## No OF geometry package?
      printf  ("shaperead: function 'clipPolygon' not found. Clip option ignored\n");
      warning ("           (OF geometry package installed and loaded?)\n");
      clip = 0;
    endif
    if (isempty (which ("distancePoints")))
      ## No OF geometry package?
      printf  ("shaperead: function 'distancePoints' not found. Clip option ignored\n");
      warning ("           (OF geometry package installed and loaded?)\n");
      clip = 0;
    endif
  endif

  ## Open .shp file
  fidp = fopen (fname, "r");

  ## ============= Preparations done, now we can start reading ============
  ## ---------------------- 2. Read .shp file proper ----------------------

  ## Start reading header
  fseek (fidp, 0, "bof");

  ## Read & check file code
  fcode = fread (fidp, 1, "int32", 20, "ieee-be");
  if (fcode != 9994)
    error ("%s is not a valid shapefile\n", fname);
  endif
  flngt = fread (fidp, 1, "int32", 0, "ieee-be") * 2;
  fvsn  = fread (fidp, 1, "int32", 0, "ieee-le");
  shpt  = fread (fidp, 1, "int32");                         ## Shape file type

  shpbox.X(1) = fread (fidp, 1, "double");
  shpbox.Y(1) = fread (fidp, 1, "double");
  shpbox.X(2) = fread (fidp, 1, "double");
  shpbox.Y(2) = fread (fidp, 1, "double");
  shpbox.Z(1) = fread (fidp, 1, "double");
  shpbox.Z(2) = fread (fidp, 1, "double");
  shpbox.M(1) = fread (fidp, 1, "double");
  shpbox.M(2) = fread (fidp, 1, "double");
  
  ## FIXME: scan shp file in advance to assess nr of XY points, to be able to
  ##        preallocate rec array. May be difficult, .shp is not favorable for
  ##        it. Initial tries showed no significant speed advantage yet over
  ##        the incremental allocation scheme implemented in Ls. 611+ & 631+
  ##        for oct/plt style output. For ml-style it is much harder as
  ##        we'd need to preallocate a potentially very heterogeneous struct
  ##        array.

  ## Prepare for unsupported (in ML output) shape types
  unsupp = 0;
  ## Echo warning if dbug was set
  ign_mz = (outopts == 0) && dbug;
  ## Buffer for record data
  BUFSIZE = 10000;
  ## Read records, 1 by 1. Initialize final array
  vals = npt = npr = bbox = idx = nullsh = [];
  ## Init nr. of shapes read
  nshp = 0;
  ## Temp pointer to keep track of array size and increase it w. BUFSIZE rows
  ivals = 1;
  ## Provisionally assume file has M (measure) values
  has_M = true;

  ## Record index number (equals struct element number)
  ir = 1;
  do
    ## If the .shx file is there, skip directly to requested record
    if (have_shx)
      fseek (fidp, ridx(s_recs(ir), 1), "bof");
    endif
    if (dbug)
      printf ("Reading record %d ...\r", ir);
    endif
    val    = NaN(1, 6);
    ## Read record index number & record length
    val(5) = fread (fidp, 1, "int32", 0, "ieee-be");
    rlen   = fread (fidp, 1, "int32", 0, "ieee-be");

    ## Here we decide if the record need be read at all, based on s_recs
    if (isempty (s_recs) || ismember (double (val(5)), s_recs))
      ## => this record # has been desired; proceed. Read shape feature type
      val(6)  = fread (fidp, 1, "int32");
      rincl = 1;                                            ## see s_bbox, below

      ## Init values for this shape item
      tbbox  = NaN(1, 8);

      switch val(6)
        case 1                                              ## Point
          val(1:2) = fread (fidp, 2, "double");
          val(3:4) = NaN;
          tnpt     = 1;
          tnpr     = 0;

        case 11                                             ## PointZ
          val(1:4) = fread (fidp, 4, "double");
          tnpt     = 1;
          tnpr     = 0;

        case 21                                             ## PointM
          val(1:2) = fread (fidp, 2, "double");
          ## "Z"
          val(3)   = 0;
          ## M
          val(4)   = fread (fidp, 1, "double");
          tnpt     = 1;
          tnpr     = 0;

        case 8                                              ## Multipoint
          ## Read bounding box
          tbbox(1:4)       = fread (fidp, 4, "double");
          tnpt             = fread (fidp, 1, "int32");
          val(1:tnpt, 1:2) = reshape (fread (fidp, tnpt*2, "double"), 2, [])';
          val(1:tnpt, 3:4) = NaN;
          ## Copy rec index & type down
          val(2:tnpt, 5:6) = repmat (val(1, 5:6), tnpt-1, 1);
          tnpr             = 0;

        case 18                                             ## MultipointZ
          tbbox(1:4)       = fread (fidp, 4, "double");
          tnpt             = fread (fidp, 1, "int32");
          val(1:tnpt, 1:2) = reshape (fread (fidp, tnpt*2, "double"), 2, [])';
          ## Z min & max values
          tbbox(5:6)       = fread (fidp, 2, "double");
          ## Augment val array with Z values
          val(1:tnpt, 3)   = fread(fidp, tnpt, "double")';
          ## M min & max values
          tbbox(7:8)       = fread (fidp, 2, "double");
          if (val(1, 5) == 1)
            has_M = checkM (fidp, val(1, 6), shpbox);
          endif
          if (has_M)
            ## Augment val array with M values
            val(1:tnpt, 4) = fread(fidp, tnpt, "double")';
            ## Copy rec index & type down
            val(2:tnpt, 5:6) = repmat (val(1, 5:6), tnpt-1, 1);
            tnpr           = 0;
          endif

        case 28                                             ## MultipointM
          tbbox(1:4)       = fread (fidp, 4, "double");
          tnpt             = fread (fidp, 1, "int32");
          val(1:tnpt, 1:2) = reshape (fread (fidp, tnpt*2, "double"), 2, [])';
          ## Insert empty column for Z
          val(1:tnpt, 3:4) = NaN;
          if (val(1, 5) == 1)
            has_M = checkM (fidp, val(1, 6), shpbox);
          endif
          if (has_M)
            ## M min & max values
            tbbox(7:8)     = fread (fidp, 2, "double");
            ## Augment val array with M values
            val(1:tnpt, 4) = fread(fidp, tnpt, "double")';
          endif
          ## Copy rec index & type down
          val(2:tnpt, 5:6) = repmat (val(1, 5:6), tnpt-1, 1);
          tnpr             = 0;

        case {3, 5}                                         ## Polyline/-gon
          tbbox(1:4)       = fread (fidp, 4, "double");
          ## Read nparts, npoints, nparts pointers
          nparts           = fread (fidp, 1, "int32");
          tnpt             = fread (fidp, 1, "int32");
          tnpr             = fread (fidp, nparts, "int32")';
          ## Read XY point coordinates
          val(1:tnpt, 1:2) = reshape (fread (fidp, tnpt*2, "double"), 2, [])';
          ## No Z or M data
          val(1:tnpt, 3:4) = NaN;
          ## Copy rec index and type down
          val(2:tnpt, 5:6) = repmat (val(1, 5:6), tnpt-1, 1);

        case {13, 15}                                       ## Polyline/-gonZ
          tbbox(1:4)       = fread (fidp, 4, "double");
          ## Read nparts, npoints, nparts pointers
          nparts           = fread (fidp, 1, "int32");
          tnpt             = fread (fidp, 1, "int32");
          tnpr             = fread (fidp, nparts, "int32")';
          ## Read XY point coordinates
          val(1:tnpt, 1:2) = reshape (fread (fidp, tnpt*2, "double"), 2, [])';
          ## Z min & max values + data
          tbbox(5:6)       = fread (fidp, 2, "double");
          val(1:tnpt, 3)   = fread(fidp, tnpt, "double")';
          if (val(1, 5) == 1)
            has_M = checkM (fidp, val(1, 6), shpbox);
          endif
          if (has_M)
            ## M min & max values + data
            tbbox(7:8)       = fread (fidp, 2, "double");
            val(1:tnpt, 4)   = fread(fidp, tnpt, "double")';
            ## Copy rec index and type down
            val(2:tnpt, 5:6) = repmat (val(1, 5:6), tnpt-1, 1);
          endif

        case {23, 25}                                       ## Polyline/-gonM
          tbbox(1:4)       = fread (fidp, 4, "double");
          ## Read nparts, npoints, nparts pointers
          nparts           = fread (fidp, 1, "int32");
          tnpt             = fread (fidp, 1, "int32");
          tnpr             = fread (fidp, nparts, "int32")';
          ## Read XY point coordinates
          val(1:tnpt, 1:2) = reshape (fread (fidp, tnpt*2, "double"), 2, [])';
          ## No Z data
          val(1:tnpt, 3)   = NaN;
          if (val(1, 5) == 1)
            has_M = checkM (fidp, val(1, 6), shpbox);
          endif
          if (has_M)
            ## M min & max values + data
            tbbox(7:8)       = fread (fidp, 2, "double");
            val(1:tnpt, 4)   = fread(fidp, tnpt, "double")';
            ## Copy rec index and type down
            val(2:tnpt, 5:6) = repmat (val(1, 5:6), tnpt-1, 1);
          endif

        case 31                                             ## Multipatch
          tbbox(1:4)       = fread (fidp, 4, "double");
          ## Read nparts, npoints, nparts pointers, npart types
          nparts           = fread (fidp, 1, "int32");
          tnpt             = fread (fidp, 1, "int32");
          ## Npart types is just another row under npart pointers => read both.
          ## Provisionally transpose, this is later on reset after NaN insertion
          tnpr             = reshape (fread (fidp, nparts*2, "int32")', [], 2)';
          ## Read XY point coordinates
          val(1:tnpt, 1:2) = reshape (fread (fidp, tnpt*2, "double"), 2, [])';
          ## Z min & max values + data. Watch out for incomplete .shp file
          EOF = (ftell (fidp) > flngt-2);
          if (! EOF)
            tbbox(5:6)       = fread (fidp, 2, "double");
            val(1:tnpt, 3)   = fread(fidp, tnpt, "double")';
          endif
          fptr = ftell (fidp);
          EOF = (fptr > flngt-2);
          if (! EOF && has_M)
            has_M = checkM (fidp, val(1, 6), shpbox);
          endif
          ## M min & max values + data. Watch out for incomplete .shp file
          if (! EOF && has_M)
            tbbox(7:8)       = fread (fidp, 2, "double");
            val(1:tnpt, 4)   = fread(fidp, tnpt, "double")';
          endif
          ## Copy rec index and type down
          val(2:tnpt, 5:6) = repmat (val(1, 5:6), tnpt-1, 1);

        otherwise                                           ## E.g., null shape (0)
          ## Keep track of null shapes to avoid reading associated attributes
          if (val(1, 6) == 0)
            nullsh = [ nullsh ir ];
            rincl = 0;
          elseif (abs(val(1, 6)) > 31)
            error ("shaperead.m: unexpected shapetype value %f for feature # %d\n\
       Looks like a faulty shape file.", val(1, 6), ir);
          endif

      endswitch

      ## Check if (X, Y) are valid coordinates
      if (any (abs (val(:, 1:2)) > 1.797e308))
        ## Probably +/- Inf
        rincl = 0;
        if (dbug)
          printf ("shape# %d has no finite XY coordinates, skipped\n", ir);
        endif

      ## Detect if shape lies (partly) within or completely out of BoundingBox.
      ## Null shapes are automatically skipped
      elseif (! isempty (s_bbox))
        ## Just check if any shape feature bounding box corner lies in s_bbox
        tbox = [tbbox(1) tbbox(2); tbbox(1) tbbox(4); ...
                tbbox(3) tbbox(4); tbbox(3) tbbox(2); tbbox(1) tbbox(2)];
        rincl = 0;
        ## For polygons, polylines & multipatches:
        if (ismember (val(1, 6), [3, 13, 23, 5, 15, 25, 31])) ## Polygon/line
          ## To avoid undue CPU time for large shapes we take a shortcut:
          ## Check if shape feature lies in BoundingBox...
          [a, b] = inpolygon (tbox(:, 1), tbox(:, 2), sbox(:, 1), sbox(:, 2));
          ## ...or BoundingBox lies in polyline/gon. Faster than indiv. points
          [c, d] = inpolygon (sbox(:, 1), sbox(:, 2), tbox(:, 1), tbox(:, 2));
          if (any ([a; b; c; d]))
            ## At least one of the shape item bbox corners lies within s_bbox
            rincl = 1;
            bb_union = unique ([bb_union val(1, 5)]);
            ## FIXME still the case of polygon edges intersecting bounding box
            ##       w/o vertices inside bounding box. Clipping required for that
          endif
        elseif (ismember (val(1, 6), [1, 11, 21]))          ## (Multi-)Point
          ## Simply select all points within or on boundaries
          [a, b] = inpolygon (val(:, 1), val(:, 2), sbox(:, 1), sbox(:, 2));
          val = val(find(a), :);
          tnpt = size (val, 1);
          if (tnpt)
            rincl = 1;
          endif
        endif
        ## If clipping has been selected, clip all parts of the shape feature
        if (rincl && clip)
          ## What to do depends on shape type. Null and MultiPatch aren't clipped
          switch val(1, 6)
            case {3, 13, 23, 5, 15, 25, 31}                   ## Polyline/gon, Multipatch
              ## Temporarily silence Octave a bit, then call clippln
              warning ("off",  "Octave:broadcast", "local");
              [val, tnpt, tnpr] = clippln (val, tnpt, tnpr, sbox, val(1, 6));
            otherwise
              warning ("shaperead: unknown shape type found (%d) - ignored\n", ...
                       val(1, 6));
              rincl = 0;
              val = [];
          endswitch
          if (isempty (val))
            ## Don't include it and remove last added bb_union entry
            bb_union(end) = [];
            rincl = 0;
          else
            ## Update bounding box
            tbbox(1) = min (val(:, 1));
            tbbox(2) = min (val(:, 2));
            tbbox(3) = max (val(:, 1));
            tbbox(4) = max (val(:, 2));
          endif
        endif
      endif

      ## What to do with the val data, if to be included
      if (rincl)
        ## Keep track of nr of shape features read
        ++nshp;

        ## M-values < -1e39 really mean absent values
        im = find (val(:, 4) < -1e39);
        val(im, 4) = NaN;

        if (outopts < 3)
          ## Prepare an Octave or Matlab style struct optimized for fast
          ## plotting by inserting a NaN row after each polyline/-gon part
          nn = size (tnpr, 2);
          valt = NaN(tnpt + nn - 1, 6);
          ipt = 1;
          ttnpr = [tnpr(1, :) tnpt];
          dtnpr = diff (ttnpr);
          for ii=2:numel (ttnpr)
            valt(ipt:ipt+dtnpr(ii-1)-1, :) = val(ttnpr(ii-1)+1:ttnpr(ii), :);
            ipt += dtnpr(ii-1) + 1;
          endfor
          val = valt;
          tnpr(1, :) = (tnpr(1, :) + [0:numel(dtnpr)-1]);
        endif

        ## Shape either included by default or it lies in requested BoundingBox
        switch outopts
          case {0, 1}
            ## Return a ML compatible mapstruct. Attributes will be added later
            if (ign_mz && val(1, 6) >= 10)
              printf ("shaperead: M and Z values ignored for ml-style output\n");
              ign_mz = 0;
            endif
            switch val(1, 6)
              case {1, 11, 21}                              ## Point
                outs(end+1, 1).Geometry = "Point";
              case {8, 18, 28}                              ## Multipoint
                outs(end+1, 1).Geometry = "MultiPoint";
              case {3, 13, 23}                              ## Polyline
                outs(end+1, 1).Geometry = "Line";
              case {5, 15, 25}                              ## Polygon
                outs(end+1, 1).Geometry = "Polygon";
              otherwise
                if (outopts == 1)
                  ## "Extended" ML-style output struct
                  if (val(1, 6) == 31)                      ## MultiPatch
                    outs(end+1, 1).Geometry = "MultiPatch";
                    outs(end, 1).Parts = tnpr;
                  endif
                else
                  if (! unsupp)
                    warning ("shaperead: shapefile contains unsupported shape \
types\n");
                    outs = oatt = [];
                    return
                  endif
                  outs(end+1, 1).Geometry = val(1, 6);
                endif
            endswitch
            ## Omit BoundingBox for Point
            if (all ([1, 11, 21] - val(1, 6)))
              outs(end).BoundingBox = reshape (tbbox(1:4), 2, [])';
            endif
            if (s_geo)
              outs(end, 1).Lon = val(:, 1)';
              outs(end, 1).Lat = val(:, 2)';
            else
              outs(end, 1).X = val(:, 1)';
              outs(end, 1).Y = val(:, 2)';
            endif
            ## (ML-incompatible) add Z- and optional M-values, if any
            if (outopts == 1 && any (isfinite (val(:, 4))))
              outs(end, 1).M = val(:, 4)';
              ## FIXME Decision needed if field Geometry should reflect the type
              ##  outs{end}.Geometry = [ outs{end}.Geometry "M"];
            endif
            if (outopts == 1 && any (isfinite (val(:, 3))))
              outs(end, 1).Z = val(:, 3)';
              ## FIXME Decision needed if field Geometry should reflect the type
              ## outs{end}.Geometry = [ outs{end}.Geometry "Z"];
            endif
            ## (ML-incompatible) add Z- and optional M-values, if any
            if (dbug)
              ## Add a field with shape feature identifier for boundingbox
              outs(end, 1).___Shape_feature_nr___ = val(1, 5);
            endif

          case {2}
            ## Return an Octave style struct.
            ## Append to vals array; keep track of appended nr. of rows
            lvals = size (vals, 1);
            lval = size (val, 1);
            if ((size (vals, 1) - ivals) < lval)
              ## Increase size of vals
              vals = [ vals; NaN(BUFSIZE, 6) ];
            endif
            vals(ivals:ivals+lval-1, :) = val;
            idx = [idx ; ivals];
            ivals += lval;
            ## Add a row of NaNs
            vals(ivals, :) = NaN(1,6);
            ++ivals;
            ## Append the other arrays
            npt = [npt; tnpt];
            if (isempty (npr))
              npr = {tnpr};
            else
              npr = [npr; {tnpr}];
            endif
            bbox = [bbox; tbbox];

          case {3}
            ## Return a compressed Octave style struct.
            ## Simply append to vals array; keep track of appended nr. of rows
            lvals = size (vals, 1);
            lval = size (val, 1);
            if ((size (vals, 1) - ivals) < lval)
              ## Increase size of vals
              vals = [ vals; NaN(BUFSIZE, 6) ];
            endif
            vals(ivals:ivals+lval-1, :) = val;
            idx = [idx ; ivals];
            ivals += lval;
            ## Append the other arrays
            npt = [npt; tnpt];
            if (isempty (npr))
              npr = {tnpr};
            else
              npr = [npr; {tnpr}];
            endif
            bbox = [bbox; tbbox];

          otherwise
        endswitch
      endif

    else
      ## Record not in s_recs list, skip rest of record
      fread (fidp, (rlen*2), "char=>char");

    endif

    ## Next .shp record
    ++ir;
    
  until (ftell (fidp) > flngt - 3 ||
         (have_shx && ir > numel (s_recs)));                ## i.e., within last
                                                            ## file bytes or all
                                                            ## req. records read
  fclose (fidp);

  ## If no shape was read, or none fitted within BoundingBox, return empty
  if (nshp < 1)
    outs = oatt = {};
    return;
  endif
  
  ## ------------------------ Post-processing ------------------------
  if (outopts > 1)
    ## Octave-style outstruct. Truncate vals to proper length
    vals = vals(1:ivals-1, :);
    if (outopts == 2 && ! isempty (vals))
      ## Remove last NaN row
      vals(end, :) = [];
    endif
    if (isempty (vals))
      s_recs = 0;
    endif
  endif

  if (dbug)
    if (outopts <= 1)
      ii = numel (outs);
    else
      ii = numel (npt);
    endif
    printf ("\n%d records read.            \n", ii);
  endif

  ## For Octave style output, add separate arrays to output struct
  if (outopts >= 2)
    outs.shpbox = shpbox;
    outs.vals   = vals;
    outs.bbox   = bbox;
    outs.npt    = npt;
    outs.npr    = npr;
    outs.idx    = idx;
    ## Clear memory (for very large shape files)
    clear shpbox vals bbox npt npr idx val;
    outs = setfield (outs, "fields", {"shpbox", "vals", "bbox", "npt", "npr"});
  endif

  ## Clean up s_recs and bb_union
  if (! isempty (s_bbox))
    s_recs = sort (unique (bb_union));
  else
    s_recs = sort (unique (s_recs));
  endif

  ## Weed out any null shape records to prevent reading their attributes
  if (! isempty (nullsh))
    if (! isempty (s_recs))
      s_recs = [1:numel (npt)];
    endif
    s_recs (ismember (nullsh, s_recs)) = [];
  endif

  ## ---------------------- 3. .dbf ----------------------

  if (iscell (s_atts) && isempty (s_atts))
    ## {} indicates no attributes to be read. Morph into ""
    return;
  endif
  ## Check if dbfread is available
  if (isempty (which ("dbfread")))
    printf ("shaperead.m: dbfread function not found. No attributes will be added.\n");
    printf ("             (io package installed and loaded?)\n");
    oatt = {};

  else
    ## Try to read the .dbf file
    try
      atts = dbfread ([ bname ".dbf" ], s_recs, s_atts);
      if (outopts < 2)
        ## First check if any attributes match fieldnames; if so, pre-/append "_"
        tags = {"Geometry", "BoundingBox", "X", "Y", "Lat", "Lon"};
        for ii=1:numel (tags)
          idx = find (strcmp (tags{ii}, atts(1, :)));
          if (! isempty (idx))
            atts(1, idx) = ["_" atts{1, idx} "_"];
          endif
        endfor
        ## Matlab style map-/geostruct. Divide attribute values over struct elems
        if (nargout < 2)
          ## Attributes appended to outs struct
          for ii=1:size (atts, 2)
            [outs.(atts{1, ii})] = deal (atts(2:end, ii){:});
          endfor
          oatt = [];
        else
          ## Attributes separately in oatt struct
          oatt(size (atts, 1) - 1).(atts{1, 1}) = [];
          for ii=1:size (atts, 2)
            [oatt.(atts{1, ii})] = deal (atts(2:end, ii){:});
          endfor
          oatt = oatt';
        endif
      else
        ## Octave output struct. Add attributes columns as struct fields
        ## Check if any attributes match fieldnames; if so, pre-/append "_"
        tags = {"shpbox", "vals", "bbox", "npt", "npr", "idx", "Geometry", ... 
                "BoundingBox", "X", "Y", "Lat", "Lon"};
        for ii=1:size (atts, 2)
          outs.fields(end+1) = atts{1, ii};
          if (islogical (atts{2, ii}) || isnumeric (atts{2, ii}))
            outs = setfield (outs, atts{1, ii}, cell2mat (atts(2:end, ii)));
          else
           outs = setfield (outs, atts{1, ii}, atts(2:end, ii));
         endif
        endfor
        atts = [];
      endif
    catch
      printf ("shaperead: file %s couldn't be read;\n", [ bname ".dbf" ]);
      printf ("           no attributes appended\n");
    end_try_catch
  endif

endfunction


function has_M = checkM (fidp, ft, shpbox)

  has_M = true;
  ## Check if we do have M-values. If so, the next 3 4-byte words 
  ## comprise the next record nr, rec length and shape type
  fptr   = ftell (fidp);
  tmp    = fread (fidp, 2, "int32", 0, "ieee-be");
  tmp(3) = fread (fidp, 1, "int32");
  if (tmp(1) == 2 && tmp(3) == ft)
    ## Looks like next record header + next shape feature type => no M
    has_M = false;
    shpbox.M(1) = 0;
    shpbox.M(2) = 0;
  endif
  fseek (fidp, fptr, "bof");

endfunction


## Only check input validation. I/O is tested in shapewrite.m
%!error <arg. #2 char value should be one of> shaperead ('tst.shp', 'j');
%!error <shaperead: arg. #2 integer value out of range> shaperead ('tst.shp', 7);
%!error <illegal value for arg. #2:> shaperead ('tst.shp', 'deb')
%!error < property name expected> shaperead ('tst.shp', "ml", []);
%!error <numeric or logical value expected> shaperead ('tst.shp', 'deb', {});

