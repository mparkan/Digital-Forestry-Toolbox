## Copyright (C) 2015-2020 Philip Nienhuis
## Copyright (C) 2018-2020 Matthew Parkan
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
## @deftypefn  {Function File} [@var{status}] = dbfwrite (@var{fname}, @var{data})
## Write data in a cell array to a dbf (xBase) file, provisionally dBase III+.
##
## @var{fname} must be a valid file name, optionally with '.dbf' suffix.
## @var{data} should be a cell array of which the top row contains column
## names (character strings, each max. 10 characters; longer column names will
## be truncated).  Each column must contain only one class of data, except of
## course the top entry (the column header).  Integers interspersed in double
## type colums will be written as doubles.  Data types that can be written are
## character (text string), numeric (integer and float, the latter with 6
## decimal places), and logical.
##
## Output argument @var{status} is 1 if the file was written successfully, -1 if
## one or more data columns were skipped, 0 otherwise.  If 0 the incomplete file
## will be deleted as well.
##
## Provisionally only dBase v. III+ files without memos can be written.
##
## @seealso{dbfread}
## @end deftypefn

## Authors: Philip Nienhuis <prnienhuis@users.sf.net>,whos
##          Matthew Parkan <matthew.parkan@gmail.com>
## Created: 2014-12-24

function [status] = dbfwrite (fname, data)

  status = 0;
  ## Input validation
  if (! ischar (fname))
    error ("dbfwrite: file name expected for argument #1\n");
  elseif (! iscell (data))
    error ("dbfwrite: cell array expected for argument #2\n");
  elseif (! iscellstr (data (1, :)))
    error ("dbfwrite: column header titles (text) expected on first row of data\n");
  endif
  ## Column headers length cannot exceed 10 characters
  toolong = [];
  for ii=1:size (data, 2)
    title = data{1, ii};
    if (length (title) > 10)
      toolong = [ toolong, ii ];
      data(1, ii) = title(1:10);
    endif
  endfor
  if (! isempty (toolong))
    ## Truncate headers if required and check for uniqueness
    warning ("dbfwrite: one or more column header(s) > 10 characters - truncated\n");
    fmt = [repmat("%d ", 1, numel (toolong))];
    printf ("Applies to column(s): %s\n", sprintf (fmt, toolong));
    if (numel (unique (data(1, :))) < numel (data(1, :)))
      error ("dbfwrite: column headers aren't unique - please fix data\n");
    endif
  endif

  ## Assess nr of records
  ## Data contains header row. Data toprow = 2
  nrecs = size (data, 1) - 1;
  tr = 2;

  ## Check file name
  [pth, fnm, ext] = fileparts (fname);
  if (isempty (ext))
    fname = [fname ".dbf"];
  elseif (! strcmpi (ext, ".dbf"))
    error ("dbfwrite: file name should have a '.dbf' suffix\n");
  endif
  ## Try to open file
  fid = fopen (fname, "w", "ieee-le");
  if (fid < 0)
    error ("dbfwrite: could not open file %s\n", fname);
  endif

  unwind_protect
    ## Start writing header
    ## Provisionally assume dbase III+ w/o memos
    fwrite (fid, 3, "uint8");

    ## Date of last update (YYMMDD), with YY the number of years since 1900
    t = now;
    upd = datevec(t) - [1900, 0, 0, 0, 0, 0];
    fwrite (fid, uint8(upd(1:3)), "uint8");

    ## Number of records in the table
    fwrite (fid, nrecs, "uint32");
    ## The next two uint16 fields are to be written later, just fill temporarily
    pos_lhdr = ftell(fid);

    fwrite (fid, 0, "uint32");

    ## Another place holder, write enough to allow next fseek to succeed
    fwrite (fid, uint32 (zeros (1, 7)), "uint32");

    ## Write record descriptors
    nfields  = size (data, 2);
    fldtyp   = "";
    fldlngs  = {};
    reclen   = 1;                                             ## "Erased" byte first
    fseek (fid, 32, "bof");

    RR = zeros (32, nfields, "uint8");
    colskipped = 0;
    for ii=1:nfields
      decpl = 0;
      recdesc = sprintf ("%d", uint32 (zeros (1, 8)));
      recdesc(1:10) = strjust (sprintf ("%10s", data{1, ii}), "left"); ## Field name
      ## Be strict on mixing char and numeric; this implies logicals
      ## interspersed in numeric type column won't be accepted either
      if (all (cellfun (@isnumeric, data(tr:end, ii), "uni", 1)))
        ## We're lax on interspersed integers, they'll be cast to double
        if (isinteger ([data{tr:end, ii}]) ||
            all ([data{tr:end, ii}] - floor([data{tr:end, ii}]) < eps))
          ftype = "N";
          decpl = 0;
        else
          ftype = "F";
          ## ML compatibility for .dbf/.shp file: 6 decimal places
          decpl = 6;
        endif
        fldlng = 20;
      elseif (all (cellfun (@ischar, data(tr:end, ii), "uni", 1)))
        ftype = "C";
        fldlng = max (cellfun (@(x) length(x), data(tr:end, ii)));
      elseif (all (cellfun (@islogical, (data(tr:end, ii)), "uni", 1)))
        ftype = "L";
        fldlng = 1;
      else
        warning ("dbfwrite: heterogeneous data types in column %d ('%s'), \
skipped.\n", ii, data{1, ii});
        RR(:, end) = [];
        nfields--;
        colskipped = 1;
        continue ;
        ## unwind_protect_cleanup takes care of closing & wiping file
      endif
      recdesc(12) = ftype;                                    ## Field type
      fldtyp      = [ fldtyp ftype ];
      recdesc(17) = uint8 (fldlng);                           ## Field length
      recdesc(18) = uint8 (decpl);                            ## Decimal places
      recdesc(32) = "\0";                                     ## Fill to byte# 32

      RR(:, ii) = recdesc';
      reclen += fldlng;
      fldlngs = [ fldlngs; sprintf("%d", fldlng) ];
    endfor

    fwrite (fid, RR, "char");

    ## Write header record terminator
    fwrite (fid, 13, "uint8");
    ## Remember position
    fpos_data = ftell (fid);
    ## Write missing data in header
    fseek (fid, pos_lhdr, "bof");
    fwrite (fid, fpos_data, "uint16");
    fwrite (fid, reclen, "uint16");

    ## Write data2
    fseek (fid, fpos_data, "bof");

    ## Determine data record format. "Erased byte" is first char of format
    fmt = "\0";
    for j = 1:nfields
      switch fldtyp(j)
          case "C"                                            ## character
            txt = ["%", fldlngs{j}, "s"];
          case "N"                                            ## numeric
            txt = ["%" fldlngs{j} "d"];
          case "L"                                            ## logical
            txt = ["%", fldlngs{j}, "c"];
          case "F"                                            ## float
            txt = ["%" fldlngs{j} "f"];
          case "D"                                            ## date; currently inactive
            ## txt = sprintf (["%" fldlngs{jj} "s"], data{ii, jj});
          otherwise
        end
        fmt = [fmt, txt];                                     ## append format
    end

    ## Convert boolean attributes to Y/N characters
    str_logical = {"N", "Y"};
    for jj = find (fldtyp == "L")
      data(2:end, jj) = str_logical (double ([data{2:end, jj}] + 1))';
    end

    ## Write data in ~100 MB chunks to avoid overflow. First find an optimal
    ## chunk size as max. nr. of records in a chunk= (= nr.of rows in data)
    chunk_sz = floor (1e8 / reclen);
    for ii=1 : chunk_sz : nrecs
      ## Reshape chunk of data matrix
      T = [data(ii+1:min (ii+chunk_sz, nrecs+1), :)'(:)];
      blob = sprintf (fmt, T{:}); 
      ## Write blob to file
      fwrite (fid, blob, "char");
    endfor
    status = 1;

  unwind_protect_cleanup
    fclose (fid);
    if (! status)
      printf ("dbfwrite: removing incomplete file %s.\n", fname);
      unlink (fname);
    elseif (colskipped)
      status = -1;
    endif
  end_unwind_protect

endfunction
