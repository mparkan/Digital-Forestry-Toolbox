## Copyright (C) 2014-2020 Philip Nienhuis
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
## @deftypefn  {Function File} [@var{data}, @var{datinfo}] = dbfread (@var{fname})
## @deftypefnx {Function File} [@var{data}, @var{datinfo}] = dbfread (@var{fname}, @var{recs})
## @deftypefnx {Function File} [@var{data}, @var{datinfo}] = dbfread (@var{fname}, @var{recs}, @var{cols})
## @deftypefnx {Function File} [@var{data}, @var{datinfo}] = dbfread (@var{fname}, @var{recs}, @var{cols}, @var{re})
## Read contents of a dbase (dbf) file, provisionally dbase III+, IV or V.
##
## @itemize
## @item
## @var{fname} should be the name of a valid dbase file; the file extension
## isn't required.
##
## @item
## @var{recs} can be an integer or logical array containing record numbers or
## record indicators for those records that need to be returned.  If omitted,
## all records are read.  Indices supplied in @var{recs} can be specified in
## any order, but the returned data are sorted in order of records in the file.
##
## @item
## @var{cols} can be a logical, integer, cellstr or character array indicating
## from which file columns the data should be returned.  If a numeric array is
## supplied, it is considered to be like a logical array if the maximum entry
## value equals 1.  Character arrays should have column names stacked in the
## vertical (first) dimension.  @var{cols} cellstr or char arrays can be
## supplied in any order, yet the returned data column order matches that of
## the columns order in the dbase file.  For dbase files containing multiple
## columns with the same name, specify a numeric or logical array to select
## columns to be returned.  If omitted, data from all file columns are
## returned.
##
## @item
## If a value of 1 or true is entered for @var{re}, dbfread also tries to
## return data from erased records.  No guarantee can be given for these data
## to be correct or consistent!  If omitted, erased records are skipped.
##
## @item
## Return value @var{data} is a N+1 x M cellstr array where the uppermost row
## contains the column names and the rest of the rows the record data.
##
## @item
## Optional return argument @var{datinfo} is a struct array containing various
## information of the dbase file and record build-up. 
## @end itemize
##
## Arguments @var{recs} and @var{cols} need not be as long as the number of
## records and columns in the file, resp.; dbfread will stop reading data if
## any of @var{recs} or @var{cols} (if supplied) is exhausted.
##
## Sometimes dbase files contain records indicated as being erased.  The data
## in such records is silently skipped, unless the @var{re} flag is set
## and/or @var{recs} is supplied and erased records happen to be present in the
## requested record numbers.
##
## Examples:
##
## @example
##   A = dbfread ("file.dbf");
##   (returns all data in file.dbf in array A) 
## @end example
##
##
## @example
##   [A, B] = dbfread ("file.dbf", [], ["colB"; "colF"]);
##   (returns all data in columns named "colB" and "colF" from
##    file.dbf in array A and information on the database
##    build-up in struct B) 
## @end example
##
## @example
##   A = dbfread ("file.dbf", [0 1 0 0 1 0 0]);
##        -or-
##   A = dbfread ("file.dbf", [2 5]);
##   (returns data from record numbers 2 and 5 in
##    file.dbf in array A)
## @end example
##
## @example
##   A = dbfread ("file", [0 1 0 0 1 0]);
##   (returns data from record numbers 2 and 5 in
##    file.dbf in array A)
## @end example
##
## @example
##   [~, B] = dbfread ("file.dbf", 0);
##   (to returns info on column names and number of
##    records, plus more info)
## @end example
##
## @example
##   [A] = dbfread ("file", [], @{"Header1", "Col5"@});
##   (returns data from columns with names (headers) 
##    Header1 and Col5, resp.)
## @end example
##
## @seealso{xlsread}
## @end deftypefn

## Author: Philip Nienhuis <prnienhuis@users.sf.net>
## Created: 2014-11-03

## References:
## http://ulisse.elettra.trieste.it/services/doc/dbase/DBFstruct.htm
## http://www.dbf2002.com/dbf-file-format.html
## http://www.dbase.com/KnowledgeBase/int/db7_file_fmt.htm

function [data, datinfo] = dbfread (fname, recs=[], cols=[], rd_erased="")

  ## Check file name
  if (ischar (fname))
    [~, fnm, ext] = fileparts (fname);
    if (isempty (ext))
      fname = [fname ".dbf"];
    endif
  else
    error ("dbfread: file name expected for arg # 1.\n");
  endif

  ## Check recs arg. If needed turn into indices
  if (! isempty (recs))
    if (! (isnumeric (recs) || islogical (recs)))
      error ("dbfread: numeric or logical array expected for arg # 2\n");
    elseif (isnumeric (recs))
      if (any (recs < 0))
        error ("dbfread: illegal record selection indices\n");
      elseif (min (recs) == 0 && max (recs) < 2)
        recs = find (recs);
        if (isempty (recs))
          recs = -1;
        endif
      endif
    elseif (islogical (recs))
      recs = find (recs);
    endif
    recs = sort (recs);
    endif

  ## Check cols arg. If needed turn into indices
  if (! isempty (cols))
    if (! (isnumeric (cols) || ischar (cols) || iscellstr (cols) || islogical (cols)))
      error ("dbfread: numeric, cellstr, logical, or character array expected for arg # 3\n");
    elseif (isnumeric (cols))
      if (any (cols < 0))
        error ("dbfread: illegal column selection indices\n");
      elseif (min (cols) == 0 && max (cols) < 2)
        cols = find (cols);
        if (isempty (cols))
          cols = -1;
        endif
      endif
    elseif (islogical (cols))
      cols = find (cols);
    endif
  endif

  ## Check rd_erased arg.
  if (! isempty (rd_erased) && ! (islogical (rd_erased) || isnumeric (rd_erased)))
    error ("dbfread: numeric or logical value expected for arg # 4\n");
  endif

  ## Open file
  fid = fopen (fname, "r");
  if (fid < 0)
    error ("dbfread: file %s couldn't be opened.\n", fname);
  endif
  ## Rewind, just to be sure
  fseek (fid, 0, "bof");

  ## First check proper type
  fbyte = uint8 (fread (fid, 1, "uint8"));
  ## Provisional type check.
  ##   3 = dbase III+ w/o memos
  ##  83 = dbase III+ w memos
  ## ... <to be cont'd>
  vsn = bitand (fbyte, 7);
  if (! ismember (vsn, [3, 4, 5]))
    error ("dbfread: unsupported file type, only dbase III[+], IV & V supported.\n");
  endif
  ## Memos present for fbyte == 83, and bits 3 and/or 7 set (1-based bit pos.)
  hasmemo = (fbyte == 83) || (bitand (fbyte, 8)) || (bitand (fbyte, uint8 (128)));

  ## Start reading header
  lasty = fread (fid, 1, "uint8") + 1900;                   ## Last dbf update
  lastm = fread (fid, 1, "uint8");                          ## month
  lastd = fread (fid, 1, "uint8");                          ## day

  nrecs = fread (fid, 1, "uint32");
  lhdr  = fread (fid, 1, "uint16");
  recl  = fread (fid, 1, "uint16");

  ## Field descriptors
  nfields  = 0;
  fseek (fid, 32, "bof");
  fdesc = fread (fid, 32, "char=>char")';
  do
    ++nfields;
    ## Get fields into struct
    dbf(nfields).fldnam = deblank (fdesc(1:11));            ## Field name
    dbf(nfields).fldtyp = fdesc(12);                        ## Field type
                                                            ## Skip field dspl.
    dbf(nfields).fldlng = int32 (fdesc(17));                ## Field length
    dbf(nfields).flddec = int32 (fdesc(18));                ## Decimal places
    dbf(nfields).fldflg = int8 (fdesc(19));                 ## Flags
    ## Get next field descriptors
    fdesc = fread (fid, 32, "char=>char")';
  until (ftell (fid) >= lhdr)
  ## Seek to position after header terminator byte
  fseek (fid, lhdr, "bof");

  ## Read rest of data. Skip if no records need be read. Turn into char array
  if (isempty (recs))
    txt = fread (fid, [recl, nrecs], "char=>char")';
  else
    txt = fread (fid, [recl, recs(end)], "char=>char")';
  endif
  ## .dbf file is no longer needed
  fclose (fid);

  ## Preallocate upper data row
  data = cell (1, numel (dbf));
  data(1, :) = {dbf.fldnam};

  ## If required, select requested records. Beware; truncate indices > nrecs
  if (any (recs > nrecs))
    recs (find (recs > nrecs)) = [];
    ## Check if we still have a selection...
    if (isempty (recs))
      ## No more, signal this to below code
      recs = -1;
    endif
  endif
  if (! isempty (recs) && ! isempty (txt))
    if (any (recs < 0))
      ## No data returned. Explore erased records anyway
      wipedrec = txt(:, 1)' == "*";
      txt = "";
      recs = 0;
      scol = numel (dbf);
    else
      txt = txt(recs, :);
      ## Preallocate data cell array for selected records
      data = [data; cell(numel (recs), numel (dbf))];
    endif
  else
    ## Preallocate data cell array for all records
    data = [data; cell(nrecs, numel (dbf))];
  endif

  if (! isempty (txt))
    ## There's something to read ;-)  First, read memo file, if any
    if (hasmemo)
      ## Also open accompanying .dbt file
      fjd = fopen ([fnm ".dbt"], "r");
      if (fjd < 0)
        warning ("dbfread: associated memo file (%s) couldn't be opened.\n", ...
                  [fnm ".dbt"]);
        printf ("(dbfread: skipping memo fields)\n");
        memos = {};
      else
        fseek (fjd,0, "bof");
        ## Read memo fields
        memos = fread (fjd, Inf, "char=>char")';
        fclose (fjd);
        ## Pimp memos: replace all non-alphanumeric chars by space
        memos(find (int8(memos) > 122)) = char(32);
        memos(find (int8(memos) <  32)) = char(32);
        switch fbyte
          case {83, 131}                                        ## Dbase III[+]
            ## Make it into a Nx512 char array. Pad unit length = multiple of 512
            pad = ceil (length (memos) / 512) * 512 - length (memos);
            memos = [memos repmat(int8 (32), 1, pad)];
            memos = cellstr (reshape (memos, 512, [])');
          otherwise
            ## FIXME: Dbase V, VII to follow
        endswitch
      endif
    endif

    txtp = 2;
    ## Init output array column pointer
    scol = 0;
    for ii=1:numel (dbf)
      ## First process selection if arg. cols was given...
      if (! isempty (cols))
        try
          ## try-catch, as cols array < nr. of cols is allowed.
          ## Switch dependent of cols input arg. type
          switch class (cols)
            case "cell"
              getcol = any (strcmp (cols{scol+1}, dbf(ii).fldnam));
            case "char"
              getcol = strcmp (strtrim (cols(scol+1, :)), dbf(ii).fldnam);
            case "double"
              getcol = ! isempty (find (cols == ii));
            otherwise
          endswitch
        catch
          getcol = 0;
        end_try_catch
      else
        ## No cols arg. was given, so we'll read this column anyway
        getcol = 1;
      endif

      if (getcol)
        ## Read column # ii
        ++scol;
        fld = txt(:, txtp : txtp+dbf(ii).fldlng - 1);
        switch (dbf(ii).fldtyp)
          case {"B", "G"}
            ## Block number into .dbt file, other than memo
            data(2:end, scol) = num2cell (str2double (fld));
          case "C"
            ## Text
            data(2:end, scol) = cellstr (fld);
          case "D"
            ## Date
            dlf = cellstr (fld);
            ## Catch empty date fields. Put ridiculous value in
            dlf(find (strcmp (cellstr (fld), "00000000"))) = "99991231";
            dlf(find (cellfun ("isempty", dlf))) = "99991231";
            dlf = datenum (dlf, "yyyymmdd");
            ## Reset temp values for empty dates
            dlf (dlf >= 3652425) = 0;
            data(2:end, scol) = num2cell (dlf);
          case "L"
            ## Logical / boolean
            data(2:end, scol) = false;
            data(regexpi (fld, "[yt]")+1, scol) = true;
          case {"F", "N"}
            ## Numeric
            data(2:end, scol) = num2cell (str2double (fld));
          case "M"
            ## Memo field pointer into .dbt file
            if (! isempty (memos))
              switch (fbyte)
                case {83, 131}                              ## Dbase III[+]
                  idx = str2double (fld);
                  idx (find (isnan (idx))) = 1;
                  data(2:end, scol) = memos(idx);
                otherwise
                  ## FIXME other .dbf file versions to be implemented...
              endswitch
            endif
          otherwise
        endswitch
      else
        ## Remove data column header & associated column from output array
        data(:, scol+1) = [];
      endif

      ## Next column of data from .dbf file
      txtp += dbf(ii).fldlng;

    endfor

    ## Only now check for erased records, to avoid false positives below
    wipedrec = txt(:, 1)' == "*";
    if (! rd_erased)
      ## No erased records in user-supplied record selection => no worries
      data (find (wiperec), :) = [];
    elseif (sum (wipedrec) && ! isempty (recs))
      ## User did request erased record # => warn
      warning ("(dbfread: %d erased records read\n", wipedrec);
    endif

  else
    wipedrec = [];
    scol = 0;

  endif

  if (isempty (data))
    data = {};
  endif

  if (nargout <= 1)
    ## Only data requested
    datinfo = [];
  else
    ## Also (or only)  dbf info struct requested
    if (isempty (recs))
      ## Infer nr. of records from file
      recs = size (data, 1) - 1;
    endif
    datinfo.type = fbyte;
    datinfo.date = datenum (lasty, lastm, lastd);
    datinfo.erasedrec = find (wipedrec);
    datinfo.nrec = nrecs;
    datinfo.srecs = recs;
    datinfo.recl = recl;
    datinfo.ncols = numel (dbf);
    datinfo.scols = scol;
    datinfo.data = dbf;
  endif

endfunction
