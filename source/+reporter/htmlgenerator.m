%HTMLGENERATOR Class for generating HTML files

classdef htmlgenerator
    
  properties(SetAccess = 'private', Hidden)
    filename = '';
    pathname = '';
    fid = [];
    title = '';
    pagewidth = 0;
  end
  
  methods
  
    %-------------
    function hg = htmlgenerator(title, pathname, filename, pagewidth)
    %-------------
    %Constructor
      if nargin < 4
        pagewidth = 1024;
        if nargin < 3
          filename = [title '.htm'];
          if nargin < 2
            pathname = pwd;
          end
        end
      end
      hg.title = title;
      hg.filename = filename;
      hg.pathname = pathname;
      hg.pagewidth = pagewidth;
    end
    
    %-------------
    function hg = start(hg)
    %-------------
    %START Method to open the HTML file and write headers
      hg.fid = fopen(fullfile(hg.pathname,hg.filename), 'w', 'n', 'UTF-8');
      write(hg, html_header(hg));
    end
    
    %-------------
    function stri = hline(hg)
    %-------------
    %HLINE Method to write a horizontal line
      stri = sprintf('<br>\n<hr>\n<br>\n');
      %These methods return a string if asked for output, otherwise write
      %to file

      if nargout < 1
        write(hg, stri);
      end
    end
    
    %-------------
    function stri = table(hg, content, boldcells, width)
    %-------------
    %TABLE Method to insert a table with content specified by a cell.
    
      if nargin < 4
        width = round(0.625*hg.pagewidth);
        if nargin < 3
          boldcells = zeros(size(content));
          if nargin < 2
            myfailed('Too few input arguments.');
          end
        end
      end
      
      stri = sprintf('<table frame=box rules=all width=%d>\n', width);
      for i = 1:size(content, 1)
        stri = [stri sprintf('<tr>\n')]; %#ok<AGROW>
        for j = 1:size(content, 2)
          if boldcells(i, j) == 2
            b1 = '<font color="red">';
            b2 = '</font>';
          elseif boldcells(i, j)
            b1 = '<b>';
            b2 = '</b>';
          else
            b1 = '';
            b2 = '';
          end
          stri = [stri sprintf('<td>%s%s%s</td>\n', b1, content{i,j}, b2)]; %#ok<AGROW>
        end
        stri = [stri sprintf('</tr>\n')]; %#ok<AGROW>
      end
      stri = [stri sprintf('</table>\n<br clear=left>\n')];
      
      if nargout < 1
        write(hg, stri);
      end
      
    end
    
    %-------------
    function stri = box(hg, text, width)
    %-------------
    %BOX Method to insert a text box
    if nargin < 3
      width = hg.pagewidth;
      if nargin < 2
        myfailed('Too few input arguments.');
      end
    end
    stri = sprintf(['<table frame=box width=%d height=200>\n'...
      '<td valign = "top">\n%s\n</td>\n</table>\n'], width, text);
    if nargout < 1
      write(hg, stri);
    end
    end
    
    %-------------
    function stri = headline(hg, headtext, prio)
    %-------------
    %HEADLINE Method to write a headline in bold
    if nargin < 3
      prio = 3;
    end
    stri = sprintf('<h%d>%s</h%d>\n', prio, headtext, prio);
    if nargout < 1
      write(hg, stri);
    end
    end
    
    %-------------
    function stri = text(hg, paragraph, alignment)
    %-------------
    %TEXT Method to write a paragraph of text
    if nargin == 3
      alstr = sprintf(' align=%s',alignment);
    else
      alstr = '';
    end
      stri = sprintf('<p%s>\n%s\n</p>\n',alstr,paragraph);
      if nargout < 1
        write(hg, stri);
      end
    end
    
    %------------
    function gsstri = ftext(hg)
    %------------
    %FTEXT Method to write a paragraph of text with formatted headlines.
    global SET
    stri = [];
    gsstri = []; %Gensvar string (unformatted)
    cr = [];
    prow = 1; %Variable to keep track of row number in current paragraph
    %Convert from Matlab string matrix to viewable text
    comment = SET(1).Report.Comments;
    newlines = [0 regexp(comment,'(\n)') numel(comment)+1];
    for row = 2:numel(newlines)
      crprev = cr;
      if ~isempty(crprev)
        crprev = regexprep(crprev,{'\*','\.'},{'\\*','\\.'});
      end
      cr = comment(newlines(row-1)+1:newlines(row)-1);
      gsstri = [gsstri cr sprintf('\n')]; %#ok<AGROW>
      if isempty(cr) %Insert new paragraph if a line is blank
        stri = [stri sprintf('\n</p>\n<p>\n')]; %#ok<AGROW>
        prow = 1;
      else
        if prow == 2 %Create headline according to rules of formatting (see manual)
          if strcmp(upper(crprev),crprev)
            stri = regexprep(stri,crprev,sprintf('<h3>%s</h3>',crprev)); %superheadline
            prow = 1;
          else
            stri = regexprep(stri,crprev,sprintf('<h4>%s</h4>',crprev)); %subheadline
          end
        elseif prow > 2
          stri = [stri sprintf('<br>\n')]; %#ok<AGROW>
        end
        stri = [stri sprintf('%s\n',cr)]; %#ok<AGROW>
        prow = prow + 1;
      end
    end
    hg.text(stri);
    end
    
    %-------------
    function stri = image(hg, imgname, imgsource)
    %-------------
    %IMAGE Method to add an image. Also stores it to disk in
    % the same folder as the HTML file. imgsource can be a 
    % location or an image matrix.
    
      if nargin < 3
        stri = sprintf('<img src="%s">\n<br clear=left>\n', imgname);
      else
        imgdest = [hg.pathname filesep imgname];
          if ischar(imgsource)
            try
              imgsource = imread(imgsource);
            catch me
              mydispexception(me);
              error('Could not read image file');
            end
          end
        try
          imwrite(imgsource, imgdest);
          stri = sprintf('<img src="%s">\n<br clear=left>\n', imgname);
        catch me
          mydispexception(me);
          error('Could not write image file');
        end
      end
        
      if nargout < 1
        write(hg, stri);
      end
    end
    
    %-------------
    function stri = columns(hg, varargin)
    %-------------
    %COLUMNS Method to insert a table without borders, used for dividing
    %input into columns
    
%     colfunix = [find(cellfun('isclass',varargin,'function_handle'))...
%       length(varargin)+1];
%     if isempty(colfunix)
%       outargs = varargin;
%     else
%       nbroffuns = length(colfunix)-1;
%       outargs = cell(1,nbroffuns);
%       for i = 1:nbroffuns
%         inargs = varargin(colfunix(i):colfunix(i+1)-1);
%         outargs{i} = feval(inargs{:});
%       end
%     end
    
      stri = sprintf('<table border="0" width=%d>\n',hg.pagewidth);
      for i = 1:length(varargin)
        stri = [stri sprintf('<td valign="top">%s</td>\n', varargin{i})]; %#ok<AGROW>
      end
      stri = [stri sprintf('</table>\n')];
      
      if nargout < 1
        write(hg, stri);
      end
    end
    
    %-------------
    function stri = link(hg, ref, text)
    %-------------
    %LINK Method to write a link to a file or URL
      if nargin < 3
        text = ref;
      end
      stri = sprintf('<a href="%s">%s</a>\n', ref, text);
      
      if nargout < 1
        write(hg, stri);
      end
    end
    
    %--------------------------
    function stri = newline(hg)
    %--------------------------
    %NEWLINE Method to insert a line break
    stri = '<br>';
    if nargout == 0
      write(hg,stri);
    end
    
    end
    
    %-------------
    function newpagepos = pagebreakchk(hg, pagepos, nextblock)
    %-------------
    %PAGEBREAKCHK Method to insert a page break
    %Checks if it is time to insert a page break using 'pagepos' parameter
    maxpageheight = 1506; %Approximated from studies of some common browser printing defaults
    newpos = pagepos + nextblock;
    if newpos > maxpageheight
      hg.text('<p style="page-break-before: always"></p>');
      newpagepos = nextblock;
    else
      newpagepos = newpos;
    end
    end
    
    %-------------
    function hg = stop(hg)
    %-------------
    %STOP Method to write necessary HTML footers and close file
      write(hg, html_footer(hg));
      fclose(hg.fid);
    end
    
  end %End of methods
  
  methods (Access = 'private')
    
    %-------------
    function write(hg, str)
    %-------------
    %Writes a string of characters to file
      fprintf(hg.fid, '%s', str);
    end
       
    %--------------------------
    function stri = html_header(hg)
    %--------------------------
    %Return string used as HTML header
      stri = [];

      stri = [stri sprintf('<html>\n')];
      stri = [stri sprintf('<head>\n')];
      stri = [stri sprintf('<title>%s</title>\n', hg.title)];
      stri = [stri sprintf('<meta http-equiv="Content-Type" content="text/html; charset=utf-8">\n')];
      stri = [stri sprintf('<style>\n@media print\n{\ntable {page-break-inside:avoid}\n}\n</style>\n')];
      stri = [stri sprintf('</head>\n\n')];
      stri = [stri sprintf('<body style="width: %dpx">\n', hg.pagewidth)];

    end
    
    %--------------------------
    function stri = html_footer(hg) %#ok<INUSD>
    %--------------------------
    %Return string used as HTML footer
      stri = sprintf('</body>\n</html>\n');
    end
    
  end %End of private methods
  
  methods(Static)
    %-------------
    function stri = conc(varargin)
    %-------------
    %Concatenate two or more function calls
    stri = [varargin{:}];
    end
  end
  
end
    