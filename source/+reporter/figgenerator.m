%FIGGENERATOR Class for generating figures containing text documents
%Position: [vänsterkant underkant bredd höjd]

classdef figgenerator < handle
  
  properties
    filename = '';
    pathname = '';
    title = '';
    allfigs = [];
    figno = [];
    lineht = 0;
    pos = 0;
    pagewidth = 0;
    stf = struct;
  end
  
  methods
    
    %---------------------------------------------------
    function fg = figgenerator(title, pathname, filename, pagewidth, lineht)
    %---------------------------------------------------
    %Constructor
    
    if nargin < 5
      lineht = 14;
      if nargin < 4
        pagewidth = 800;
        if nargin < 3
          filename = [title '_%04.0f.jpg'];
          if nargin < 2
            pathname = pwd;
          end
        end
      end
    end
    fg.title = title;
    fg.filename = filename;
    fg.pathname = pathname;
    fg.pagewidth = pagewidth;
    fg.lineht = lineht;
    fg.stf.fontname = 'courier new';
    end
    
    %---------------------------
    function start(fg, figtitle)
    %---------------------------
    %START Method to open a new figure
    
    fg.figno = figure;
    fg.allfigs = [fg.allfigs fg.figno];
    axis ij
    %axis off
    hold on
    figpos = get(fg.figno,'position');
    figpos(3) = fg.pagewidth;
    set(fg.figno,'position',figpos);
    fg.pos = fg.lineht/2;
    if nargin == 2
      title(figtitle);
    end
    end
    
    %----------------
    function zoom(fg)
    %----------------
    %ZOOM Zoom to current page size
    pageht = fg.pos+fg.lineht;
    axis off
    axis([0 fg.pagewidth 0 pageht]);
    figpos = get(fg.figno, 'position');
    figpos(4) = pageht + 2*fg.lineht;
    set(fg.figno, 'position', figpos);
    end
    
    %-----------------
    function hline(fg)
    %-----------------
    %HLINE Method to write a horizontal line
    
    return %Might not get this to work properly
    figure(fg.figno)
    fg.newline;
    plot([0 fg.pagewidth],[fg.pos fg.pos],'k');
    fg.newline;
    
    end
    
    %---------------------------------------------------------
    function call = table(fg, content, boldcells, width, lpos)
    %---------------------------------------------------------
    %TABLE Method to insert a table with content specified by a cell

    if nargin < 5
      lpos = 0;
      if nargin < 4
        width = round(0.625*fg.pagewidth);
        if nargin < 3
          boldcells = zeros(size(content));
          if nargin < 2
            myfailed('Too few input arguments.');
          end
        end
      end
    end
    
    if nargout > 0
      call = @(lpos)table(fg,content,boldcells,width,lpos);
      return
    end
    
    startpos = fg.pos;
    
    for j = 1:size(content, 2)
      maxwidth = max(cellfun(@length,content(:,j)));
      lpos = [lpos lpos(end)+(maxwidth+3)*fg.lineht*10/14];
    end
    
    tpos = startpos;
    for i = 1:size(content, 1)
      for j = 1:size(content, 2)
        if ~isempty(content{i,j})
          st = fg.stf;
          if boldcells(i, j) == 2
            st.color = [1 0 0];
          elseif boldcells(i, j)
            st.fontweight = 'bold';
          end
          if ~isempty(regexp(content{i,j},'(&sup2;)', 'once'))
            constri = regexprep(content{i,j},'(&sup2;)','^2');
            tpos = tpos + fg.lineht/4;
            text(lpos(j), tpos, constri, st);
          else
            text(lpos(j), tpos, content{i,j}, st);
          end
        end
      end
      tpos = tpos + fg.lineht;
      fg.newline;
    end
    
    end
    
    %-------------------------------------------
    function call = headline(fg, headtext, lpos)
    %-------------------------------------------
    %HEADLINE Method to write a headline in bold
    if nargin < 3
      lpos = 0;
    end
    
    if nargout > 0
      call = @(lpos)headline(fg, headtext, lpos);
      return
    end
    
    fg.newline;
    text(lpos, fg.pos, headtext, 'fontsize', 16, 'fontweight', 'bold');
    fg.newline;
    fg.newline;
    end
    
    %-----------------------------------------------------------
    function call = text(fg, paragraph, alignment, lpos, weight)
    %-----------------------------------------------------------
    %TEXT Method to write a paragraph of text

    if nargin < 5
      weight = 'normal';
      if nargin < 4
        lpos = 0;
        if nargin < 3
          alignment = 'left';
        end
      end
    end
    
    if nargout > 0
      call = @(lpos)text(fg, paragraph, alignment, lpos);
      return
    end
    
    st = fg.stf;
    if nargin == 3
      st.horizontalalignment = alignment;
      if strcmp(alignment,'right')
        lpos = fg.pagewidth;
      end
    end
    h = text(lpos, fg.pos, paragraph, st);
    set(h,'FontWeight',weight);
    fg.newline;
    end
    
    
    %------------
    function gsstri = ftext(fg)
    %------------
    %FTEXT Method to write a paragraph of text with formatted headlines.
    global SET
    %stri = [];
    gsstri = []; %Gensvar string (unformatted)
    cr = [];
    prow = 1; %Variable to keep track of row number in current paragraph
    %Convert from Matlab string matrix to viewable text
    for row = 1:size(SET(1).Report.Comments,1)
      crprev = cr;
      cr = deblank(SET(1).Report.Comments(row,:));
      gsstri = [gsstri cr sprintf('\n')]; %#ok<AGROW>
      if isempty(cr) %Insert new paragraph if a line is blank
        %stri = [stri sprintf('\n</p>\n<p>\n')]; %#ok<AGROW>
        if ~isempty(crprev)
          fg.text(crprev);
        end
        fg.newline;
        prow = 1;
      else
        if prow == 2 %Create headline according to rules of formatting (see manual)
          if strcmp(upper(crprev),crprev)
            %stri = regexprep(stri,crprev,sprintf('<h3>%s</h3>',crprev)); %superheadline
            fg.text(crprev,'left',0,'bold');
            %fg.text(cr);
            prow = 1;
          else
            %stri = regexprep(stri,crprev,sprintf('<h4>%s</h4>',crprev)); %subheadline
            fg.text(crprev,'left',0,'bold');
            %fg.text(cr);
          end
        elseif prow > 2
          fg.text(crprev);
          fg.newline;
          %stri = [stri sprintf('<br>\n')]; %#ok<AGROW>
        end
        %stri = [stri sprintf('%s\n',cr)]; %#ok<AGROW>
        prow = prow + 1;
      end
    end
    fg.text(cr);
      
    end
       
    %--------------------------------------------------
    function call = image(fg, imgname, imgsource, lpos)
    %--------------------------------------------------
    %IMAGE Method to add an image. imgsource can be a 
    %location or an image matrix.
    
    if nargout > 0
      call = @(lpos)image(fg, imgname, imgsource, lpos);
      return
    end
    
    if nargin < 4
      lpos = 0;
    end
    
    if ischar(imgsource)
      try
        imgsource = imread(imgsource);
      catch me
        mydispexception(me);
        error('Could not read image file');
      end
    end
    sz = size(imgsource);
    newpos = fg.pos + sz(1);
    fg.newline;
    imagesc([lpos lpos+sz(2)], [fg.pos newpos], imgsource);
    fg.pos = newpos;
    fg.newline;
    end
    
    %-----------------------------
    function columns(fg, varargin)
    %-----------------------------
    %COLUMNS Method for dividing input into columns
    startpos = fg.pos;
    endpos = startpos;
    nbrcols = length(varargin);
    step = fg.pagewidth / nbrcols;
    for i = 1:nbrcols
      fg.pos = startpos;
      lpos = (i-1)*step;
      if iscell(varargin{i})
        for j = 1:length(varargin{i})
          if iscell(varargin{i}{j})
            for k = 1:length(varargin{i}{j})
              feval(varargin{i}{j}{k}, lpos)
            end
          elseif ~isempty(varargin{i}{j})
            feval(varargin{i}{j}, lpos);
          end
        end
      elseif ~isempty(varargin{i})
        feval(varargin{i}, lpos)
      end
      endpos = max(fg.pos, endpos);
    end
    fg.pos = endpos;
    
    end
    
    %--------------------------
    function call = newline(fg)
    %--------------------------
    %NEWLINE Method to insert a line break
    if nargout > 0
      call = @()newline(fg);
      return
    end
    
    fg.pos = fg.pos + fg.lineht;
    end
    
    %---------------------------------------------------------
    function newpagepos = pagebreakchk(fg, pagepos, nextblock)
    %---------------------------------------------------------
    %PAGEBREAKCHK Method to insert a page break
    %Checks if it is time to insert a page break using 'pagepos' parameter
    maxpageheight = 1506; %Approximated from studies of some common browser printing defaults
    newpos = pagepos + nextblock;
    if newpos > maxpageheight
      fg.zoom;
      fg.start;
      fg.lineht=15;
      axis equal
      set(gcf,'units','normalized','outerposition',[0 0 1 1])
      %xlim([0,1000])
      newpagepos = nextblock;
      %axis equal
    else
      newpagepos = newpos;
    end
    end
    
    %----------------
    function stop(fg)
    %----------------
    %STOP Save all figures to image files, then close them
    fg.zoom;
    for i = 1:numel(fg.allfigs)
      saveas(fg.allfigs(i),fullfile(fg.pathname,sprintf(fg.filename,i)))
    end
    close(fg.allfigs)
    end
  
  end
  
  methods(Access = 'private')
    
    
  end
  
  methods(Static)
    %-------------
    function call = conc(varargin)
    %-------------
    %Concatenate two or more function calls
    call = varargin;
    end
  end
  
end