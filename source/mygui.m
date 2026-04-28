classdef mygui < handle 
  %Inherits from handles to get persistent objects.
  %MYGUI Class for handeling GUI's in a efficeint manner
%
%G = MYGUI(FIGFILENAME,<'BLOCKING'>);
%
%The purpose of this class is to facilitate to generate 
%GUI's than can also work as a temporary containter data  
%for that so that data transfer between subfunctions becomes 
%simple. If blocking is specified then segment internal blocking 
%figure mechanism will be used.
%
%The function also keeps track of position for the gui to 
%be saved later and to be able to set position of messageboxes, 
%called with for example mywarning, mymsgbox and yesno when 
%calling with last input argument as the mygui object.

%Example:
% First, reserve space for it in DATA Structure:
%
%  DATA.GUI.MyGUIExample = []; %in maingui.m
%
%In your code to initialize:
%  gui = mygui('myguiexample.fig');
%
%In your subfunctions:
%  global DATA
%
%  gui = DATA.GUI.MyGuiExample;
%
%Now you can use gui as an ordinary variable, and everthing 
%is stored in the GUI, and will be available to other 
%subfunctions until the GUI is deleted. Some examples:
%
%  gui.myfieldname = somearray;
%  gui.myfieldname1.myfieldname2.myfieldname3 = 134;
%  temp = gui.myfieldname(2);
%  temp = get(gui.handles.text2,'string');
%
%You can also use the gui to get good alignment of the 
%msgboxes mywarning, myfaield, mymsgbox, mymenu, mywaitbarstart and yesno.
%
% myfailed('An error has occured.',gui);
%
%When closing the gui the close function can be called with
% DATA.GUI.MyGUIExample=close(DATA.GUI.MyGUIExample);
% when closing by using this syntax the position of the gui 
% will be saved in the struct DATA.GUIPositions
%
%The position of the gui can also by saved by 
% saveguiposition(gui);
% this could be used in a resize function
% 
%See also SUBSREF, SUBSASGN.

%Jane Sjogren

%#ok<*GVMIS> 

properties (SetAccess = 'private', Hidden)
  filename = '';
  blocking = false;
  handles = [];
  fig = [];
  ispseudomaximized = false;
  isWindows11 = false;
end
properties (Hidden, Constant)
  mainfignames = {...
    'segment.fig',...
    ['+segmentct' filesep 'segmentct.fig'],...
    ['+segmentmr' filesep 'segmentmr.fig'],...
    };
    %Segment 3DPrint does not need it  ['+segment3dp' filesep 'segment3dp.fig'],...
end

methods

    %--------------
    function g = mygui(filename,blocking)

      global DATA 
      
      %Constructor
      isinvisible = false;
      if nargin>1
        if isequal(blocking,'blocking')
          blocking = true;
        elseif isequal(blocking,'invisible')
          isinvisible = true;
        else
          error('Expected ''blocking''');
        end
      end

      if nargin<2
        blocking = false;
      end

      g.fig = []; %save both as internal class and also as a appdata (later)
      g.filename = filename; %OK to have here since read only
      g.blocking = blocking; %OK to have here since read only

      g.fig = openfig(g.filename,'reuse','invisible'); %set visble off until position is set and handles are stored
      if g.blocking
        blockfig(g.fig);
      end

      if ~isempty(DATA)
        set(g.fig,'Color',DATA.GUISettings.BackgroundColor); %get(0,'defaultUicontrolBackgroundColor'));
      else
      end
      
      setupicon(g.fig); %set up Segment icon
      
      %subsasgn(g,'handles',guihandles(g.fig)); %Store handles as application data
      g = subsasgn(g,'fig',g.fig);
      %g.handles.fig = g.fig;
      g.checkwindowsversion;

      set(g.fig, 'SizeChangedFcn', @g.sizechanged);

      if ~isempty(DATA)
        if ~DATA.issegment3dp
          setguiposition(g);%set saved position must be set after handles have been set
        end
      end
      
      g.handles = guihandles(g.fig); %Store handles as application data %Moved this line down. JS, EH:
      if ~strcmpi(g.fig.Name,'PACS preferences')
        translation.translatealllabels(g.fig); %Translate labels into preferred language
      end
      if not(isempty(DATA)) && ~isempty(DATA.Pref) && ~strcmpi(g.fig.Name,'preffig') % if non-existent then the default color is use anyway
        %don't use background color for preffig. 
        setinterfacecolor(g.fig); %set background color and text color for all objects in alla interfaces
      end
%       lgcolor = [0.94 0.94 0.94];
%       try set(DATA.Handles.iconuipanel,'BackgroundColor',lgcolor); catch, end
%       try set(DATA.Handles.ribbonuipanel,'BackgroundColor',[0.94 0.94 0.94]); catch, end
            
      if ~isinvisible
        set(g.fig,'visible','on');%set visble off until position is set and handles are stored
      end
      flushlog;
    end

    %------------------
    function display(g)
      %DISPLAY method
      fprintf('MYGUI object %s\n',g.filename);
      try
        c = fieldnames(getappdata(g.fig));
      catch
        fprintf('Deleted object\n')
        return
      end

      for loop = 1:length(c)
        if ~ismember(c{loop},{'GUIDEOptions','lastValidTag','Listeners','SavedVisible'})
          type = class(getappdata(g.fig,c{loop}));
          fprintf('  %s: %s\n',c{loop},type);
        end
      end
    end

    %---------------------------
    function g = subsasgn(g,s,v)
      %SUBSASGN Method to set properties in MYGUI object
      %
      %G = SUBSASGN(G,S,V)    Where G is GUI-object, S subsasgn struct, V is value
      %
      %G = SUBSASGN(G,STRI,V) Where G i s GUI-object, STRI is string, V is value
      %
      %See also MYGUI.

      %Einar Heiberg with great help from Jane with debugging.

      if isstruct(s)
        %STRUCT

        if length(s)==1
          %Only one, i.e gui.field = value
          if isequal(s.type,'.')
            if isequal(s.subs,'handles')
              g.handles = v;
            else
              setappdata(g.fig,s.subs,v);
            end
          else
            error('Invalid syntax, expected ''gui.field = value''.');
          end
        else
          %Longer than one more complex interaction
          if ~isequal(s(1).type,'.')
            error('Invalid syntax, expected ''gui.field...''');
          end
          if isequal(s(1).subs,'handles')
            g = builtin('subsasgn',g,s,v);
          else
            %Normal case
            if isappdata(g.fig,s(1).subs)
              %Already exists
              temp = getappdata(g.fig,s(1).subs); %Extract it
              temp = subsasgn(temp,s(2:end),v);
              setappdata(g.fig,s(1).subs,temp); %Store it
            else
              %Does not exist
              temp = [];
              temp = subsasgn(temp,s(2:end),v);
              setappdata(g.fig,s(1).subs,temp); %Store it
            end
          end
        end
      else
        %STRING
        setappdata(g.fig,s,v);
      end
    end
    
    %-----------------------------
    function p = subsref(g,s)
      %SUBSREF Method to extract parameters from MYGUI object

      %Einar Heiberg, with great help from Jane with debugging.
      p = [];

      %Check if first is . (has to be).
      if ~isequal(s(1).type,'.')
        error('Invalid syntax. Valid syntax is gui.name or gui.handles.handle');
      end

      switch s(1).subs
        case 'handles'
          if length(s)>1
            p = subsref(g.handles,s(2:end));
          else
            p = g.handles;
          end
        otherwise
          if isstruct(s)
            %struct
            if ~ishandle(g.fig)
              return
            end

            if length(s)==1
              %Only one => has to be gui.field
              if isappdata(g.fig,s.subs)
                p = getappdata(g.fig,s.subs);
              else
                error('Field ''%s'' does not exist.',s.subs);
              end
            else
              %Complex
              temp = getappdata(g.fig,s(1).subs);

              p = subsref(temp,s(2:end));
            end
          else
            %char
            p = getappdata(g.fig,s);
          end
      end
    end
    
    %-----------------------------
    function g = close(g)
      %CLOSE Method to close MYGUI objects

      %Einar Heiberg
      %also saves the guiposition to DATA
      saveguiposition(g);
      delete(g.fig);
      g = [];
      flushlog;
    end
    
    %-----------------------------
    function myfailed(stri,g)

      myfailed(stri,g.fig);
    end
    
    %----------------------------
    function mywarning(stri,g)

      mywarning(stri,g.fig)
    end
    %-------------------------------------
    function mymsgbox(stri,title,g)

      mymsgbox(stri,title,g.fig);
    end
    
    %---------------------------------------------
    function temp = mymenu(header,varargin)
      n= length(varargin);
      g = varargin{n};
      n = length(varargin{1});
      temp=cell(n+1,1);
      for loop=1:n
        temp{loop}=varargin{1}{loop};
      end
      temp{n+1}=g.fig;
      
      temp=mymenu(header,temp);
    end

    %----------------------------------------------
    function myadjust(h,g)

      myadjust(h,g.fig);
    end
    
    %---------------------------------------------
    function h = mywaitbarstart(iter,stri,ignorelimit,g)
      
      h = mywaitbarstart(iter,stri,ignorelimit,g.fig);
    end

    %-------------------------------------------------
    function result = yesno(s,arg,g)

      result = yesno(s,arg,g.fig);
    end

    %--------------------------------------------
    function mymaximise(g)
      %--------------------------------------------
      % Maximise GUI
      g.fig.WindowState = 'maximized';
%       saveguiposition(g);
    end
    
    %--------------------------------------------
    function sizechanged(g,~,~)
      %--------------------------------------------
      % adjust figure size to pseudomaximized for Windows 11 
     
      if g.isWindows11

        if ~strcmp(version('-release'),'2022a')
          return
        end

        warning off MATLAB:ui:javaframe:PropertyToBeRemoved
        jFrame = get(handle(g.fig), 'JavaFrame'); %#ok<JAVFM>
        javafig = jFrame.getFigurePanelContainer;
        if ~javafig.isShowing
          return % do not do settings when not visible on screen
        end

        origresizefunc = g.fig.SizeChangedFcn;

        % Check if the figure is maximized
        ismaximized = jFrame.isMaximized;

        startpos = javafig.getLocationOnScreen;
        x = startpos.getX + 1;
        y = startpos.getY;
        figwidth = javafig.getWidth;
        figheight = javafig.getHeight;

        currentscreensize = g.getcurrentscreensize(x,y,figwidth,figheight);

        % Define the new size for the figure when maximized
        if ~g.ispseudomaximized &&  ismaximized 
          drawnow limitrate;
          guipos = get(g.fig,'position');
          set(g.fig,'SizeChangedFcn',[]);
          g.ispseudomaximized = true;        
          set(g.fig,'position',[guipos(1) guipos(2) guipos(3)-1 guipos(4)-1]);
          drawnow;
          set(g.fig,'SizeChangedFcn',origresizefunc);
        elseif g.ispseudomaximized &&  ismaximized
          drawnow limitrate;
          set(g.fig,'SizeChangedFcn',[]);
          g.ispseudomaximized = false;

          guipos = g.getsavedguipositions;
          
          newx = currentscreensize(1)+guipos(1)*currentscreensize(3);
          newy = currentscreensize(2)+guipos(2)*currentscreensize(4);
          newwidth = guipos(3)*currentscreensize(3);
          newheight = guipos(4)*currentscreensize(4);
          set(g.fig,'position',[newx newy newwidth newheight]);
          drawnow;
          set(g.fig,'SizeChangedFcn',origresizefunc);
        else
          g.saveguiposition;
        end
        if ismember(g.filename,g.mainfignames)
          segment_main('rendericonholders');
          drawnow limitrate; 
        end
        warning on MATLAB:ui:javaframe:PropertyToBeRemoved

      end %End of Windows 11

    end

    %----------------------------------
    function screensize = getcurrentscreensize(~,x,y,figwidth,figheight) %#ok<INUSL> 
      %----------------------------------
      % get screensize of the current monitor
      
      monitorPositions = get(0, 'MonitorPositions'); %positions of all monitors
      y = figheight-y; % figheight - y is conversion from java coordinates into matlab's

      % Check which monitor contains the figure
      screennum = [];
      for snum = 1:size(monitorPositions, 1)
        screenpos = monitorPositions(snum, :);
        if x >= screenpos(1) && ...
            y >= screenpos(2) && ... 
            x  <= screenpos(1) + screenpos(3) && ...
            y  <= screenpos(2) + screenpos(4)
          screennum = snum;
          break;
        end
      end

      if ~isempty(screennum)
        screensize = monitorPositions(screennum, :);
      else
        screensize = get(0, 'ScreenSize');
      end
    end

    %----------------------------------
    function guipos = getsavedguipositions(g)
      %----------------------------------
      % get positions saved for the current figure

      global DATA
      guipos = [];
      fname = g.filename;
      for loop = 1:length(DATA.GUIPositions)
        if isequal(DATA.GUIPositions(loop).FileName,fname)
          guipos = DATA.GUIPositions(loop).Position;
          break
        end        
      end
      maxpos = 0.95;
      minpos = 0.05;
      if ~isempty(guipos)
        tcypos = guipos(2)+guipos(4);
        tcxpos = guipos(1)+guipos(3);
        if guipos(2) > maxpos || tcypos > maxpos || guipos(2) < minpos || tcypos < minpos || ...
            guipos(1) > maxpos || tcxpos > maxpos || guipos(1) < minpos || tcxpos < minpos
          guipos = [0.1 0.1 0.8 0.7];
        end
      else
        guipos = [0.1 0.1 0.8 0.7];
      end
    end
    
    %----------------------------------
    function saveguiposition(g)
      global DATA
      if g.ispseudomaximized || isempty(DATA)
        return
      end

      units = get(g.fig,'units');
      set(g.fig,'units','normalized');

      guipos = get(g.fig,'position');

      %search for gui and save Position
      found=0;
      for loop=1:length(DATA.GUIPositions)
        if isequal(DATA.GUIPositions(loop).FileName,g.filename)
          DATA.GUIPositions(loop).Position = guipos;
          found=1;
        end
      end
      if not(found)
        numberguis=length(DATA.GUIPositions);
        DATA.GUIPositions(numberguis+1).FileName = g.filename;
        DATA.GUIPositions(numberguis+1).Position = guipos;
      end

      set(g.fig,'position',guipos);
      set(g.fig,'units',units);
    end
    
    %----------------------------
    function setguiposition(g)
      global DATA

      set(g.fig,'units','normalized');

      if isempty(DATA)
        return;
      end
      
      fname = g.filename;

      %search for saved position
      guipos = [];
      segmentpos = [];
            
      for loop=1:length(DATA.GUIPositions)
        if isequal(DATA.GUIPositions(loop).FileName,fname)
          guipos = DATA.GUIPositions(loop).Position;
        end
        if ismember(DATA.GUIPositions(loop).FileName,g.mainfignames)
          segmentpos = DATA.GUIPositions(loop).Position;
        end
      end
      
      if isempty(guipos) && ~isempty(segmentpos)
        guipos = get(g.fig,'position');
        guipos(1:2) = segmentpos(1:2)+(segmentpos(3:4)-guipos(3:4))/2;
      end
      
      if isempty(guipos)
        g.saveguiposition();
        set(g.fig,'units','pixels');
        return
      end
      
      %Make sure the GUI do not end up off screen
      tcypos = guipos(2)+guipos(4);
      tcxpos = guipos(1)+guipos(3);
      maxvalue = 0.95;
      minvalue = 0.05;

      if guipos(2) > maxvalue || tcypos > maxvalue || guipos(2) < minvalue || tcypos < minvalue || ...
          guipos(1) > maxvalue || tcxpos > maxvalue || guipos(1) < minvalue || tcxpos < minvalue
        if ismember(DATA.GUIPositions(loop).FileName,g.mainfignames) || isempty(segmentpos)
          guipos = [0.1 0.1 0.8 0.7];
          set(g.fig,'position',guipos);
        else
          guipos = get(g.fig,'position');
          guipos(1:2) = segmentpos(1:2)+(segmentpos(3:4)-guipos(3:4))/2;
          % ensure new figure coordinates inside main segment figure
          guipos(1) = max([guipos(1), segmentpos(1)]); % maximal x pos
          guipos(2) = max([guipos(2), segmentpos(2)]); % maximal y pos
          guipos(3) = min([guipos(3), segmentpos(3)]); % minimal height
          guipos(4) = min([guipos(4), segmentpos(4)]); % minimal width
          set(g.fig,'position',guipos);
        end
      else
        set(g.fig,'position',guipos);
      end
      set(g.fig,'units','pixels');

      saveguiposition(g);
    end

    %-------------------------------------
    function checkwindowsversion(g)
      %-------------------------------------
      g.isWindows11 = helperfunctions('iswindows11');
    end

    %-----------------------
    function requestfocus(g)
      %---------------------
      % request focus for figure, useful when using arrows that they do not
      % trigger callback all the time
      warning off
      javafig = get(g.fig,'JavaFrame'); %#ok<JAVFM>
      javafig.requestFocus;
      warning on
    end

    %-------------------------------------
    function didclose = closeifhasno(g,no)
      %-------------------------------------
      %Close GUI if it is associated with image stack no
      didclose = false;
      p = getappdata(g.fig);
      props = fieldnames(p);
      for i = 1:numel(props)
        if length(props{i}) >= 2 && ...
            (strcmpi(props{i}(end-1:end),'no') || ...
            strcmpi(props{i}(1:2),'no') && length(props{i}) == 3) ...
            && isequal(p.(props{i}),no) || ~isempty(regexp(props{i},'taggroup'))
          close(g);
          didclose = true;
          return
        end
      end
    end

  end%end of methods
  
  %------------------------------
  methods (Access = 'private')

    %----------------------------
    function setdlgposition(g,h)
      gunits=get(g.fig,'units');
      hunits=get(h,'units');
      set(h,'units','pixels');
      set(g.fig,'units','pixels');

      %get position for gui and the handle h
      gpos=get(g.fig,'position');%guiposition [x(g) y(g) width(g) height(g)]
      hpos=get(h,'position');%handleposition [x(h) y(h) width(h) height(h)]

      %calculate new position for the handle h
      temp=gpos(1:2)+(gpos(3:4)-hpos(3:4))/2;%[x y](g) + ([width height](g) -[width height](h))/2
      hpos(1:2)=temp;

      %set position
      set(h,'position',hpos);

      set(g.fig,'units',gunits);
      set(h,'units',hunits);
    end
    
    
  end %End private methods

end
