classdef mywaitbar < handle
  %H = MYWAITBAR(STRING,FIGHANDLE,INTERUPTABLE)
  %
  %New waitbar class for Segment software packages. If fighandles is added (recommended)
  %then it looks in that figure for an axes called waitbaraxes and uses that as holder
  %for the waitbar. Important difference from old waitbars is that the old
  %broken waitbars are cleared if they are present.
  %
  %Methods
  %
  %set(h,fraction); %sets the waitbar
  %
  %close(h); %Close the waitbar
  %
  %isinterrupted(h); %returns true if interrupted

  %Einar Heiberg

  %#ok<*GVMIS>

  properties
    parentfighandle = [];
    fighandle = [];
    docked = false;
    waitbaraxes = [];
    waitbarpatch = [];
    waitbarline = [];
    stoppushbutton = [];
    texthandle = [];
    stri = '';
    fraction = 0;
    interruptable = false;
    interrupted = false;
    visible = true;
  end

  methods

    %---------------------------------------------------------
    function h = mywaitbar(stri,parentfighandle,interruptable)
      %-------------------------------------------------------

      %Store
      h.stri = stri;
      h.fraction = 0;
      h.waitbaraxes = []; %Parse below

      %If stri is not a char then make it invisble
      if ~ischar(stri)
        h.visible = false;
      end

      if (nargin > 1)
        h.parentfighandle = parentfighandle; %default is []
      end

      if nargin > 2
        h.interruptable = interruptable;
      else 
        h.interruptable = false;
      end

      %Parse to find parent axes
      h.waitbaraxes = findwaitbaraxes(get(h.parentfighandle,'Children'));
      h.texthandle = findwaitbartext(get(h.parentfighandle,'Children'));

      %Check if need graphical update
      if ~h.visible
        if ~isempty(h.waitbaraxes)
          h.docked = true;
        end
        return
      end

      if ~isempty(h.waitbaraxes)
        h.docked = true;

        try
          h.waitbaraxes.Parent.Visible = 'on';
        catch
        end

        cla(h.waitbaraxes); %Clear the axes
        h.fighandle = h.parentfighandle;
        
        %Enable that it is working
        myworkon(h.parentfighandle);

        %Create a new patch
        h.waitbarpatch = createpatch(h.waitbaraxes);

      else
        %Create a new figure and waitbar
        h.createnew();
      end

      %Get foreground and background colors
      [bc,fc] = getcolors();

      if ~isempty(h.texthandle)
        %Update text

        if h.interruptable
          set(h.texthandle,...
            'Style','pushbutton',...
            'String',dprintf('Cancel'),...
            'BackgroundColor',bc,...
            'ForegroundColor',fc,...
            'Visible','on',...
            'Callback',@(hself,event) h.stopit());
        else
          set(h.texthandle,...
            'Style','text',...
            'String',h.stri,...
            'Visible','on',...
            'Callback','');
        end
      end

      %Force update
      drawnow();

    end

    %--------------------
    function createnew(h)
      %------------------
      %Creates a new function

      logdisp('Create new waitbar figure')


      % set possible height and width, so always the whole figure is visible
      reqw = 600;
      reqh = 50;

      [px, py, reqw, reqh] = helperfunctions('getfigposition',reqw,reqh,h.parentfighandle);

      %Get foreground and background colors
      [bc,fc] = getcolors();

      %Create the figure
      h.fighandle = figure('Name',h.stri,...
        'Color',bc,...
        'MenuBar','none',...
        'NumberTitle','off',...
        'Position',[px py reqw reqh],...
        'CloseRequestFcn',@(hself,event) h.stopit()); 
      setupicon(h.fighandle);

      %Create axes
      if h.interruptable
        h.waitbaraxes = axes(h.fighandle,'Position',[0.25 0.25 0.7 0.5],'Color',bc);
        h.stoppushbutton = uicontrol(...
          h.fighandle,...
          'Style','pushbutton',...
          'Position',[15 15 90 25],...
          'String',dprintf('Cancel'),...
          'BackgroundColor',bc,...
          'ForegroundColor',fc,...
          'Visible','on',...
          'Callback',@(hself,event) h.stopit());
      else
        h.waitbaraxes = axes(h.fighandle,'Position',[0.05 0.25 0.9 0.5],'Color',bc);
      end

      %--- Add patch etc to axes
      h.waitbarpatch = createpatch(h.waitbaraxes);

      %Make it current
      figure(h.fighandle);

    end

    %----------------------------
    function i = isinterrupted(h)
      %--------------------------
      %Returns true if interrupted
      i = h.interrupted;
    end

    %---------------------
    function set(h,f,stri)
      %-------------------
      %Sets the waitbar with the fraction f, and optionally also updates
      %the string

      if ~h.interrupted
        if nargin < 2
          f = 0;
        end
        f = min(f,1);

        %Store string if presented
        if nargin > 2
          h.stri = stri;
          set(h.texthandle,'String',h.stri);
        end

        %Update patch coordinates
        if isgraphics(h.waitbarpatch)
          set(h.waitbarpatch,'XData',[0;0;f;f]);
          drawnow;
        else
          h.waitbarpatch = createpatch(h.waitbaraxes);
          set(h.texthandle,'Visible','on','Style','text')
          drawnow;
        end
      end

    end

    %----------------
    function close(h)
      %--------------
      %Close the waitbar

      if h.docked
        cla(h.waitbaraxes); %Clear the axes
        set(h.texthandle,... %Clear the text and ensure it is a text
          'Style','text',...
          'Visible','off');

        %Disable that it is working
        myworkoff(h.parentfighandle);

      else
        h.interruptable = true; %now we force it
        h.stopit();
      end

    end

    %-----------------
    function stopit(h)
      %---------------

      try
        if h.interruptable
          %--- It is an interruptable type
          h.interrupted = true;
          if ~h.docked
            delete(h.fighandle); %Delete it if docked
          else
            cla(h.waitbaraxes); %Clear axis if not docked
            set(h.texthandle,... %Clear texthandle
              'Style','text',...
              'Visible','off');

            %Disable that it is working
            myworkoff(h.parentfighandle);

          end
        else

          if yesno('This operation is not interruptible, are you sure you want to force stop?')
            logdisp('Forced interruption of waitbar');
            h.interrupted = true;

            %Disable that it is working
            try
              myworkoff(h.parentfighandle);
            catch
            end
            
            delete(h.fighandle); %Delete it
          end

        end
      catch me
        mydispexception(me)
        return
      end

    end

  end %End of methods

end %End of classdef

%--------------------------------
function ax = findwaitbaraxes(hs)
%--------------------------------
%Find waitbaraxes, if not found returns empty

ax = []; %Mark not found

for loop = 1:length(hs)
  type = get(hs(loop),'Type');
  switch type
    case 'axes'
      if isequal(get(hs(loop),'Tag'),'waitbaraxes')
        ax = hs(loop);
        break; %Hurray found it
      end
    case 'uipanel'
      ax = findwaitbaraxes(get(hs(loop),'Children'));
      if ~isempty(ax)
        break; %Hurray found it
      end
  end
end

end %end of function findwaitbaraxes

%-----------------------------------
function texth = findwaitbartext(hs)
%---------------------------------
%Find waitbaraxes, if not found returns empty

texth = []; %Mark not found

for loop = 1:length(hs)
  type = get(hs(loop),'Type');
  switch type
    case 'uicontrol'
      if isequal(get(hs(loop),'Tag'),'waitbartext')
        texth = hs(loop);
        break; %Hurray found it
      end
    case 'uipanel'
      texth = findwaitbartext(get(hs(loop),'Children'));
      if ~isempty(texth)
        break; %Hurray found it
      end
  end
end

end %end of function findwaitbaraxes

%----------------------------
function ph = createpatch(ax)
%----------------------------
%Create patch

axis(ax,'normal');
waitbarcolor = 'g'; %green
ph = patch([0 0 0 0],[0 1 1 0],[1 1 1 1],'Parent',ax,'FaceColor',waitbarcolor);
hold(ax,'on');
rectangle(ax,'Position',[0 0 1 1],'Curvature',0.1,'LineWidth',4,'EdgeColor',waitbarcolor);
hold(ax,'off');
xlim(ax,[0 1]);
axis(ax,'off');

end

%---------------------------
function [bc,fc] = getcolors
%---------------------------
global DATA

%get color settings
if  ~isempty(DATA)
  bc = DATA.GUISettings.BackgroundColor;
  fc = DATA.GUISettings.ForegroundColor;
else
  % default values
  bc = [0.9400 0.9400 0.9400];
  fc = [0 0 0];
end

end
