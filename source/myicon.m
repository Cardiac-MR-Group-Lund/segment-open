classdef myicon < handle 
  %Inherits from handles to get persistent objects.
  %MYICON add icon to MYICONPLACEHOLDER
% Constructor takes (axes,image,text,exestring,dentstick (Says if the button stays indented afterpress),indentimage)

%Klas Berggren & Einar Heiberg
%#ok<*GVMIS> 
  
  properties (SetAccess = public)
    name='';
    dropdowniconcell={};
    dropdowniconholder=[];
    dropdownpanel = []; %dropdown axes are placed in panels
    dropdownaxes=[];    
    parentobj=[];
    mouseovertext=[];
    cdataOriginal = [];
    cdata = [];
    cdataDisplay=[];
    cdataClicked=[];
    cdataIndent=[];
    cdataDisabled=[];
    enabled=1;
    execute=[];
    isindented=0;
    type=0
    group=1;
    isFDACleared = false;
  end
  
  methods
    function g = myicon(name,parentobj,cdata,text,execute,type,group,cdataIndent,isfdacleared,themebackground) %dropdowniconcell)%cdataClicked,cdataDisplay)
      %Constructor
      %-name is name of the icon
      %-parent object is icon
      %-cdata is icon cdata
      %-text is the tooltip text
      %-execute is the callback
      %-type is the type. 0=can not be indented (example zoom). 
      %                   1=indented, but only one in the group (this is the default)
      %                   2=can be indented no group behaviour
      %-group is the group number, only used when type is 1
      %-cdataIndent if it is set to indented
      %-isfdacleared if false then it is not visible when fda is on.
      
      if nargin==1
        text=[];
      end
      
      if nargin<6
        g.type=1;
      else
        g.type=type;
      end
      
      if nargin<7
        g.group=1;
      else
        g.group=group;
      end

      if nargin < 10
        themebackground = parentobj.themebackground;
      end
      
      %Set properties
      g.name=name;
      g.parentobj=parentobj;
      g.mouseovertext=text;
      g.cdataOriginal = cdata;
      if themebackground
        g.cdata = g.generatebackground(cdata);
      else
        g.cdata = cdata;
      end
      g.cdataDisplay = g.cdata;
      g.execute=execute;
      
      if nargin<8 || isempty(cdataIndent)
        g.generateclickeddisabledandindent(themebackground)
      else
        g.generateclickeddisabledandindent(themebackground,cdataIndent)
      end

      if nargin > 8
        g.isFDACleared = isfdacleared;
      else
        global DATA %#ok<TLEV> 
        if ~strcmp(DATA.ProgramName,'Segment') %FDA cleared icons are implemented only for Segment
          g.isFDACleared = true;
        end
      end
      %generate all stuff needed for dropdowniconholders
%       if nargin==9
%         g.dropdowniconcell = dropdowniconcell;
%         fig = g.parentobj.axeshandle.Parent.Parent;
%         p = uipanel('parent',fig,'units','normalized','position',[0 0 1 1],'BorderType','none','Visible','off');
%         ax = axes('parent',p,'units','normalized','position',[0 0 1 1]);
%         g.dropdowniconholder = myiconplaceholder(ax);
%         add(g.dropdowniconholder,g.dropdowniconcell);
%         g.dropdownpanel.Visible = 'off';
%         g.dropdownpanel = p;
%         g.dropdownaxes = ax; 
        
% % %         fig = g.parentobj.axeshandle.Parent.Parent;
% % %         p = uipanel('parent',fig,'units','normalized','position',[0 0 0 0],'BorderType','none','Visible','off');
% % %         ax = axes('parent',p,'units','normalized','position',[0 0 1 1]);
       
% %         add(g.dropdowniconholder,g.dropdowniconcell)
%         g.dropdownpanel.Visible = 'off';
%       end
      if g.type ~= 0 && g.type ~= 1 && g.type ~= 2 && g.type ~= 3
        g.enabled = 0;
        g.cdataDisplay=g.cdataDisabled;
      end
      
    end
    
    %--------------------------------------------------
    function generateclickeddisabledandindent(varargin)
      global DATA
      g = varargin{1};
      dobackground = varargin{2};
      
      % set background mask when all RGB channels have value  = 240
      backgroundmask = all(g.cdataOriginal == 240,3);
      msk = repmat(backgroundmask, [1 1 3]);

      %Generate CLICKED icon
      tmp = g.cdata;
      tmp(msk) = 100;
      g.cdataClicked = uint8(tmp);


      %Generate DISABLE icon
      if ~dobackground
        %Gray out the icon
        disabledicon = rgb2gray(g.cdataOriginal)/2;
        %Restore background
        disabledicon = disabledicon + 120; %was (240 - max(max(disabledicon)));
        %Assign
        g.cdataDisabled = uint8(repmat(disabledicon,1,1,3));
      else

        % convert the non-background region to grayscale. First, get
        % parameters param1 and param2 for a better scaling based on the GUI background.
        backgroundbrightness = colorfunctions.getbrightness(DATA.GUISettings.BackgroundColor);
        if backgroundbrightness < 0.5
          %dark background
          param1 = 3;
          param2 = 10;
        else
          %light background
          param1 = 2;
          param2 = 120;
        end
        grayimage = repmat(rgb2gray(g.cdataOriginal)/param1 + param2, [1 1 3]);

        % convert mask to uint8 and to RGB
        rgbbackgroundmask = repmat(uint8(backgroundmask), [1 1 3]);
        nonbackgroundmask = repmat(uint8(~backgroundmask), [1 1 3]);

        % apply masks to the image
        g.cdataDisabled = g.cdata.* rgbbackgroundmask + grayimage.* nonbackgroundmask;
      end

      %Generate INDENT icon
      if nargin < 3
        tmp = g.cdataOriginal;
        if dobackground
          newcolour = 120;
        else
          newcolour = 200;
        end

        tmp(msk) = newcolour;
        g.cdataIndent = uint8(tmp);
      else
        %Allows user to specify how indented button should look like
        if dobackground
          usericon = g.generatebackground(varargin{3});
        else
          usericon = varargin{3};
        end
        g.cdataIndent = usericon;
      end
    end

    %--------------------------------------------------
    function cdata = generatebackground(varargin)
      g = varargin{1};
      cdata = varargin{2};
      global DATA

      tmp = cdata;
      %Extract the individual red, green and blue channels
      redchannel = tmp(:, :, 1);
      greenchannel = tmp(:, :, 2);
      bluechannel = tmp(:, :, 3);

      %Change gray background to Segment's background colour
      backgroundcolour = DATA.GUISettings.BackgroundColor*255;
      % Get the gray mask
      graymask = redchannel == 240 & greenchannel == 240 & bluechannel == 240;
      %Change mask to desired colour
      redchannel(graymask) = backgroundcolour(1);
      greenchannel(graymask) = backgroundcolour(2);
      bluechannel(graymask) = backgroundcolour(3);

      %Change black text to Segment's foreground colour
      foregroundcolour = DATA.GUISettings.ForegroundColor*255;
      % Get the black mask
      blackmask = (redchannel == 40 & greenchannel == 40 & bluechannel == 40);
      %Change mask to desired colour
      redchannel(blackmask) = foregroundcolour(1);%54;
      greenchannel(blackmask) = foregroundcolour(2);%60;
      bluechannel(blackmask) = foregroundcolour(3);%72;

      %Recombine the separate channels into a single RGB image
      cdata = cat(3, redchannel, greenchannel, bluechannel);
    end

    %---------------------------
    function highlight(varargin)
      g=varargin{1};
      if g.enabled
        g.cdataDisplay=g.cdataClicked;
      end
    end
    
    %-----------------------------
    function unhighlight(varargin)
      global DATA
      
      g=varargin{1};
      %First check so that button is enabled
      if g.enabled
        %Check that if the button we are over is the same as the the
        %clicked button, if so indent it. If not we want to unhighlight the previous button
        if DATA.Testing
          crit1 = true;
          crit2 = true;
        else
          currenticon=g.parentobj.geticon;
          crit1 = isequal(currenticon,g.parentobj.clickedicon);
          crit2 = isequal(hittest(g.parentobj.fig),g.parentobj.imagehandle); 
        end
        
        if crit1 && crit2
          switch g.type
            case 0 %one click button
              g.cdataDisplay=g.cdata;
            case 1 %dentstick button
              g.cdataDisplay=g.cdataIndent;
              g.isindented=1;
            case 2 %toggle button
              if g.isindented
                g.cdataDisplay=g.cdata;
                g.isindented=0;
              else
                g.cdataDisplay=g.cdataIndent;
                g.isindented=1;
              end
            case 3 %dropdown button
              if g.isindented
                g.cdataDisplay=g.cdata;
                g.isindented=0;
              else
                g.cdataDisplay=g.cdataIndent;
                g.isindented=1;
              end
          end
          %Redundant for most cases but assures that play button gets
          %proper appearance
          g.parentobj.render
          
          %Kill timer in parent so we dont display text
          if ~isempty(g.parentobj.texttimer)
            g.parentobj.displayinfo=0;
            stop(g.parentobj.texttimer);
            delete(g.parentobj.texttimer);
            g.parentobj.texttimer=[];
          end
          if g.type~=3
            feval(g.execute);
          else
            %g.setdropdown
            %               if g.isindented
            %               %do dropdown i.e place new temporary iconaxes which contains
            %               %a dropdown iconholder
            %
            %               %Number of rows and cols in iconholder divide this with the height and width of
            %               %the axes of iconholder also calculate lower left corner of
            %                [iconcornerx,iconcornery] = buttoncorner(g.parentobj,g.name);
            %
            %                rows = size(g.parentobj.iconCell,1);
            %                cols = size(g.parentobj.iconCell,2);
            %                pos =plotboxpos(g.parentobj.axeshandle);%get(g.parentobj.axeshandle,'position');
            %                iconheight = pos(4)/rows; %in gui
            %                iconwidth = pos(3)/cols; %in gui
            %                offsety = size(g.dropdowniconcell,1)*iconheight;
            %
            %               %place new axes on bottom left corner of dropdown icon
            %               newaxespos = [iconcornerx,iconcornery-offsety,iconwidth,offsety+iconheight/2];
            %               fig = get(g.parentobj.axeshandle,'parent');
            %               ax = axes('parent',fig,'units','normalized','position',newaxespos);
            %               g.dropdownaxes=ax; %add to dropdown axes field
            %               %add iconholder to axes and render
            %               g.dropdowniconholder = myiconplaceholder(ax,0,3);
            %               g.parentobj.dropdowniconholders{end+1} = g.dropdowniconholder;
            %               add(g.dropdowniconholder,g.dropdowniconcell)
            %
            %               for i = 1:numel(g.dropdowniconcell)
            %                 g.dropdowniconcell{i}.parentobj.axeshandle
            %               end
            %
            %               else
            %                 g.parentobj.dropdowniconholders([g.parentobj.dropdowniconholders{:}]==g.dropdowniconholder)=[];
            %                 delete(g.dropdownaxes);
            %                 delete(g.dropdowniconholder);
            %                 g.dropdowniconholder = [];
            %                 g.dropdownaxes=[];
            %               end
          end
        else
          g.cdataDisplay=g.cdata;
        end
      end
    end
    
    %-----------------------------
    function showdropdown(varargin)
      g = varargin{1};
      if g.isindented
        if isequal(g.dropdownpanel.Visible,'off')
          %uistack(g.dropdownpanel,'top');
          g.dropdownpanel.Visible = 'on';
          add(g.dropdowniconholder,g.dropdowniconcell)
        end
      else
        if isequal(g.dropdownpanel.Visible,'on')
          g.dropdownpanel.Visible = 'off';          
        end
      end
    end

    %-----------------------------
    function placedropdown(varargin)
      g = varargin{1};
      
        %do dropdown i.e place new temporary iconaxes which contains
        %a dropdown iconholder
        
        %Number of rows and cols in iconholder divide this with the height and width of
        %the axes of iconholder also calculate lower left corner of
        [iconcornerx,iconcornery] = buttoncorner(g.parentobj,g.name);
        
        rows = size(g.parentobj.iconCell,1);
        cols = size(g.parentobj.iconCell,2);
        pos = plotboxpos(g.parentobj.axeshandle);%get(g.parentobj.axeshandle,'position');
        
        panelpos = g.parentobj.axeshandle.Parent.Position;
        %snap buttoncorner to panelbottom
        iconcornery = panelpos(2);
        
        
        iconheight = pos(4)*panelpos(4)/rows; %in gui
        iconwidth = pos(3)*panelpos(3)/cols*1.2; %in gui
        offsety = (size(g.dropdowniconcell,1))*iconheight-0*iconheight/2;
        newaxespos = [iconcornerx,iconcornery-offsety-iconheight/2,iconwidth,offsety+iconheight/2];
        g.dropdownpanel.Position = newaxespos;
       
    end

      
      
%       if isindented
%         %do dropdown i.e place new temporary iconaxes which contains
%         %a dropdown iconholder
%         
%         %Number of rows and cols in iconholder divide this with the height and width of
%         %the axes of iconholder also calculate lower left corner of
%         [iconcornerx,iconcornery] = buttoncorner(g.parentobj,g.name);
%         
%         rows = size(g.parentobj.iconCell,1);
%         cols = size(g.parentobj.iconCell,2);
%         pos = plotboxpos(g.parentobj.axeshandle);%get(g.parentobj.axeshandle,'position');
%         
%         if strcmp(g.parentobj.axeshandle.Parent.Type,'uipanel')
%           panelpos = g.parentobj.axeshandle.Parent.Position;
%           fig = g.parentobj.axeshandle.Parent.Parent;
%           %snap buttoncorner to panelbottom
%           iconcornery = panelpos(2);
%         else
%           panelpos = [0,0,1,1];
%           fig = g.parentobj.axeshandle.Parent;
%         end
%         
%         iconheight = pos(4)*panelpos(4)/rows; %in gui
%         iconwidth = pos(3)*panelpos(3)/cols*1.2; %in gui
%         offsety = (size(g.dropdowniconcell,1))*iconheight-0*iconheight/2;
%         
%         %place new axes on bottom left corner of dropdown icon
%         newaxespos = [iconcornerx,iconcornery-offsety-iconheight/2,iconwidth,offsety+iconheight/2];
%         p = uipanel('parent',fig,'units','normalized','position',newaxespos,'BorderType','none');
%         ax = axes('parent',p,'units','normalized','position',[0 0 1 1]);        
%         g.dropdownpanel = p;
%         g.dropdownaxes = ax; %add to dropdown axes field
%         
%         %add iconholder to axes and render
%         g.dropdowniconholder = myiconplaceholder(ax,0,3);
%         if isempty(g.parentobj.dropdowniconholders)
%           g.parentobj.dropdowniconholders = g.dropdowniconholder;
%         else
%           g.parentobj.dropdowniconholders(end+1) = g.dropdowniconholder;
%         end
%         add(g.dropdowniconholder,g.dropdowniconcell)
%       else
%         
%         if (~isempty(g.parentobj.dropdowniconholders)) && (~isempty(g.dropdowniconholder))
%           g.parentobj.dropdowniconholders(g.parentobj.dropdowniconholders==g.dropdowniconholder)=[];
%           delete(g.dropdownaxes);
%           delete(g.dropdownpanel);
%           delete(g.dropdowniconholder);
%           g.dropdowniconholder = [];
%           g.dropdownpanel = [];
%           g.dropdownaxes = [];          
%         end
%         
%       end
%     end
    
    %------------------------
    function undent(varargin)
      g=varargin{1};
      g.isindented=0;
      if g.enabled
        g.cdataDisplay=g.cdata;
      else
        g.cdataDisplay=g.cdataDisabled;
      end
    end
    
    %------------------------
    function enable(varargin)
      g=varargin{1};
      runFDAVersion = varargin{2}; %check if the GUI is run in FDA mode
      if nargin > 2
        dorendering = varargin{3};
      else
        dorendering = true;
      end
      if ((runFDAVersion && g.isFDACleared) || ~runFDAVersion) && g.type ~= 100 % icon type = 100 means this icon is the whole time disabled, since the user does not have the license
        g.enabled=1;
        g.cdataDisplay=g.cdata;
        if dorendering
          g.parentobj.render
        end
      end
    end
    
    %-------------------------
    function disable(varargin)
      g=varargin{1};
      if nargin > 1
        dorendering = varargin{2};
      else
        dorendering = true;
      end
      g.enabled=0;
      g.cdataDisplay=g.cdataDisabled;
      if dorendering
        g.parentobj.render
      end
    end
    
    %-------------------------
    function setIcon(varargin)
      g=varargin{1};
      g.cdata=varargin{2};
    end
    end
end