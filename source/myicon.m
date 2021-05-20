%MYICON add icon to MYICONPLACEHOLDER
% Constructor takes (axes,image,text,exestring,dentstick (Says if the button stays indented afterpress),indentimage)
%Klas Berggren & Einar Heiberg

classdef myicon < handle %Inherits from handles to get persistent objects.
  
  properties (SetAccess = public)
    name='';
    dropdowniconcell={};
    dropdowniconholder=[];
    dropdownpanel = []; %dropdown axes are placed in panels
    dropdownaxes=[];    
    parentobj=[];
    mouseovertext=[];
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
  end
  
  methods
    function g = myicon(name,parentobj,cdata,text,execute,type,group,cdataIndent,dropdowniconcell)%cdataClicked,cdataDisplay)
      
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
      
      %Set properties
      g.name=name;
      g.parentobj=parentobj;
      g.mouseovertext=text;
      g.cdata = cdata;
      g.cdataDisplay=cdata;
      g.execute=execute;
      
      if nargin<8 || isempty(cdataIndent)
        g.generateclickeddisabledandindent(cdata)
      else
        g.generateclickeddisabledandindent(cdata,cdataIndent)
      end
      
      %generate all stuff needed for dropdowniconholders
      if nargin==9
        g.dropdowniconcell = dropdowniconcell;
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
      end
      if g.type ~= 0 && g.type ~= 1 && g.type ~= 2 && g.type ~= 3
        g.enabled = 0;
        g.cdataDisplay=g.cdataDisabled;
      end
      
    end
    
    %--------------------------------------------------
    function generateclickeddisabledandindent(varargin)
      g = varargin{1};
      cdata = varargin{2};
      tmp = cdata;
     tmp(cdata==240) = 100;
      g.cdataClicked = uint8(tmp);
      cdataDisabled = rgb2gray(cdata)/5;
      cdataDisabled = cdataDisabled+(240-max(max(cdataDisabled)));%(cdataDisabled>240)=240;
      g.cdataDisabled=uint8(repmat(cdataDisabled,1,1,3));
      
      %Allows user to specify how indented button should look like
      if nargin<3
        tmp=g.cdata;
         tmp(g.cdata==240)=160;
        g.cdataIndent=uint8(tmp);
      else
        g.cdataIndent=varargin{3};
      end
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
      
      g=varargin{1};
      %First check so that button is enabled
      if g.enabled
        %Check that if the button we are over is the same as the the
        %clicked button, if so indent it. If not we want to unhighlight the previous button
        currenticon=g.parentobj.geticon;
        crit1 = isequal(currenticon,g.parentobj.clickedicon);
        crit2 = isequal(hittest(g.parentobj.fig),g.parentobj.imagehandle);
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
      g.enabled=1;
      g.cdataDisplay=g.cdata;
      g.parentobj.render
    end
    
    %-------------------------
    function disable(varargin)
      g=varargin{1};
      g.enabled=0;
      g.cdataDisplay=g.cdataDisabled;
      g.parentobj.render
    end
    
    %-------------------------
    function setIcon(varargin)
      g=varargin{1};
      g.cdata=varargin{2};
    end
  end
end