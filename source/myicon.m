%MYICON add icon to MYICONPLACEHOLDER
% Constructor takes (axes,image,text,exestring,dentstick (Says if the button stays indented afterpress),indentimage)
%Klas Berggren & Einar Heiberg

classdef myicon < handle %Inherits from handles to get persistent objects.

  properties (SetAccess = public)
    name='';
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
    function g = myicon(name,parentobj,cdata,text,execute,type,group,cdataIndent)%cdataClicked,cdataDisplay)
       
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
      if nargin<8
        g.generateclickeddisabledandindent(cdata)
      else
        g.generateclickeddisabledandindent(cdata,cdataIndent)
      end
    end
    
    function generateclickeddisabledandindent(varargin)
      g=varargin{1};
      cdata=varargin{2};
      tmp=cdata;
      tmp(cdata==240)=100;
      g.cdataClicked=uint8(tmp);
      cdataDisabled=rgb2gray(cdata)/5;
      cdataDisabled=cdataDisabled+(240-max(max(cdataDisabled)));%(cdataDisabled>240)=240;
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
    
    function highlight(varargin)  
      g=varargin{1};
       if g.enabled
         g.cdataDisplay=g.cdataClicked;
       end
      end
    
      function unhighlight(varargin)
        global DATA
      g=varargin{1};
      %global DATA
      %First check so that button is enabled
       if g.enabled
         %Check that if the button we are over is the same as the the
         %clicked button, if so indent it. If not we want to unhighlight the previous button
         currenticon=g.parentobj.geticon;
         %if isequal(g,g.parentobj.clickedicon) && isequal(hittest(DATA.fig),g.parentobj.imagehandle)
        if isequal(currenticon,g.parentobj.clickedicon) && isequal(hittest(DATA.fig),g.parentobj.imagehandle)
         switch g.type
            case 0
             g.cdataDisplay=g.cdata;
            case 1 %|| ~isempty(g.parentobj.geticon)
             g.cdataDisplay=g.cdataIndent;
             g.isindented=1;
            case 2
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
            if ~isempty(g.parentobj.timer)
              g.parentobj.displayinfo=0;
              stop(g.parentobj.timer);
              delete(g.parentobj.timer);
              g.parentobj.timer=[];
            end
              %profile report
         %if isequal(hittest(get(g.parentobj.axeshandle,'Parent')),g.parentobj.imagehandle)
           %segment(g.execute);
           feval(g.execute)%run(g.execute)
         %end
         else
           g.cdataDisplay=g.cdata;
         end
       end
      end
      
      function undent(varargin)
      g=varargin{1};
      g.isindented=0;
      if g.enabled
        g.cdataDisplay=g.cdata;
      else
        g.cdataDisplay=g.cdataDisabled;
      end
      end
      
      function enable(varargin)
      g=varargin{1};
      g.enabled=1;
      g.cdataDisplay=g.cdata;
      g.parentobj.render
      end
      
      function disable(varargin)
      g=varargin{1};
      g.enabled=0;
      g.cdataDisplay=g.cdataDisabled;
      g.parentobj.render
      end
    
    function setIcon(varargin)
      g=varargin{1};
      g.cdata=varargin{2};
    end
  end
end