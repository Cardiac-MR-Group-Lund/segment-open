%MYICONPLACEHOLDER Class to handle icons

%Klas Berggren & Einar Heiberg

classdef myiconplaceholder < handle %Inherits from handles to get persistent objects.

  properties
    cdata = [];
    numberoficons = 0;
    iconCell={};
    axeshandle = NaN;
    imagehandle = NaN;
    texthandle=[];
    displayinfo=0;
    timer=[];
    clickedicon=[];
    buttonisdown=0;
    isribbon=0;
    disablepad=0;
    configatclick=[];
    killtimer=[];
    pull2=[];
  end
  
  methods
    
    function g = myiconplaceholder(axes,isribbon,pull2)%,numberoficons)
      %Constructor takes icons 
      %Initialize properties
      g.axeshandle = axes;
      g.numberoficons =0;
      
      %Added possibility to design ribbon buttons
      if nargin == 2
        g.isribbon=isribbon;
      end
      
      %Set pull 2 always places axes to left = 1 or right = 2 
      if nargin < 3
        g.pull2=1;
      else
        g.pull2 =pull2;
      end
       
    end
    
    %function position = add(varargin)
    function add(varargin)
      %Adds an iconcell. Varagin contains [iconplaceholder,icon/cell of
      %icons] procedure is to overwrite current content of cell.
    g=varargin{1};
    
    icons=varargin{2};
    g.numberoficons=numel(icons);
    
    %Makes the axes usable
    %graphics
      g.imagehandle = image(g.cdata,'parent',g.axeshandle);
      
      %tooltip handle
      if g.pull2==2
        g.texthandle=text(1,1,'','Parent',g.axeshandle,'Background','white','VerticalAlignment','baseline','HorizontalAlignment','Right','HitTest','off');
      else
        g.texthandle=text(1,1,'','Parent',g.axeshandle,'Background','white','VerticalAlignment','baseline','HorizontalAlignment','Left','HitTest','off');
      end
      set(g.texthandle,'visible','off')
     
      %Adjust axis
      axis(g.axeshandle,'image','off');
      axis(g.axeshandle,'ij') 
      set(g.axeshandle,'DataAspectRatio',[1,1,1]);
      
      %Set clickcallback
       set(g.imagehandle,'ButtonDownFcn',@g.click);
      
      
    if iscell(icons)
      %try to add entire cell position is linear indices
      g.iconCell=icons;
    else
      ME = MException('FaultyInput','Please give icons in cell of your desired format of the toolbar.');
      throw(ME)
    end
    %render image
    g.render
    end
    
    function render(varargin)
      g=varargin{1};
      cdataCell=cell(size(g.iconCell));
      for i =1:g.numberoficons
        icon=g.iconCell{i};
        cdataCell{i}=icon.cdataDisplay;
      end
      cdata=cell2mat(cdataCell);
      
      %idea 2
      %get axes height then interpolate so that rows= axes height and
      %length is interpolated so that we have the same aspect ratio.
      set(g.axeshandle,'Units','pixels');
      axpos = get(g.axeshandle,'Position'); 
      height=axpos(4);
      cdata=imresize(cdata,height/size(cdata,1));%[height,height*g.numberoficons,3]);
      
      %Here we add a line that fills out the remaining right part of the
      %axes if we are dealing with a ribbon interface
      pad=[];
      if g.isribbon
        ribbonline=240*ones(size(cdata,1),1);
        ribbonline(end-1:end)=0;
       pad=repmat(imresize(ribbonline,height/size(cdata,1)),1,int64(axpos(3)-size(cdata,2)+100),3);
       %We get some weird interpolation so find make image binary by
       %setting everything but the zeros to 0
       pad(pad~=0)=240;
       
       if g.disablepad
         %We know that pad is binary colorwise so lets pick the darkest for
         %the line and brightest for rest
         pad(pad==0)=min(min(min(cdata)));
%         pad=rgb2gray(pad)/4;
%         pad=pad+(240-max(max(pad)));%(cdataDisabled>240)=240;
%         pad=uint8(repmat(pad,1,1,3));
       end
        %pad=zeros(int64([size(cdata,1),axpos(3)-size(cdata,2)+100,3]));%the plus 100 is so that we always have slightly more to avoid edge at right end.   
      end
      
      g.cdata=cdata;
      set(g.axeshandle,'Units','normalized');
      set(g.axeshandle,'Visible','off');
      set(g.imagehandle,'cdata',[g.cdata,pad]);
      
      %switch that renders images.
      axis(g.axeshandle,'image','on');
      axis(g.axeshandle,'image','off');
      
      %Assert position oficon placeholders.
      if g.pull2==1
        pos=plotboxpos(g.axeshandle);
        currentpos=get(g.axeshandle,'position');
        set(g.axeshandle,'position',currentpos-[pos(1),0,0,0]);
      else
        pos=plotboxpos(g.axeshandle);
        currentpos=get(g.axeshandle,'position');
        currentpos(1)=currentpos(1)+1-(pos(1)+pos(3));
        set(g.axeshandle,'position',currentpos);
      end
      
      if ~isempty(g.timer)
        g.displayinfo=0;
        stop(g.timer);
        delete(g.timer);
        g.timer=[];
      end
    end
    
    function iconsandstates=iconson(varargin)
      %Returns all type 2 icons and there states.
      g=varargin{1};
      
      iconcell=g.iconCell;
      iconsandstates=cell(g.numberoficons,2);
      for i=1:g.numberoficons
        if iconcell{i}.type==2
          iconsandstates{i,1}=iconcell{i}.name;
          iconsandstates{i,2}=iconcell{i}.isindented;
        end
      end
      
      %clean output
      emptyCells = cellfun('isempty', iconsandstates); 
      iconsandstates(all(emptyCells,2),:) = [];
    end
    
    function overicon=geticon(varargin)
      g=varargin{1};
      try
      %find which icon we are over then return it.
      clicked=get(g.axeshandle, 'CurrentPoint');
      x_over=clicked(1,1);
      y_over=clicked(1,2);
      [rows,cols]=size(g.iconCell);
      
      %No longer true as we rescale the image constantly
      %[yres,xres,~] = size(g.iconCell{1}.cdata);
      yres = size(g.cdata,1);
      xres = size(g.cdata,2);
      %xres = size(g.cdata,2)/g.numberoficons;%Not true cant handle variable size buttons
      
      %get normalized length of icons in iconcell inorder to set up xres
      totXsize=0;
      totYsize=0;
      for i=1:cols
        totXsize=totXsize+size(g.iconCell{1,i}.cdata,2);
      end
      
      for i=1:rows
        totYsize=totYsize+size(g.iconCell{i,1}.cdata,1);
      end
      
      X_norm=zeros(size(g.iconCell));
      Y_norm=zeros(size(g.iconCell));
      
      
      for i=1:numel(g.iconCell)
        X_norm(i)=size(g.iconCell{i}.cdata,2)/totXsize;
        Y_norm(i)=size(g.iconCell{i}.cdata,1)/totYsize;
      end
      
      X_norm=reshape(cumsum(X_norm,2),1,numel(g.iconCell));
      Y_norm=reshape(cumsum(Y_norm,1),1,numel(g.iconCell));
      
      %X=0:xres:xres*cols;
      X=[0,X_norm.*xres];
      Y=[0,Y_norm.*yres];
      %Y=0:yres:yres*rows;
      x_ind=find(x_over>X,1,'last');
      y_ind=find(y_over>Y,1,'last');
      %x_ind=find(x_over>=X,1,'last');
      %y_ind=find(y_over>=Y,1,'last');
      overicon=g.iconCell{y_ind,x_ind};
      catch
        %there is a possibility that when we query for location in this
        %function we nolonger are over a icon
        overicon=[];
      end
    end
    
    function indent(varargin)
      
      g=varargin{1};
      name=varargin{2};
      runicon=varargin{3};
      icon=[];
      for i=1:numel(g.iconCell)
          if strcmp(g.iconCell{i}.name,name)
            icon=g.iconCell{i};
            break;
          end
      end
      
      if isempty(icon)
      return;
      end
      %if type toggle tool then traverse list again and make sure all togglers in same group are undented
      if icon.type==1
        for i=1:g.numberoficons
          if g.iconCell{i}.type==1 && g.iconCell{i}.group==icon.group
            g.iconCell{i}.cdataDisplay=g.iconCell{i}.cdata;
            g.iconCell{i}.isindented=0;
          end
        end
      end
      icon.cdataDisplay=icon.cdataIndent;
      icon.isindented=1;
      if runicon
        feval(icon.execute)
      end
      g.render;
    end
    
    function click(varargin)
      global DATA
      %profile on;
      g=varargin{1};
      g.buttonisdown=1;
      
      %Save the current click configuration if slide off icon axes
      g.configatclick=g.findindented;
      
      if isequal(hittest(DATA.fig),g.imagehandle) 
        set(DATA.fig,'WindowButtonUpFcn',@g.upclick);
        g.clickedicon=g.geticon;
         if ~isempty(g.clickedicon)
        g.clickedicon.highlight
        set(g.texthandle,'visible','off')
        g.displayinfo=0;
        if ~isempty(g.timer)
          stop(g.timer);
          delete(g.timer);
          g.timer=[];
        end
      g.render
      %pause(0.1)
       %inorder to trigger text again
       %g.motion;
         end
      end
    end
    
    function indented=findindented(varargin)
      %returns toggle buttons (icon.type=1) buttons that are indented
      g=varargin{1};
      indented=[];
       for i=1:(g.numberoficons)
          if g.iconCell{i}.isindented %&& g.iconCell{i}.type==1
            indented=[indented,i];
            %return;
          end
       end
    end
    
    function upclick(varargin)
      global DATA
      g=varargin{1};
      
      g.buttonisdown=0;
      if isequal(hittest(DATA.fig),g.imagehandle) && ~isempty(g.clickedicon) %&& isequal(g.clickedicon
        indented=g.findindented;
        for k= indented
          if ~isempty(g.clickedicon) && g.iconCell{k}.group==g.clickedicon.group && g.iconCell{k}.type==1 && g.clickedicon.type==1
            g.iconCell{k}.undent;
          end
        end
        g.clickedicon.unhighlight
      else
        g.notover;
      end
    end

    function motion(varargin)
      global DATA
      g=varargin{1};
      
%      if isequal(hittest(get(g.axeshandle,'Parent')),g.imagehandle)
        current_icon=g.geticon;
     
        
        if ~isempty(current_icon)
        %If the icon we are over currently isnt the one we clicked
        %originally we want to unhighlight the prior icon and highlight
        %the new one
        
        %No sliding over icons allowed
        if ~isequal(current_icon,g.clickedicon) && ~isempty(g.clickedicon)...
            && g.buttonisdown && ~isempty(current_icon)
          g.notover
        end
        
        %Start timer so that text doesnt display immediately
        if isempty(g.timer) && ~g.displayinfo
            g.timer=timer('ExecutionMode','singleShot',...
              'TimerFcn',@g.displayinfoYN, ...
              'StartDelay',1);
            start(g.timer);
        end 
        
        %if timer has triggered and we are over the toolbar show text
        if g.displayinfo
          set(g.texthandle,'visible','on')
          xy = get(g.axeshandle, 'CurrentPoint');
          x=xy(1,1);
          y=xy(1,2);
          overicon=g.geticon;
          set(g.texthandle,'Position',[x,y,0],'String',translation.dictionary(overicon.mouseovertext));%text(x,y,'blablabla','Parent',g.axeshandle);
        end
        else 
          g.notover
        end
    end
    
    function killtext(varargin)
  global DATA
      g=varargin{1};
   % Check mouse location on screen
   mouseloc=get(0, 'PointerLocation');
   guiloc=get(DATA.fig, 'Position');
  state= inpolygon(mouseloc(1),mouseloc(2),[guiloc(1),guiloc(1)+guiloc(3)],[guiloc(1),guiloc(2)+guiloc(4)]) ;
   
      %if nolonger over parentaxes kill text
      if ~state%isempty(overobj(DATA.fig))%~isequal(hittest(DATA.fig),g.imagehandle) %&& overobj(DATA.fig)
        %if we no longer are over the toolbar, kill timer.
        set(g.texthandle,'visible','off')
        g.displayinfo=0;
        if ~isempty(g.timer)
          stop(g.timer)
          delete(g.timer)
          g.timer=[];
          
          %remove killtimer also
          stop(g.killtimer)
          delete(g.killtimer)
          g.killtimer=[];
        end
      else
        stop(g.killtimer)
        delete(g.killtimer)
        g.killtimer=[];
        g.killtimer=timer('ExecutionMode','singleShot',...
          'TimerFcn',@g.killtext, ...
          'StartDelay',1);
        start(g.killtimer);
      end
      
      
    end
    
    %      function notover(varargin)
%       global DATA
%       g=varargin{1};
%       
%       if ~isempty(g.clickedicon) && g.buttonisdown
%         g.clickedicon.unhighlight;
%         g.clickedicon=[];
%         g.render
%       end
%       
%       %if a user clicks and holds and drags of the iconholder and no button
%       %is clicked we clear all indentations and remove buttonupfcn.
%       if isempty(g.clickedicon) && ~isempty(g.findindented)  && g.buttonisdown
%         if g.axeshandle==DATA.Handles.configiconholder.axeshandle;
%         for i=1:g.findindented
%           
%         %undent all but current tool 
%           if g.iconCell{i}.type==1 && ~strcmp(DATA.CurrentTool,g.iconCell{i}.name)
%             g.iconCell{i}.undent;
%           else
%             g.iconCell{i}.isindented=1;
%             g.iconCell{i}.cdataDisplay=g.iconCell{i}.cdataIndent;
%           end
%         end
%         else
%               %First the panel group is found
%     str=sprintf('%d%d',DATA.ViewMatrix(1),DATA.ViewMatrix(2));
%     switch str
%       case '11'
%       name='panel1';
%       case '12'
%       name='panel2';
%       case '21'
%       name='panel2x1';
%       case '13'
%       name='panel3';
%       case '31'
%       name='panel3x1';
%       case '22'
%         name='panel4';
%       case '23'
%         name='panel6';
%     end
%     
%     iconcell=g.iconCell;
%     for i= 1:g.numberoficons  
%       if strcmp(iconcell{i}.name,name)
%         ind=i;
%       end
%     end
%     
%     iconcell{ind}.cdataDisplay=iconcell{ind}.cdataIndent;
%     iconcell{ind}.isindented=1;
%     
%     %Then view mode
%     switch DATA.ViewPanelsType{DATA.CurrentPanel}
%       case 'montage'
%         name='viewall';
%       case 'one'
%         name='viewone';
%     end
%     
%     for i= 1:g.numberoficons  
%       if strcmp(iconcell{i}.name,name)
%         ind=i;
%       end
%     end
%     
%     iconcell{ind}.cdataDisplay=iconcell{ind}.cdataIndent;
%     iconcell{ind}.isindented=1;
%         end
%         
%         set(DATA.fig,'WindowButtonUpFcn',[]);
%         g.render
%       end
%       
%       %if we no longer are over the toolbar, kill timer.
%       set(g.texthandle,'visible','off')
%       g.displayinfo=0;
%       if ~isempty(g.timer)
%         stop(g.timer)
%         delete(g.timer)
%         g.timer=[];
%       end
%     end
%     
%     function displayinfoYN(varargin)
%       %Should we display info now that timer has triggered
%       global DATA
%       g=varargin{1};
%       if  isequal(hittest(DATA.fig),g.imagehandle) && ~g.buttonisdown && ~g.isribbon
%         g.displayinfo=1;
%         set(g.texthandle,'visible','on')
%         overicon=g.geticon;
%         xy = get(g.axeshandle, 'CurrentPoint');
%         x=xy(1,1);
%         y=xy(1,2);
%         set(g.texthandle,'Position',[x,y,0],'String',overicon.mouseovertext);%text(x,y,'blablabla','Parent',g.axeshandle);
%       else
%         g.displayinfo=0;
%       end
%     end
%     
    function notover(varargin)
      % DATA
      g=varargin{1};
      
      if ~isempty(g.clickedicon) && g.buttonisdown
        g.clickedicon.undent;%unhighlight;
        g.clickedicon=[];
        for i=g.configatclick
          g.iconCell{i}.cdataDisplay=g.iconCell{i}.cdataIndent;
          g.iconCell{i}.isindented=1;
        end
          g.render
      end
      
      %if we no longer are over the toolbar, kill timer.
      set(g.texthandle,'visible','off')
      g.displayinfo=0;
      if ~isempty(g.timer)
        stop(g.timer)
        delete(g.timer)
        g.timer=[];
      end
    end
    
    function textoff(g)
       set(g.texthandle,'visible','off')
      g.displayinfo=0;
      if ~isempty(g.timer)
        stop(g.timer)
        delete(g.timer)
        g.timer=[];
      end
    end
    
    function displayinfoYN(varargin)
      %Should we display info now that timer has triggered
      global DATA
      g=varargin{1};
      if  isequal(hittest(DATA.fig),g.imagehandle) && ~g.buttonisdown && ~g.isribbon
        g.displayinfo=1;
        set(g.texthandle,'visible','on')
        overicon=g.geticon;
        xy = get(g.axeshandle, 'CurrentPoint');
        x=xy(1,1);
        y=xy(1,2);
        set(g.texthandle,'Position',[x,y,0],'String',translation.dictionary(overicon.mouseovertext));%text(x,y,'blablabla','Parent',g.axeshandle);
        
        %Start kill text timer that checks for a while if we are over the axis if not kill text
        %if isempty(g.killtimer)
       
          g.killtimer=timer('ExecutionMode','singleShot',...
            'TimerFcn',@g.killtext, ...
            'StartDelay',0.2);
          start(g.killtimer);
       %end
      else
        g.displayinfo=0;
      end
    end
    
    function kill
      %Function which empties the iconcell and renders
    end
    
  end
end
  