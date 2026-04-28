classdef myiconplaceholder < handle
  %Inherits from handles to get persistent objects.
  %MYICONPLACEHOLDER Class to handle icons

  %Klas Berggren & Einar Heiberg

  properties
    cdata = [];
    numberoficons = 0;
    iconCell={};
    axeshandle = NaN;
    imagehandle = NaN;
    dropdowniconholders=[];
    texthandle=[];
    displayinfo=0;
    texttimer=[];
    clickedicon=[];
    buttonisdown=0;
    isribbon=0;
    disablepad=0;
    configatclick=[];
    killtimer=[];
    pull2=[];
    fig=[];
    iconarray = zeros(128,128,3,'uint8');
    themebackground = false; %use Segment's background colour
  end

  methods
    %--------------------------------
    function g = myiconplaceholder(axes,isribbon,pull2,fig,dobackground)%,numberoficons)
      %--------------------------------
      %Constructor takes icons
      %Initialize properties
      g.axeshandle = axes;
      g.numberoficons =0;

      %Added possibility to design ribbon buttons
      if nargin == 2
        g.isribbon = isribbon;
      end

      %Set pull 2 always places axes to left = 1 or right = 2
      if nargin < 3
        g.pull2 = 1;
      else
        g.pull2 = pull2;
      end

      if nargin<4
        global DATA
        g.fig = DATA.fig;
      else
        g.fig = fig;
      end

      if nargin < 5
        dobackground = true;
      end

      if dobackground
        g.themebackground = true;
        g.applybackground;
      end
    end

    %---------------------------------------------
    function applybackground(varargin)
      %---------------------------------------------
      %Apply theme colour to the panels containing the iconholder
      g = varargin{1};
      global DATA
    
      panel = get(g.axeshandle,'Parent');
      children = get(panel,'Children');
      allcomponents = [panel;children];
      for loop = 1:numel(allcomponents)
        if isa(allcomponents(loop),'matlab.ui.container.Panel') || ...
            isa(allcomponents(loop),'matlab.ui.container.ButtonGroup')
          set(allcomponents(loop),'BackgroundColor',DATA.GUISettings.BackgroundColor);
        end
        if isa(allcomponents(loop),'matlab.ui.container.ButtonGroup')
          buttongroupcomponents = get(allcomponents(loop),'Children');
          for button = 1:numel(buttongroupcomponents)
            set(buttongroupcomponents(button),'BackgroundColor',DATA.GUISettings.BackgroundColor,...
              'ForegroundColor',DATA.GUISettings.ForegroundColor);
          end
        end
      end
    end

    %function position = add(varargin)
    %--------------------------------
    function add(varargin)
      %--------------------------------
      %Adds an iconcell. Varagin contains [iconplaceholder,icon/cell of
      %icons] procedure is to overwrite current content of cell.
      g=varargin{1};

      icons=varargin{2};
      g.numberoficons=numel(icons);

      %clear any dropdown if adding new icons to placeholder
      for iconholder = g.dropdowniconholders
        %       if ishandle(iconholder)
        delete(iconholder.axeshandle.Parent)
        delete(iconholder.axeshandle)
        iconholder.axeshandle = NaN;
        delete(iconholder);
        %       end
      end
      g.dropdowniconholders=[];

      %if icons are added to a new placeholder this will be their new parent.
      for i = 1:g.numberoficons
        if not(isempty(icons{i}))
          icons{i}.parentobj=g;

          %add iconholder to axes and render
          if not(isempty(icons{i}.type))&& icons{i}.type == 3 % this is the dropdown type
            %         if ~ishandle(icons{i}.dropdowniconholder)
            p = uipanel('parent',g.fig,'units','normalized','position',[0 0 1 1],'BorderType','none','Visible','off','Clipping','off');
            ax = axes('parent',p,'units','normalized','position',[0 0 1 1],'Clipping','off');
            icons{i}.dropdownpanel = p;
            icons{i}.dropdownaxes = ax;
            icons{i}.dropdowniconholder = myiconplaceholder(ax,0,3);
            add(icons{i}.dropdowniconholder,icons{i}.dropdowniconcell)
            %         end

            if isempty(g.dropdowniconholders)
              g.dropdowniconholders = icons{i}.dropdowniconholder;
            else
              g.dropdowniconholders(end+1) = icons{i}.dropdowniconholder;
            end
          end
        end

      end
      %Adjust axis
      %       axis(g.axeshandle,'image','off');
      %       axis(g.axeshandle,'ij')
      set(g.axeshandle,'DataAspectRatio',[1,1,1],'fontsmoothing','off','Clipping','off');
      %       set(g.axeshandle,);
      %       set(g.axeshandle,);
      %set(g.axeshandle,'graphicssmoothing','off');

      if iscell(icons)
        %try to add entire cell position is linear indices
        g.iconCell=icons;
      else
        ME = MException('FaultyInput','Please give icons in cell of your desired format of the toolbar.');
        throw(ME)
      end
      %render image
      g.render
      set(g.imagehandle,'ButtonDownFcn',@g.click);
    end

    %--------------------------------
    function enable(varargin)
      %--------------------------------
      %Enable icons in the iconholder. Disable is not FDA cleared.
      g=varargin{1};
      if nargin > 1
        names = varargin{2};
      end
      if nargin > 2
        runFDAVersion = varargin{3}; %check if the GUI is run in FDA mode
      else
        runFDAVersion = false; %FDA mode is explicitly given from Segment, if implicit then not FDA mode
      end

      dorendering = false; %do not render on icon level, but on iconholder level

      if isempty(names)
        for i = 1:g.numberoficons
          currenticon = g.iconCell{i};
          if (runFDAVersion && currenticon.isFDACleared) || ~runFDAVersion
            if ~currenticon.enabled
              currenticon.enable(runFDAVersion,dorendering);
            end
          else
            if currenticon.enabled
              currenticon.disable(dorendering);
            end
          end
        end
      else
        for i = 1:g.numberoficons
          if any(strcmp(g.iconCell{i}.name,names)) && ~g.iconCell{i}.enabled
            g.iconCell{i}.enable(runFDAVersion,dorendering);
          end
        end
      end
      g.render
    end

    %--------------------------------
    function disable(varargin)
      %--------------------------------
      g = varargin{1};

      if nargin==1
        for i=1:g.numberoficons
          g.iconCell{i}.enabled=0;
          if not(strcmp(g.iconCell{i}.name,'verticallinewide')) && not(strcmp(g.iconCell{i}.name,'verticalline'))
            g.iconCell{i}.cdataDisplay = g.iconCell{i}.cdataDisabled;
          end
        end
      else
        names = varargin{2};
        for i=1:g.numberoficons
          if any(strcmp(g.iconCell{i}.name,names))
            %g.iconCell{i}.disable
            g.iconCell{i}.enabled=0;
            if not(strcmp(g.iconCell{i}.name,'verticallinewide')) && not(strcmp(g.iconCell{i}.name,'verticalline'))
              g.iconCell{i}.cdataDisplay = g.iconCell{i}.cdataDisabled;
            end
          end
        end
      end
      g.render
    end

    %--------------------------------
    function render(varargin)
      %--------------------------------
      g = varargin{1};
      global DATA

      cdataCell = cell(size(g.iconCell));
      for i =1:g.numberoficons
        icon = g.iconCell{i};
        cdataCell{i} = icon.cdataDisplay;
      end

      numrows = size(cdataCell,1);
      numcols = size(cdataCell,2);
      if numrows > 1 ...&& numcols > 1
        % find maximum length of an icon in the column
        maxlencol = zeros(1,numcols);
        for col = 1:numcols
          for row = 1:numrows
            maxlencol(1,col) = max(maxlencol(1,col),size(cdataCell{row,col},2));
          end
        end
        for col = 1:numcols
          for row = 1:numrows
            maxlencol = max(maxlencol,size(cdataCell{row,col},2));
            actysize = size(cdataCell{row,col},2);
            if actysize < maxlencol(col) && actysize ~= 14
              % padding with zeros
              cdataCell{row,col} = padarray(cdataCell{row,col},[0,maxlencol(col)-actysize],240,'post');
            end
          end
        end
        %           tempdata = cell(numrows,1);
        %           maxysize = 0;
        %           for row = 1:numrows
        %             tempdata{row,1} = cell2mat(cdataCell(row,:));
        %             maxysize = max(maxysize,size(tempdata{row,1},2));
        %           end
        %           for row = 1:numrows
        %             actysize = size(tempdata{row,1},2);
        %             if actysize < maxysize
        %               % padding with zeros
        %               tempdata{row,1} = padarray(tempdata{row,1},[0,maxysize-actysize],240,'post');
        %             end
        %           end
        cdata = cell2mat(cdataCell);
      else
        cdata = cell2mat(cdataCell);
      end

      if isempty(cdata)
          return
      end


      %idea 2
      %get axes height then interpolate so that rows= axes height and
      %length is interpolated so that we have the same aspect ratio.
      set(g.axeshandle,'Units','pixels');
      axpos = get(g.axeshandle,'Position');
      height = axpos(4);
      cdata = imresize(cdata,height/size(cdata,1));%imresize(cdata,height/size(cdata,1));%[height,height*g.numberoficons,3]);

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

         tmp = pad;
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
      %Recombine the separate channels into a single RGB image
      pad = cat(3, redchannel, greenchannel, bluechannel);


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
      set(g.axeshandle,'Units','normalized','Visible','off');

      %graphics
      if ~ishandle(g.imagehandle)
        g.imagehandle = imagesc(g.cdata,'parent',g.axeshandle); %If handle lost or not yet created recreate it and the texthandle
        axis(g.axeshandle,'image','off');
        %tooltip handle
        if g.pull2==2
          g.texthandle=text(1,1,'','Parent',g.axeshandle,'Background','white','VerticalAlignment','bottom','HorizontalAlignment','Right','HitTest','off');
        else
          g.texthandle=text(1,1,'','Parent',g.axeshandle,'Background','white','VerticalAlignment','bottom','HorizontalAlignment','Left','HitTest','off');
        end
        set(g.texthandle,'visible','off')
      else
        set(g.imagehandle,'cdata',[g.cdata,pad]);
      end

      %switch that renders images.
      %       axis(g.axeshandle,'image','on');
      %       axis(g.axeshandle,'image','off');

      %Assert position oficon placeholders.
      switch g.pull2
        case 1
          pos=plotboxpos(g.axeshandle);
          currentpos=get(g.axeshandle,'position');
          set(g.axeshandle,'position',currentpos-[pos(1),0,0,0]);
        case 2
          pos=plotboxpos(g.axeshandle);
          currentpos=get(g.axeshandle,'position');
          currentpos(1)=currentpos(1)+1-(pos(1)+pos(3));
          set(g.axeshandle,'position',currentpos);
        case 3
          %Stay where you are
      end

      %if any dropdown iconholders there adjust their position
      for i = 1:numel(g.iconCell)
        if g.iconCell{i}.type==3 %&& g.iconCell{i}.isindented %this is the dropdown type
          %The below procedure reboots the dropdowns
          g.iconCell{i}.placedropdown%(0)
          g.iconCell{i}.showdropdown
          %g.iconCell{i}.setdropdown(1)
          %elseif g.iconCell{i}.type==3 && ~g.iconCell{i}.isindented && ~isempty(g.iconCell{i}.dropdownaxes)
          %g.iconCell{i}.setdropdown(0)
        end
      end

      if ~isempty(g.texttimer)
        g.displayinfo=0;
        stop(g.texttimer);
        delete(g.texttimer);
        g.texttimer=[];
      end
    end

    %--------------------------------
    function iconsandstates=iconson(varargin)
      %--------------------------------
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

    %--------------------------------
    function overicon = geticon(varargin)
      %--------------------------------
      g = varargin{1};
      try
        %find which icon we are over then return it.
        clicked = get(g.axeshandle, 'CurrentPoint');
        x_over = clicked(1,1);
        y_over = clicked(1,2);
        [rows,cols] = size(g.iconCell);

        %No longer true as we rescale the image constantly
        %[yres,xres,~] = size(g.iconCell{1}.cdata);
        yres = size(g.cdata,1);
        xres = size(g.cdata,2);
        %xres = size(g.cdata,2)/g.numberoficons;%Not true cant handle variable size buttons

        %get normalized length of icons in iconcell inorder to set up xres
        totXsize = 0;
        totYsize = 0;
        for i=1:cols
          totXsize = totXsize+size(g.iconCell{1,i}.cdata,2);
        end

        for i = 1:rows
          totYsize = totYsize + size(g.iconCell{i,1}.cdata,1);
        end

        X_norm = zeros(size(g.iconCell));
        Y_norm = zeros(size(g.iconCell));

        for i = 1:numel(g.iconCell)
          X_norm(i) = size(g.iconCell{i}.cdata,2)/totXsize;
          Y_norm(i) = size(g.iconCell{i}.cdata,1)/totYsize;
        end

        if rows > 1 && cols >1
          X_norm = reshape(cumsum(X_norm(1,:),2),1,cols);
          Y_norm = reshape(cumsum(Y_norm(:,1),1),1,rows);
        else
          X_norm = reshape(cumsum(X_norm,2),1,numel(g.iconCell));
          Y_norm = reshape(cumsum(Y_norm,1),1,numel(g.iconCell));
        end
        X = [0,X_norm.*xres];
        Y = [0,Y_norm.*yres];

        x_ind = find(x_over>X,1,'last');
        y_ind = find(y_over>Y,1,'last');
        overicon = g.iconCell{y_ind,x_ind};
      catch
        %there is a possibility that when we query for location in this
        %function we nolonger are over a icon
        overicon=[];
      end
    end

    %--------------------------------
    function undent(varargin)
      %--------------------------------
      %g.undent(name,runicon)
      
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
      icon.cdataDisplay=icon.cdata;
      icon.isindented=0;
      g.render;
      if runicon
        feval(icon.execute);
      end
    end

    %--------------------------------
    function indent(varargin)
      %--------------------------------
      %g.indent(name,runicon)
      g = varargin{1};
      name = varargin{2};
      runicon = varargin{3};

      icon=[];
      %all the iconholders including dropdowniconholders
      iconholders = [g,g.dropdowniconholders];

      %extract the iconCells from the above iconholders then do name check
      shouldskip = false;
      parenticonname = [];
      for ich = 1:numel(iconholders)
        g = iconholders(1,ich);
        for i = 1:numel(g.iconCell)
          if strcmp(g.iconCell{i}.name,name)
            icon = g.iconCell{i};
            if ich > 1
              % index larger than 1 indicates that the icon is a dropdown
              % find tha corresponding index of parent icon
              parentind = cellfun(@(x) isequal(x.dropdowniconholder, icon.parentobj), iconholders(1,1).iconCell);
              parenticonname = iconholders(1,1).iconCell{parentind}.name;
            end
            shouldskip = true;
            break;
          end
        end
        if shouldskip
          break;% go further as sson as the name was found
        end
      end

      if isempty(icon) || ~icon.enabled
        return
      end

      %if type toggle tool then traverse list again and make sure all togglers in same group are undented
      if icon.type==1
        for i = 1:g.numberoficons
          if g.iconCell{i}.type==1 && g.iconCell{i}.group==icon.group && g.iconCell{i}.enabled==1
            g.iconCell{i}.cdataDisplay=g.iconCell{i}.cdata;
            g.iconCell{i}.isindented=0;
          end
        end
        if ~isempty(parenticonname)
          iconholders(1,1).indent(parenticonname,1)
        end
      end

      %If icon can not be indented (for example "zoom"), do nothing. 
      if ~icon.type==0        
        icon.cdataDisplay=icon.cdataIndent;
        icon.isindented=1;
        g.render;
      end

      if runicon
        feval(icon.execute);
      end
    end

    %--------------------------------
    function click(varargin)
      %--------------------------------
      g = varargin{1};
      didreset = g.resetbuttonupfcn;
      if didreset
        return
      end

      g.buttonisdown = 1;
      %Save the current click configuration if slide off icon axes
      g.configatclick = g.findindented;

      if isequal(hittest(g.fig),g.imagehandle)
        set(g.fig,'WindowButtonUpFcn',@g.upclick);
        g.clickedicon = g.geticon;
        if ~isempty(g.clickedicon)
          iconcells = g.clickedicon.parentobj.iconCell;
          if ~contains(g.clickedicon.name, 'verticalline')
            g.clickedicon.highlight
          end

          set(g.texthandle,'visible','off')
          g.displayinfo=0;
          if ~isempty(g.texttimer)
            stop(g.texttimer);
            delete(g.texttimer);
            g.texttimer=[];
          end
          g.render
          %pause(0.1)
          %inorder to trigger text again
          %g.motion;
        end
      end
    end

    %-------------------------------------
    function isreset = resetbuttonupfcn(g)
      %-----------------------------------
      % function to reset WindowButtonUpFcn if it fits criteria
      isreset = false;
      icon = g.clickedicon;

      if isvalid(g.fig) && ~isempty(g.fig.WindowButtonUpFcn)
        if ~isempty(icon) && icon.type ~= 2

          if isa(g.fig.WindowButtonUpFcn,"function_handle")
            fcn = func2str(g.fig.WindowButtonUpFcn);
          else
            fcn = g.fig.WindowButtonUpFcn;
          end

          if contains(fcn,'g.upclick')
            isreset = true;
            set(g.fig,'WindowButtonUpFcn','');
          end
        end
      end
    end

    %--------------------------------
    function indented=findindented(varargin)
      %--------------------------------
      %returns toggle buttons (icon.type=1) buttons that are indented
      g=varargin{1};
      if nargin ==1
        indented=[];
        for i=1:(g.numberoficons)
          if g.iconCell{i}.isindented %&& g.iconCell{i}.type==1
            indented=[indented,i];
            %return;
          end
        end
      else
        name=varargin{2};
        indented = 0;
        for i=1:(g.numberoficons)
          if g.iconCell{i}.isindented && strcmp(g.iconCell{i}.name,name)
            indented=1;
            return;
          end
        end
      end

    end

    %--------------------------------
    function [x,y]=buttoncorner(varargin)
      %--------------------------------
      g=varargin{1};
      name=varargin{2};
      pos =plotboxpos(g.axeshandle);%
      %pos = get(g.axeshandle,'position');

      if strcmp(g.axeshandle.Parent.Type,'figure')
        panelpos = ones(1,4);
        panelpos(1:2)=0;
      else
        panelpos = g.axeshandle.Parent.Position;
      end

      cdataCell=cell(size(g.iconCell));
      for i =1:g.numberoficons
        cdataCell{i}=g.iconCell{i}.cdataDisplay;
      end
      sz =  size(cell2mat(cdataCell));
      y = zeros(1,size(g.iconCell,1)+1);
      x = zeros(1,size(g.iconCell,2)+1);

      for i= 1:size(g.iconCell,1)
        for j = 1:size(g.iconCell,2)
          if strcmp(g.iconCell{i,j}.name,name)
            row = i;
            col = j;
          end
          x(j+1)=size(g.iconCell{i,j}.cdata,2);
        end
        y(i+1)=size(g.iconCell{i,j}.cdata,1);
      end
      x = cumsum(x(1:end-1));
      y = cumsum(y(1:end-1));

      %normalized coordinate system
      wp= sz(2);
      wn=pos(3)*panelpos(3);
      hp= sz(1);
      hn=pos(4)*panelpos(4);
      x = pos(1)+panelpos(1)+x(col)/wp*wn;
      y = pos(2)+panelpos(2)+y(row)/hp*hn;
    end

    %--------------------------------
    function upclick(varargin)
      %--------------------------------
      global DATA
      g = varargin{1};
      
      if DATA.Testing
        crit1 = true;
      else
        crit1 = isequal(hittest(g.fig),g.imagehandle);
      end
      crit2 = ~isempty(g.clickedicon);

      %      fig=get(g.axeshandle,'parent');
      g.buttonisdown=0;
      if crit1 && crit2
        indented=g.findindented;
        for k = indented
          if ~isempty(g.clickedicon) && g.iconCell{k}.group==g.clickedicon.group && g.iconCell{k}.type==1 && g.clickedicon.type==1
            if g.clickedicon.enabled
              g.iconCell{k}.undent;
            end
          end
        end
        iconcells = g.clickedicon.parentobj.iconCell;
        g.clickedicon.unhighlight
        g.resetbuttonupfcn();
      else
        g.notover;
      end

      %if any dropdown iconholders there show them
      %       for i = 1:numel(g.iconCell)
      %         if g.iconCell{i}.type==3 %&& g.iconCell{i}.isindented %this is the dropdown type
      %           g.iconCell{i}.showdropdown%(0)
      %         end
      %       end

    end

    %--------------------------------
    function motion(varargin)
      %--------------------------------
      %global DATA
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

        %Start texttimer so that text doesnt display immediately
        if isempty(g.texttimer) && ~g.displayinfo
          g.texttimer=timer('ExecutionMode','singleShot',...
            'TimerFcn',@g.displayinfoYN, ...
            'Name','TextTimer',...
            'StartDelay',1);
          start(g.texttimer);
        end

        %         %if timer has triggered and we are over the toolbar show text
        if g.displayinfo

          xy = get(g.axeshandle, 'CurrentPoint');
          x=xy(1,1);
          y=xy(1,2);

          xl = xlim(g.axeshandle);
          yl = ylim(g.axeshandle);

          if y/yl(2)<0.4
            verticalalignment = 'top';
          else
            verticalalignment = 'bottom';
          end

          if g.pull2==2 || x/xl(2)>0.8%pull2
            horizontalalignment = 'right';
            offset = -10;
          else
            horizontalalignment = 'left';
            offset = 10;
          end
          txt = gettoolinfo(current_icon.name);
          if isempty(txt)
            txt = translation.dictionary(current_icon.mouseovertext);
          end
          set(g.texthandle,'Position',[x+offset,y,0],'String',sprintf(txt),...
            'verticalalignment',verticalalignment,'horizontalalignment',...
            horizontalalignment,'Clipping','off','Interpreter','none');%text(x,y,'blablabla','Parent',g.axeshandle);
        end
      else
        g.notover
      end
    end

    %--------------------------------
    function killtext(varargin)
      %--------------------------------
      g=varargin{1};
      % Check mouse location on screen
      mouseloc=get(0, 'PointerLocation');
      guiloc=get(g.fig, 'Position');
      state= inpolygon(mouseloc(1),mouseloc(2),[guiloc(1),guiloc(1)+guiloc(3)],[guiloc(2),guiloc(2)+guiloc(4)]) ;%test

      %if nolonger over parentaxes kill text
      if ~state%isempty(overobj(DATA.fig))%~isequal(hittest(DATA.fig),g.imagehandle) %&& overobj(DATA.fig)
        %if we no longer are over the toolbar, kill timer.
        set(g.texthandle,'visible','off')
        g.displayinfo=0;
        if ~isempty(g.texttimer)
          stop(g.texttimer)
          delete(g.texttimer)
          g.texttimer=[];

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
          'Name','KillTimer',...
          'StartDelay',1);
        start(g.killtimer);
      end
    end

    %--------------------------------
    function notover(varargin)
      %--------------------------------
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
      if ~isempty(g.texttimer)
        stop(g.texttimer)
        delete(g.texttimer)
        g.texttimer=[];
      end
    end

    %--------------------------------
    function textoff(g)
      %--------------------------------
      set(g.texthandle,'visible','off')
      g.displayinfo=0;
      if ~isempty(g.texttimer)
        stop(g.texttimer)
        delete(g.texttimer)
        g.texttimer=[];
      end
    end

    %--------------------------------
    function displayinfoYN(varargin)
      %--------------------------------
      %Should we display info now that timer has triggered

      g=varargin{1};
      if  isequal(hittest(g.fig),g.imagehandle) && ~g.buttonisdown && ~g.isribbon
        g.displayinfo=1;
        set(g.texthandle,'visible','on')
        overicon=g.geticon;
        xy = get(g.axeshandle, 'CurrentPoint');
        %         %we need to place the the text better as it is in the way of
        %         %clicking try adding a small offset also added so that it
        %         switches direction in the last 20 % of axes.
        x=xy(1,1);
        y=xy(1,2);

        xl = xlim(g.axeshandle);
        yl = ylim(g.axeshandle);

        if y/yl(2)<0.4
          verticalalignment = 'top';
        else
          verticalalignment = 'bottom';
        end


        if g.pull2==2 || x/xl(2)>0.8%pull2
          horizontalalignment = 'right';
          offset = -10;
        else
          horizontalalignment = 'left';
          offset = 10;
        end
        txt = gettoolinfo(overicon.name);
        if isempty(txt)
          txt = translation.dictionary(overicon.mouseovertext);
        end
        set(g.texthandle,'Position',[x+offset,y,0],'String',txt,'verticalalignment',verticalalignment,'horizontalalignment',horizontalalignment);%text(x,y,'blablabla','Parent',g.axeshandle);
        
        %         ht = text(0,0,'*text', 'tag', 'rollover');
        %         pointerBehavior.enterFcn = ...
        %                 @(hfig, currentPoint)set(findobj(hfig, 'tag', 'rollover'), ...
        %                 'string', 'text - here''s some more info...');
        %               pointerBehavior.traverseFcn = [];
        % % The exitFcn is similar to the enterFcn, but it changes the string back to
        % % the shorter version...
        %       pointerBehavior.exitFcn = ...
        %           @(hfig, currentPoint)set(findobj(hfig, 'tag', 'rollover'), ...
        %           'string', 'text*');
        %       % Now, I need to link the pointer behavior to the object (the text box):
        %       iptSetPointerBehavior(ht, pointerBehavior);
        %       % Now, I need to enable pointer management for the figure:
        %       iptPointerManager(gcf, 'enable');
        %Start kill text timer that checks for a while if we are over the axis if not kill text
        %if isempty(g.killtimer)

        g.killtimer=timer('ExecutionMode','singleShot',...
          'TimerFcn',@g.killtext, ...
          'Name','KillText',...
          'StartDelay',0.2);
        start(g.killtimer);
        %end
      else
        g.displayinfo=0;
      end
    end

    %--------------------------------
    function kill(varargin)
      %--------------------------------
      g = varargin{1};
      if ~isempty(g.dropdowniconholders)
        delete(g.dropdowniconholders)
      end
      delete(g)
    end

    %--------------------------------
    function hidedropdown(g)
      %--------------------------------
      for loop = 1:length(g.iconCell)
        if ~isempty(g.iconCell{loop}.dropdownpanel)
          try
            g.iconCell{loop}.dropdownpanel.Visible = 'off';
          catch
          end
          %g.iconCell{loop}.isindented = 0;
          %g.iconCell{loop}.undent;
        end
      end
    end

  end
end