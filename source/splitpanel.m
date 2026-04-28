classdef splitpanel < handle
  % graphical object to split the figure into several areas

  % 2024 Jelena Bock, Medviso
  properties
    isDragging = false;
    StartPoint = [];
    StartOffset = [];
    FigureHandle = [];
    Handle = [];
    HandleType = '';
    LeftHandle = [];
    RightHandle = [];
    TopHandle = [];
    BottomHandle = [];
    MotionFcn = '';
    ButtonUpFcn = '';
    MaxHeight = 0.9;
    MinHeight = 0.1;
    MaxWidth = 0.9;
    MinWidth = 0.1
  end
  properties (SetObservable)
    Position = []; % normalized value
  end

  methods
    function self = splitpanel(figureHandle, position,varargin)
      self.FigureHandle = figureHandle;
      self.Position = position; % set original position

      % store original functions, in order to reset them
      self.MotionFcn = figureHandle.WindowButtonMotionFcn;
      self.ButtonUpFcn = figureHandle.WindowButtonUpFcn;

      % parse name-value pair arguments
      defaultmax = 0.9;
      defaultmin = 0.1;
      validationfcn = @(x) assert(isnumeric(x) && (x >= defaultmin) && (x <= defaultmax));
      p = inputParser;
      addParameter(p, 'HandleType', 'vertical', @(x) any(validatestring(x, {'vertical', 'horizontal'})));
      addParameter(p, 'MaxHeight', defaultmax, validationfcn);
      addParameter(p, 'MinHeight', defaultmin, validationfcn);
      addParameter(p, 'MaxWidth', defaultmax, validationfcn);
      addParameter(p, 'MinWidth', defaultmin, validationfcn);
      addParameter(p, 'LeftHandle', []);
      addParameter(p, 'RightHandle', []);
      addParameter(p, 'TopHandle', []);
      addParameter(p, 'BottomHandle', []);
      parse(p, varargin{:});

      % assign handles based on input arguments
      self.LeftHandle = p.Results.LeftHandle;
      self.RightHandle = p.Results.RightHandle;
      self.TopHandle = p.Results.TopHandle;
      self.BottomHandle = p.Results.BottomHandle;
      self.HandleType = p.Results.HandleType;
      self.MinHeight = p.Results.MinHeight;
      self.MaxHeight = p.Results.MaxHeight;
      self.MinWidth = p.Results.MinWidth;
      self.MaxWidth = p.Results.MaxWidth;
      % draw splitpanel
      self.draw;
      % add property listeners
      self.addproplisteners;
    end

    %-------------------------------
    function addproplisteners(self)
      %----------------------------
      %Attach listener to SectorRotation to properly update sectors
      addlistener(self,'Position','PostSet',@splitpanel.handlePropEvents);
    end

    %------------------
    function draw(self)
      %----------------
      global DATA %#ok<GVMIS> 
      clr = DATA.GUISettings.ForegroundColor;

      set(self.FigureHandle,'Units','normalized');

      self.Position = calculateposition(self);
      set(self.FigureHandle,'Units','pixels');
      self.Handle = uicontrol(self.FigureHandle, ...
        'Style','Text', ...
        ...'String','H',...
        ...'FontSize',20,...
        ...'FontWeight','bold','HorizontalAlignment','center','ForegroundColor',[1 0 0],...
        'Units','normalized', ...
        'UserData',['splitpanel-',lower(self.HandleType)],... % this can be used to change pointer in parent figure motion function
        'Enable','inactive', ... % inactive makes uicontrol text react to ButtonDownFcn
        'Position', self.Position, ...
        'BackgroundColor', clr ...
        );
      self.Handle.ButtonDownFcn = @self.startDrag;
    end

    %-----------------------------------------
    function position = calculateposition(self)
      %---------------------------------------
      position = self.Position;
      switch self.HandleType
        case 'horizontal'
          if ~isempty(self.LeftHandle)
            if isa(self.LeftHandle,'splitpanel')
              pos = getposition(self.LeftHandle);
              leftpos = pos(1);
            else
              leftpos = self.LeftHandle.Position(3);
            end
          else
            leftpos = self.FigureHandle.Position(3);
          end
          if ~isempty(self.RightHandle)
            if isa(self.RightHandle,'splitpanel')
              pos = getposition(self.RightHandle);
              rightpos = pos(3);
            else
              rightpos = self.RightHandle.Position(3);
            end
          else
            rightpos = 1;% obj.FigureHandle.Position(3);
          end
          width = abs(leftpos - rightpos);
          position(1) = leftpos;
          if width > 0
            position(3) = width;
          end
        case 'vertical'
          if ~isempty(self.TopHandle)
            if isa(self.TopHandle,'splitpanel')
              pos = getposition(self.TopHandle);
              toppos = pos(4);
            else
              toppos = self.TopHandle.Position(4);
            end
          else
            toppos = self.FigureHandle.Position(2);
          end

          if ~isempty(self.BottomHandle)
            if isa(self.BottomHandle,'splitpanel')
              bpos = getposition(self.BottomHandle);
            else
              bpos = self.BottomHandle.Position;
            end
          else
            bpos = self.FigureHandle.Position;
          end
          bottompos = bpos(2) + bpos(4);
          height = abs(toppos-bottompos);
          if height > 0
            position(4) = 1;%height;
          end
          position(2) = 0;
      end
    end

    %--------------------------------------------------
    function position = validateposition(self,position,ind)
      %-----------------------------------------------
      % check width, use x position
      switch ind
        case 1
          position = min(max(position,self.MinWidth),self.MaxWidth);
        case 2
          position = min(max(position,self.MinHeight),self.MaxHeight);
      end
    end

    %---------------------------
    function updateposition(self)
      %-------------------------

      position = calculateposition(self);
      set(self.Handle,'Position',position);
    end

    %---------------------------
    function startDrag(self, ~, ~)
      %---------------------------
      set(self.FigureHandle,'Units','normalized');
      self.isDragging = true;
      self.StartPoint = get(self.FigureHandle, 'CurrentPoint');
      ind = self.getpositionind;
      pos = self.Handle.Position;
      self.StartOffset = pos(ind);

      set(self.FigureHandle, 'WindowButtonMotionFcn', @self.drag);
      set(self.FigureHandle, 'WindowButtonUpFcn', @self.stopDrag);
    end

    %-----------------------
    function drag(self, ~, ~)
      %---------------------
      persistent updatecounter;
      if isempty(updatecounter)
        updatecounter = 0; % init counter
      end
      updatecounter = updatecounter + 1;

      % reduce graphical update rate to every third drag
      if mod(updatecounter, 3) == 0
        if self.isDragging
          ind = self.getpositionind;
          currentPoint = get(self.FigureHandle, 'CurrentPoint');
          offset = currentPoint(ind) - self.StartPoint;
          newposition = self.StartOffset + offset(ind);
          newposition = self.validateposition(newposition,ind);
          self.Handle.Position(ind) = newposition;
          self.Position = self.Handle.Position;
        end
      end
    end

    %---------------------------
    function stopDrag(obj, ~, ~)
      %---------------------------
      set(obj.FigureHandle,'Units','pixels');
      obj.isDragging = false;
      set(obj.FigureHandle, 'WindowButtonMotionFcn', obj.MotionFcn);
      set(obj.FigureHandle, 'WindowButtonUpFcn', obj.ButtonUpFcn);
    end

    %------------------------------
    function pos = getposition(self)
      %---------------------------
      pos = self.Handle.Position;
    end

    %----------------------------------
    function ind = getpositionind(self)
      %--------------------------------
      switch self.HandleType
        case 'horizontal'
          ind = 2;
        otherwise
          ind = 1;
      end
    end

    %----------------------------------------------------------------
    function [pixwidth, pixheight] = getnormalizedpixels(self,numpix)
      %--------------------------------------------------------------
      % figure size in pixels
      figposinpixels = getpixelposition(self.FigureHandle);

      % figure size in normalized units
      figposnormalized = get(self.FigureHandle, 'Position');

      % width and height conversion factors
      widthfactor = figposnormalized(3) / figposinpixels(3);
      heightfactor = figposnormalized(4) / figposinpixels(4);

      % convert number of pixels to normalized units
      pixwidth = numpix * widthfactor;
      pixheight = numpix * heightfactor;
    end
  end

  methods (Static)
    %-------------------------------
    function handlePropEvents(~,evnt)
      %-------------------------------
      % handling of the listeners events
      self = evnt.AffectedObject;
      % check if function handle
      fcn = self.FigureHandle.ResizeFcn;
      if isa(fcn,'function_handle')
        fcn();
      else
        eval(fcn);
      end
    end
  end
end