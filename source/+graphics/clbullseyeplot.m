classdef (Abstract) clbullseyeplot < graphics.clgeneralgraphics
  %Abstract class to implement classes that store an axes in which a
  % bullseye can be plotted. The classes also store the bullseye values, as
  % well as default parameters for number of slices and sectors, values
  % range, colormap, etc.
  %
  % Justine Le Douaron, Medviso, 2024

  properties
    DataToPlot
    Unit
    NumSlices
    NumSectors
    NumSectorsPerSlice
    Resolution
    Scale
    MinValue
    hasSymmetricRange
    Colormap
    Indices
    Rad %polar coordinates
    FrameColor
  end

  properties (SetObservable)
    MaxValue = [];
  end

  methods
    %--------------------------------------------------------------------
    function self = clbullseyeplot(axesobject,datatoplot,resolution,numslices)
      %--------------------------------------------------------------------
      %Constructor
      global DATA

      %Call constructor of superclass
      self = self@graphics.clgeneralgraphics(axesobject);

      %Assign
      self.PlotAxes = axesobject;
      self.DataToPlot = datatoplot;
      self.Resolution = resolution;
      self.NumSlices = numslices;
      self.Scale = self.Resolution/(self.NumSlices + 1);
      self.LineWidth = 1.5;
      self.LineStyle = '-';
      self.TextHorizontalAlignment = 'center';
      self.TextWeight = 'bold';
      self.TextColor = DATA.GUISettings.ForegroundColor;
      self.FrameColor = DATA.GUISettings.VolumeColorGraph;

      %Init property listeners
      self.initlisteners;
    end

    %--------------------------------------------------------------------
    function initlisteners(self)
      %--------------------------------------------------------------------
      %Initialise property listerners
      addlistener(self,'MaxValue','PostSet',@graphics.clbullseyeplot.handlePropEvents);
    end

    %--------------------------------------------------------------------
    function self = updatemaxvalue(self)
      %--------------------------------------------------------------------
      %Perform needed operations when max value is changed

      %Update min value if the colorscale is supposed to be symmetrical
      if self.hasSymmetricRange
        self.MinValue = -self.MaxValue;
      end

      %Update colorbar accordingly
      self.setcolorbarlimits;
    end

    %--------------------------------------------------------------------
    function self = initbullseye(self)
      %--------------------------------------------------------------------
      %Initialise the plot and display the bullseye values and colormap

      %Generate the colormap used for Strain analysis
      self = self.generatebullseyecolormap;

      %Draw bullseye values and plot
      self.updatebullseye;

      %Add a colorbar next to the plot
      self.addcolorbar;
    end

    %--------------------------------------------------------------------
    function self = updatebullseye(self)
      %--------------------------------------------------------------------
      %Initialise the plot and display the bullseye values and colormap

      %Update max/min values and colorbar limits
      self = self.assignmaxminvalues;

      %Draw bullseye
      self.drawbullseye;

      %Draw the strain values in the bullseye
      self.drawbullseyevalues;
    end


    %--------------------------------------------------------------------
    function self = updatedata(self,data)
      %--------------------------------------------------------------------
      %Update the bullseye's data and display the bullseye values

      %Update bullseye's data
      self.DataToPlot = data;

      %Update bullseye values and plot
      self.updatebullseye;
    end

    %--------------------------------------------------------------------
    function drawbullseyevalues(self)
      %--------------------------------------------------------------------
      %Draw the strain values in the bullseye

      %Generate circular data for plotting
      [xcoord,ycoord] = self.generatecircularcoordinates;

      %Get number of sectors for each slice
      numsectors = self.NumSectorsPerSlice;

      %Loop over each slice in the bullseye
      for currentslice = 1:(self.NumSlices + 1)
        %Loop over each sector in the slice
        for currentsector = 1:numsectors(currentslice)

          %Find the indices for the current sector in the current slice
          switch numsectors(currentslice)
            case 6 %basal and mid
              indices = 1 + (currentsector-1)* 4 : 4 + (currentsector-1)*4;
            case 4 %apical
              if currentsector == 4
                indices = [1:2 21:24];
              else
                indices = 3 + (currentsector-1)*6 : 8 + (currentsector-1)*6;
              end
            otherwise %apex
              indices = 1:24;
          end

          %Calculate the average value for the current sector in the
          %current slice
          value = (mynanmean(self.DataToPlot(indices,currentslice)));

          %Display the value as text on the bullseye plot
          if ~isnan(value)
            %Find the coordinates of the current sector
            if currentslice > 2 %basal and mid
              factor = 2/6;
            else %apical
              factor = 1/2;
            end
            ind = max(1,...
              mod(round(currentsector*length(xcoord)/numsectors(currentslice) + factor*length(xcoord)), ...
              length(xcoord)));
            offset = (self.Scale/2 * min(1,max(0,currentslice-1)) + self.Scale*(currentslice-1));
            xposition = self.Resolution + 1 + offset * xcoord(ind);
            yposition = self.Resolution + 1 + offset * ycoord(ind);
            %Get value as a string
            valuestr = self.getvaluestring(value);
            %Create text object
            text(xposition,yposition,valuestr,'Parent',self.PlotAxes)
          end
        end
      end

      %Apply default settings to text in bullseye
      self.applydefaulttextsettings;
    end

    %--------------------------------------------------------------------
    function valuestring = getvaluestring(self,value)
      %--------------------------------------------------------------------
      %Helper function to cast the strain value to a string in order to be
      % displayed in the bullseye plot

      %Get format for casting the value
      formatstr = self.getformatstring(value);

      %Cast the value to a string
      valuestring = sprintf(formatstr,value);

      %Correction to display "0"
      if strcmp(valuestring,'0.00') || strcmp(valuestring,'0.0')
        valuestring = '0';
      end
    end

    %--------------------------------------------------------------------
    function drawbullseye(self)
      %--------------------------------------------------------------------
      %Draw bullseye plot based on the data to plot

      %Update the bullseye properties
      self = self.getbullseyeindices;

      %Generate the image mask and its corresponding transparency mask
      [scaledimagemask,transparencymask] = self.getimagemasks;

      %Display bullseye
      self.drawbullseye_helper(scaledimagemask,transparencymask);
    end

    %--------------------------------------------------------------------
    function [scaledimagemask,transparencymask] = getimagemasks(self)
      %--------------------------------------------------------------------
      %Generate the image mask and its corresponding transparency mask

      %Create an image by mapping bullseye's indices to the data to plot
      scaledimagemask = self.DataToPlot(self.Indices);

      %Temporary assign NaN values for creating the transparency mask
      scaledimagemask(self.Rad > (self.NumSlices+1)) = NaN;

      %Generate a transparency mask that is used later in "imagesc"
      transparencymask = double(not(isnan(scaledimagemask)));

      %Reset NaN values
      scaledimagemask(isnan(scaledimagemask)) = 0;
    end

    %--------------------------------------------------------------------
    function drawbullseye_helper(self,scaledimagemask,transparencymask)
      %--------------------------------------------------------------------
      %Helper function to either draw the bullseye plot from scratch or to
      % update it.

      if isempty(self.PlotAxes.Children)
        %------Create bullseye
        imageobject = imagesc(scaledimagemask,'parent',self.PlotAxes);
        %Set transparancy
        set(imageobject,'AlphaData',transparencymask,'AlphaDataMapping','scaled');
        %Use same length along each axis and fit axes box tightly around data
        axis(self.PlotAxes,'image');
        %Assign colormap to axes
        colormap(self.PlotAxes,self.Colormap);
        %Draw bullseye's frame
        self.drawbullseyeframe;
      else
        %%------Update bullseye
        %Delete all text
        textobjects = findobj(self.PlotAxes.Children,'Type','text');
        delete(textobjects);
        %Update CDATA only
        imageobject = findobj(self.PlotAxes.Children,'Type','image');
        imageobject.CData = scaledimagemask;
        imageobject.AlphaData = transparencymask;
      end

      %Hide axis
      set(self.PlotAxes,'Visible','off');
    end

    %--------------------------------------------------------------------
    function addcolorbar(self)
      %--------------------------------------------------------------------
      %Add colorbar to plot

      %Add colorbar to axes
      colorbarobject = colorbar(self.PlotAxes);

      %Set axis and ticks colors
      set(colorbarobject,'Color',self.TextColor);

      %Set unit of plotted data
      title(colorbarobject,self.Unit,...
        'Color',self.TextColor,'Units','normalized',...
        'Position',[2.5 0.4775 0]);

      %Adjust colorscal to min and max values
      self.setcolorbarlimits;
    end

    %--------------------------------------------------------------------
    function setcolorbarlimits(self)
      %--------------------------------------------------------------------
      %Adjust colorscal to min and max values
      clim(self.PlotAxes,[self.MinValue self.MaxValue]);
    end

    %--------------------------------------------------------------------
    function self = assignmaxminvalues(self)
      %--------------------------------------------------------------------
      %Function to assgin th properties MinValue and MaxValue based on the
      % data to plot

      %Find maxima and choose the maximal absolute value
      tempminvalue = min(self.DataToPlot(:));
      tempmaxvalue = max(max(self.DataToPlot(:)),abs(tempminvalue));
      if isnan(tempmaxvalue) || tempmaxvalue == 0
        self.MaxValue = 1;
      else
        self.MaxValue = tempmaxvalue;
      end

      %Min value is opposite of absolute max value is the data range has to
      % be symmetrical
      if self.hasSymmetricRange
        self.MinValue = -self.MaxValue;
      else
        self.MinValue = -tempmaxvalue;
      end
    end

    %--------------------------------------------------------------------
    function formatstring = getformatstring(value)
      %--------------------------------------------------------------------
      %Return the format to be used based on the value range

      if abs(value) >= 10
        formatstring = '%0.0f';
      else
        formatstring = '%0.1f';
      end
    end

    %--------------------------------------------------------------------
    function drawbullseyeframe(self)
      %--------------------------------------------------------------------
      %Draw the bullseye frame

      linecolor = self.FrameColor;
      scale = self.Scale;

      %---Draw circles
      %Generate circular data for plotting
      [xcoord,ycoord] = self.generatecircularcoordinates(scale);
      %Draw circles
      hold(self.PlotAxes,'on');
      for radius = 1:(self.NumSlices + 1) %number of slices + center (apex)
        %Plot circles with varying radius
        plot(self.PlotAxes, ...
          self.Resolution + 1 + xcoord * radius, ...
          self.Resolution + 1 + ycoord * radius, ...
          'Color',linecolor);
      end

      %---Draw lines
      b = sqrt(0.75);
      a = 0.5;
      c = 1/sqrt(2);
      xcoord = scale * [...
        [0 2];...
        [6 8];...
        [4-c 4-2*c];...
        [4+c 4+2*c];...
        [4-c 4-2*c];...
        [4+c 4+2*c];...
        [4-4*a 4-2*a];...
        [4-4*a 4-2*a];...
        [4+4*a 4+2*a];...
        [4+4*a 4+2*a];...
        ];

      ycoord = scale * [...
        [4 4];...
        [4 4];...
        [4-c 4-2*c];...
        [4+c 4+2*c];...
        [4+c 4+2*c];...
        [4-c 4-2*c];...
        [4-4*b 4-2*b];...
        [4+4*b 4+2*b];...
        [4+4*b 4+2*b];...
        [4-4*b 4-2*b];...
        ];

      for loop = 1:length(xcoord)
        plot(self.PlotAxes,xcoord(loop,:),ycoord(loop,:),'Color',linecolor);
      end
      hold(self.PlotAxes,'off');

      %Apply settings to Line pbjects
      self.applydefaultlinesettings;
    end

    %--------------------------------------------------------------------
    function self = getbullseyeindices(self)
      %--------------------------------------------------------------------
      %Get indices for bullseye plot

      %Generate a row vector
      gridvector = linspace(-self.NumSlices - 1, self.NumSlices + 1, 2*self.Resolution + 1);
      %Generate a 2-D grid
      [x,y] = ndgrid(gridvector);
      %Convert from Cartesian coordinates to polar coordinates
      self.Rad = sqrt(x.*x+y.*y);

      %Create idx outer
      ang = angle(complex(y,x)) + pi;
      idxouter = self.getindices(ang);

      %Create idx inner
      ang = mod(angle(complex(y,x)) + pi + pi/4, 2*pi);
      idxinner = self.getindices(ang);

      %Assign
      self.Indices = idxouter;
      self.Indices(self.Rad < 2) = idxinner(self.Rad < 2);
    end

    %--------------------------------------------------------------------
    function indices = getindices(self,ang)
      %--------------------------------------------------------------------
      %Helper function fot getbullseyeindices

      ang = self.NumSectors * ang / (2*pi);

      %orient it according to sectors in AHA 17-segment model
      ang = mod(-(ang - self.NumSectors/3), self.NumSectors);

      indices = 1 + min(floor(ang),(self.NumSectors - 1)) + ...
        (self.NumSectors) * min(floor(self.Rad),self.NumSlices);
    end

    %--------------------------------------------------------------------
    function self = generatebullseyecolormap(self)
      %--------------------------------------------------------------------
      %Stores the colormap. Overloaded in subclasses if a specific colormap
      % needs to be generated.

      self.Colormap = 'jet';
    end
  end

  methods (Static)
    %--------------------------------------------------------------------
    function handlePropEvents(src,event)
      %--------------------------------------------------------------------
      %Handling of the listeners events based on the affected property
      propname = src.Name;
      self = event.AffectedObject;
      switch propname
        case 'MaxValue'
          self.updatemaxvalue;
      end
    end

    %--------------------------------------------------------------------
    function cmap = createinterpolatedcolormap(color1,color2,n)
      %--------------------------------------------------------------------
      %Creates a colormap with upper value as color1 and lower value as color2.
      %The colormap values are then interpolated to values between the max
      %values.

      %temporary colormap with the max values
      cmap(1,:) = color1;
      cmap(2,:) = color2;

      %mesh of indices
      [X,Y] = meshgrid(1:3,1:n);
      %interpolate colormap
      cmap = interp2(X([1,n],:),Y([1,n],:),cmap,X,Y);
    end

    %--------------------------------------------------------------------
    function [xcoord,ycoord] = generatecircularcoordinates(scale)
      %--------------------------------------------------------------------
      %Generate circular data for plotting in a bullseye
      arguments
        scale = 1
      end

      theta = linspace(0,2*pi,100); % 360 degrees in radians = a circle
      xcoord = scale * sin(theta); % cartesian x coordinates of the circle
      ycoord = scale * cos(theta); % cartesian y coordinates of the circle
    end
  end
end