classdef (Abstract) clgeneralgraphics < matlab.mixin.SetGet
  %Abstract class to implement classes that store an axes in which data can
  % be plotted. The classes also store default parameters for line style,
  % markers, etc. that are best suited for a certain type of plot (graph,
  % bullseye, etc).
  %
  % Justine Le Douaron, Medviso, 2024

  properties
    PlotAxes %parent axes
    ForegroundColor
    BackgroundColor
    LineStyle
    LineWidth
    MarkerStyle
    MarkerSize
    TextHorizontalAlignment
    TextWeight
    TextColor
    FontUnits
    FontSize
  end

  methods
    %--------------------------------------------------------------------
    function self = clgeneralgraphics(axesobject)
      %--------------------------------------------------------------------
      %Constructor
      global DATA

      %Assign parent axes
      self.PlotAxes = axesobject;

      %Use Segment's theme colors
      self.ForegroundColor = DATA.GUISettings.ForegroundColor;
      self.BackgroundColor = DATA.GUISettings.VolumeColorGraph;

      %Set font size
      self.FontUnits = 'pixels';
      self.FontSize = 12;
    end

    %--------------------------------------------------------------------
    function setlinesettings(self,linestyle,linewidth)
      %--------------------------------------------------------------------
      %Assign value to the properties related to Line objects

      self.LineStyle = linestyle;
      self.LineWidth = linewidth;
    end

    %--------------------------------------------------------------------
    function setmarkersettings(self,markerstyle,markersize)
      %--------------------------------------------------------------------
      %Assign value to the properties related to Markers

      self.MarkerStyle = markerstyle;
      self.MarkerSize = markersize;
    end

    %--------------------------------------------------------------------
    function applydefaultlinesettings(self)
      %--------------------------------------------------------------------
      %Apply default line settings to Line objects in the axes

      %Find objects
      lines = self.findlines;

      %Apply settings
      set(lines,'LineWidth',self.LineWidth,'LineStyle',self.LineStyle);
    end

    %--------------------------------------------------------------------
    function applydefaultmarkersettings(self)
      %--------------------------------------------------------------------
      %Apply default marker settings to Line objects in the axes

      %Find objects
      lines = self.findlines;

      %Apply settings
      set(lines,'Marker',self.MarkerStyle,'MarkerSize',self.MarkerSize);
    end

    %--------------------------------------------------------------------
    function lines = findlines(self)
      %--------------------------------------------------------------------
      %Return Line objects found in the axes
      lines = findobj(self.PlotAxes.Children,'Type','Line');
    end

    %--------------------------------------------------------------------
    function applydefaulttextsettings(self)
      %--------------------------------------------------------------------
      %Apply default text settings to Text objects in the axes

      %Find objects
      texts = self.findtexts;

      %Apply settings
      set(texts,'HorizontalAlignment',self.TextHorizontalAlignment, ...
        'Color',self.TextColor,'FontWeight',self.TextWeight);
    end

    %--------------------------------------------------------------------
    function texts = findtexts(self)
      %--------------------------------------------------------------------
      %Return Text objects found in the axes
      texts = findobj(self.PlotAxes.Children,'Type','Text');
    end

    %--------------------------------------------------------------------
    function linestyle = getlinestyle(self)
      %--------------------------------------------------------------------
      %Return the property LineStyle
      linestyle = self.LineStyle;
    end

    %--------------------------------------------------------------------
    function markerstyle = getmarkerstyle(self)
      %--------------------------------------------------------------------
      %Return the property MarkerStyle
      markerstyle = self.MarkerStyle;
    end

    %--------------------------------------------------------------------
    function markersize = getmarkersize(self)
      %--------------------------------------------------------------------
      %Return the property MarkerSize
      markersize = self.MarkerSize;
    end
  end

  methods (Static)
    %--------------------------------------------------------------------
    function handlePropEvents(source,event)
      %--------------------------------------------------------------------
      %Handling of the listeners events based on the affected property
      propname = source.Name;
      self = event.AffectedObject;
      switch propname
        case 'FontSize'
          self.updatefontsize;
      end
    end
  end
end