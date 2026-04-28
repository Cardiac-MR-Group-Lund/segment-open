classdef (Abstract) clplot < graphics.clgeneralgraphics
  %Abstract class to implement classes that store an axes in which graphs
  % can be plotted. The classes also store default parameters for line 
  % style, markers, etc. that are best suited for a certain type of graph 
  % (volume, strain, flow, etc).
  %
  % Justine Le Douaron, Medviso, 2024

  properties
  end

  methods
    %--------------------------------------------------------------------
    function self = clplot(axesobject,tvec)
      %--------------------------------------------------------------------
      %Constructor

      %Call constructor of superclass
      self = self@graphics.clgeneralgraphics(axesobject);

      %Initialise the graph plot
      self.init(tvec);
    end

    %--------------------------------------------------------------------
    function self = init(self,tvec)
      %--------------------------------------------------------------------
      %Initialise plot background, grid and axis

      %plot background
      set(self.PlotAxes,'Color',self.BackgroundColor);

      %axis use foreground color
      set(self.PlotAxes,'XColor',self.ForegroundColor,'YColor',self.ForegroundColor);
      set(self.PlotAxes,'YTickLabelMode','auto','YTickMode','auto');

      %show axes on top and bottom, and left and right side of the plot
      set(self.PlotAxes,'Visible','on','Box',true);

      %set x and y axis limits
      set(self.PlotAxes,'XLim',[min(tvec)-2 max(tvec)*1.005]);
      set(self.PlotAxes,'YLimMode','auto');

      %grid is on and uses foreground color
      grid(self.PlotAxes,'on');
      set(self.PlotAxes,'GridColor',self.ForegroundColor);

      %set x and y labels settings
      xlabelstr = makeunitstring(dprintf('Time'),'ms');
      set(self.PlotAxes.XLabel,'String',xlabelstr);
      set(self.PlotAxes.YLabel,'String','');
      set([self.PlotAxes.YLabel self.PlotAxes.XLabel],...
        'Visible','on','Color',self.ForegroundColor,...
        'FontUnits',self.FontUnits,'FontSize',self.FontSize,...
        'Interactions',[]);
    end
  end
end