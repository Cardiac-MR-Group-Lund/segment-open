classdef clstrainplot < graphics.clplot
  %Class that stores an axes in which data can be plotted.
  % The other parameters stored in the class are typical for plotting
  % relevant results of strain analysis.
  %
  % Justine Le Douaron, Medviso, 2024

  properties
    isStrainRate
  end

  methods
    %--------------------------------------------------------------------
    function self = clstrainplot(axesobject,tvec,isstrainrate)
      %--------------------------------------------------------------------
      %Constructor
      arguments
        axesobject
        tvec
        isstrainrate = false
      end

      %Call constructor of superclass
      self = self@graphics.clplot(axesobject,tvec);
      self.isStrainRate = isstrainrate;

      %Assign typical parameters for strain plots
      self.LineStyle = '-';
      self.LineWidth = 1.2;
      self.MarkerStyle = '.';
      self.MarkerSize = 13;

      %Assign axes labels
      if isstrainrate
        ylabelstr = makeunitstring(dprintf('Strain rate'),'1/s');
      else
        ylabelstr = makeunitstring(dprintf('Strain'),'%');
      end
      set(self.PlotAxes.YLabel,'String',ylabelstr);
    end
  end
end