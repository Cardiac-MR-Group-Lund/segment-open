classdef clstrainbullseyeplot < graphics.clbullseyeplot
  %Class that stores an axes in which data can be plotted in a bullseye.
  % The other parameters stored in the class are typical for plotting
  % relevant results of strain analysis.
  %
  % Justine Le Douaron, Medviso, 2024

  properties
    isStrainRate
  end

  methods
    %--------------------------------------------------------------------
    function self = clstrainbullseyeplot(axesobject,datatoplot,isstrainrate)
      %--------------------------------------------------------------------
      %Constructor
      arguments
        axesobject
        datatoplot
        isstrainrate = false
      end

      numslices = 3;
      resolution = 200;

      %Call constructor of main class
      self = self@graphics.clbullseyeplot(axesobject,datatoplot,resolution,numslices);

      self.hasSymmetricRange = true;
      self.NumSectors = 24;
      self.NumSectorsPerSlice = [1 4 6 6];
      self.isStrainRate = isstrainrate;

      if isstrainrate
        self.Unit = '[1/s]';
      else
        self.Unit = '[%]';
      end

      self = self.initbullseye;
    end

    %--------------------------------------------------------------------
    function formatstring = getformatstring(self,value)
      %--------------------------------------------------------------------
      if self.isStrainRate
        formatstring = '%0.2f';
      else
        if abs(value) >= 10
          formatstring = '%0.0f';
        else
          formatstring = '%0.1f';
        end
      end
    end

    %--------------------------------------------------------------------
    function self = generatebullseyecolormap(self)
      %--------------------------------------------------------------------
      %Helper function to get a customised colormap for strain bullseye
      %The colormap is red for negative strain values and blue for positive
      %strain values. Between 0 and maximum values, the colormap is interpolated
      %with n points.

      n = 30; %number of interpolation points

      %Create red colormap
      darkred = [220 0 0]/255;
      lightred = [255 200 200]/255;
      cmapred = self.createinterpolatedcolormap(darkred,lightred,n);

      %Create blue colormap
      darkblue = [0 0 220]/255;
      lightblue = [200 200 255]/255;
      cmapblue = self.createinterpolatedcolormap(lightblue,darkblue,n);

      %Concatenate colormaps
      self.Colormap = [cmapred; cmapblue];
    end
  end
end