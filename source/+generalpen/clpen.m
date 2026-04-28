classdef (Abstract) clpen < handle
  % This abstract class provides an interface to implement pen tools used
  % for segmentation. The abstract class needs to be implemented using
  % subclasses designed for specific applications, for example a subclass
  % that implements a general segmentation pen.
  % One instance of the class corresponds to one segmentated object.
  %
  % Justine Le Douaron, Medviso, 2023

  properties (SetAccess = public, GetAccess = public)
    X
    Y
    InterpX
    InterpY
    InterpOngoing = false
  end
  properties (SetAccess = protected, GetAccess = public)
    NumPoints {mustBeInteger,mustBePositive}
    Label %{mustBeText}
    Color
    Visible {mustBeMember(Visible,{'on','off'})} = "on"
    Selected {mustBeNumericOrLogical} = false
  end

  %   methods (Abstract, Access = public) %Methods that needs concrete implementation in subclasses
  %     penbuttonup
  %   end
  methods (Access = public)
    % Public access to call these methods from other functions,
    % can be overriden in subclasses, since they are NOT sealed

    %------------------------------------
    function self = deletetimeframes(self,ind2keep)
      %----------------------------------
      % Removes time frames, called from tools('removetimeframes'),
      % provided indexes are indexes that are kept

      self.X = self.X(:,ind2keep,:);
      self.Y = self.Y(:,ind2keep,:);
      self.InterpX = self.InterpX(ind2keep,:);
      self.InterpY = self.InterpY(ind2keep,:);
    end

    %------------------------------------
    function self = deleteslices(self,ind2keep)
      %----------------------------------
      % Removes slices, called from tools('removesliceshelper'),
      % provided indexes are indexes that are kept

      self.X = self.X(:,:,ind2keep);
      self.Y = self.Y(:,:,ind2keep);
      self.InterpX = self.InterpX(:,ind2keep);
      self.InterpY = self.InterpY(:,ind2keep);
    end
  end

  methods (Access = public, Sealed = true)
    % Public access to call these methods from other functions,
    % SEALED can not be overriden in subclasses
    %----------------------------------------------
    function self = clpen
      %----------------------------------------------
      %Constructor. Needed for concrete classes' constructors.

    end
 
    %GET methods
    %----------------------------------------------
    function numpoints = getnumpoints(self)
      %----------------------------------------------
      numpoints = self.NumPoints;
    end
    
    %----------------------------------------------
    function label = getlabel(self)
      %----------------------------------------------
      label = self.Label;
    end

    %----------------------------------------------
    function color = getcolor(self)
      %----------------------------------------------
      color = self.Color;
    end

    %----------------------------------------------
    function rgbcolor = getcolorRGB(self)
      %----------------------------------------------
      %Convert the property "Color" from hexadecimal to RGB value      
      rgbcolor = colorfunctions.hex2rgb(self.Color);
    end

    %SET methods
     %----------------------------------------------
     function setnumpoints(self,numpoints)
      %----------------------------------------------
      self.NumPoints = numpoints;
    end

    %----------------------------------------------
    function setlabel(self,label)
      %----------------------------------------------
      self.Label = label;
    end

    %----------------------------------------------
    function setcolor(self,color)
      %----------------------------------------------
      self.Color = color;
    end

    %----------------------------------------------
    function setvisibility(self,state)
      %----------------------------------------------
      arguments
        self
        state {mustBeMember(state,{'on','off'})}
      end
      self.Visible = state;
    end

    %----------------------------------------------
    function select(self,bool)
      %----------------------------------------------
      arguments
        self
        bool {mustBeNumericOrLogical} = true;
      end
      self.Selected = bool;
    end

  end %end sealed public methods

  methods (Access = protected, Sealed = true)
  end %end sealed protected methods

  methods (Access = protected)
  end %end protected methods

  methods (Static, Access = protected)
  end %end protected static methods

  methods (Static, Sealed = true)
  end %end sealed static methods
end