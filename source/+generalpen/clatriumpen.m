classdef clatriumpen < generalpen.clpen
  % This class implements a pen for atrial segmentation. 
  % One instance of the class corresponds to one segmentated object.
  %
  % Justine Le Douaron, Medviso, 2023

  properties (SetAccess = public, GetAccess = public)
    isManualContourinES {mustBeNumericOrLogical} = false %flag used in Strain MITT
  end

  methods (Access = public, Sealed = true)
    %----------------------------------------------
    function self = clatriumpen(label,color)
      %----------------------------------------------
      %Constructor

      %Object instantiation
      self = self@generalpen.clpen; %call superclass' constructor
      if nargin > 0
        %Assign default properties
        self.setlabel(label);
        self.setcolor(color);
      end
      self.setnumpoints(160);
    end

    %----------------------------------------------
    function objCopy = copy(self)
      %----------------------------------------------
      % Custom copy method

      label = self.Label;
      color = self.Color;

      % generate new instance of the object
      objCopy = generalpen.clatriumpen(label,color);

      % Copy relevant properties from original object to object copy
      objCopy.X = self.X;
      objCopy.Y = self.Y;
      objCopy.InterpX = self.InterpX;
      objCopy.InterpY = self.InterpY;
    end
  end
end