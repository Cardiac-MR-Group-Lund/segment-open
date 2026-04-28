classdef clpensettings < handle
  % This class stores general settings for pen objects.
  %
  % Justine Le Douaron, Medviso, 2023

  properties (SetAccess = protected, GetAccess = protected)
    MaxNumObjects {mustBeInteger,mustBePositive} = 10
    StandardLabels containers.Map
    StandardColors containers.Map
  end

  methods (Access = public, Sealed = true)
    %----------------------------------------------
    function self = clpensettings
      %----------------------------------------------
      %Constructor.

      %Initialise the object
      self.StandardLabels = generalpen.clpensettings.getlabelsdictionary;
      self.StandardColors = generalpen.clpensettings.getcolorsdictionary;
    end

    %----------------------------------------------
    function maxnumobjects = getmaxnumobjects(self)
      %----------------------------------------------
      maxnumobjects = self.MaxNumObjects;
    end

    %----------------------------------------------
    function color = getdefaultcolor(self)
      %----------------------------------------------
      %Default color for general pen object
      color = self.StandardColors("Tangerine Yellow");
    end

    %----------------------------------------------
    function argout = getatriumcolor(self,atriumtype,getkey)
      %----------------------------------------------
      %Colors for atrial contours
      arguments
        self
        atriumtype {mustBeMember(atriumtype,{'la','ra'})}
        getkey {mustBeNumericOrLogical} = false
      end

      switch atriumtype
        case 'la'
          key = "Mango Tango";
          color = self.StandardColors(key);
        case 'ra'
          key = "Purple";
          color = self.StandardColors(key);
      end

      if getkey
        argout = key;
      else
        argout = color;
      end
    end

    %----------------------------------------------
    function label = getatriumlabel(self,atriumtype)
      %----------------------------------------------
      %Labels for atrial contours
      arguments
        self
        atriumtype {mustBeMember(atriumtype,{'la','ra'})}
      end

      label = self.StandardLabels(upper(atriumtype));
    end

    %----------------------------------------------
    function colordictionary = getgeneralpencolors(self)
      %----------------------------------------------
      %Return a dictionnary containing the colors available for general pen
      %objects.
      heartparts = {'la','ra'};
      reservedcolors = cell(1,numel(heartparts));
      getkey = true; %Get the keys to be able to compare dictionaries
      for loop = 1:numel(heartparts)
        reservedcolors{loop} = self.getatriumcolor(heartparts{loop},getkey);
      end
      colordictionary = self.getcolorsdictionary;
      remove(colordictionary,reservedcolors);
    end

    %----------------------------------------------
    function labels = getstandardlabels(self)
      %----------------------------------------------
      labels = self.StandardLabels;
    end

    %------------------------------------------------------------
    function varargout = createnewobject(self,no,label,color)
      %------------------------------------------------------------
      %Create a new general pen object and store it in SET
      global SET NO

      numobjects = self.getnumobjects;
      if numobjects == self.getmaxnumobjects
        %TODO warning to user that max number of objects is reached
        return
      end

      if nargin < 2
        no = NO;
      end

      if nargin < 3
        label = self.getdefaultlabel;
      end

      if nargin < 4
        color = self.getanyunusedcolor;
      end

      %Create a new general pen object and store it in SET
      ind = numobjects + 1;
      newobject = generalpen.clgeneralpen(label,color);
      if ind == 1
        SET(no).GeneralPenObjects = newobject;
      else
        SET(no).GeneralPenObjects(ind) = newobject;
      end
      %Set the new object as current
      self.setcurrentobject(ind);

      if nargout == 1
        varargout{1} = ind;
      end
    end

    %------------------------------------------------------------
    function createatriumobject(self,atriumtype,no)
      %------------------------------------------------------------
      %Create a new atrial object and store it in SET
      arguments
        self
        atriumtype {mustBeMember(atriumtype,{'la','ra'})}
        no
      end
      global SET

      %Create a new atrial object and store it in SET
      label = self.getatriumlabel(atriumtype);
      color = self.getatriumcolor(atriumtype);
      newobject = generalpen.clatriumpen(label,color);
      SET(no).(upper(atriumtype)) = newobject;
    end
  end %end sealed public methods

  methods (Access = protected, Sealed = true)
    %----------------------------------------------
    function color = getanyunusedcolor(self)
      %----------------------------------------------
      %Return a standard color in hexadecimal format that is not already
      %used by any of the general pen objects. If all the colors are
      %already used, return the default color.
      global SET NO

      %Define the default color
      color = self.getdefaultcolor;

      %If no other object exists, use the default color
      numobjects = self.getnumobjects;
      if numobjects == 0
        return
      end

      %Find the color that are already used by the other objects
      colors = arrayfun(@(x) x.getcolor, SET(NO).GeneralPenObjects(:), 'UniformOutput', false);

      %Compare them to the standard available colors to find one that is
      %not used
      standardcolors = self.getgeneralpencolors;
      standardcolors = standardcolors.values;
      indmatches = matches(standardcolors,colors);
      standardcolors(indmatches') = [];
      if ~isempty(standardcolors)
        color = standardcolors{1};
      end

    end
  end %end sealed protected methods

  methods (Access = protected)
  end %end protected methods

  methods (Static, Access = private)
    %----------------------------------------------
    function labelsdictionnary = getlabelsdictionary
      %----------------------------------------------
      %TODO: handle translation of labels name
      keys = {...
        'RA', ...
        'LA', ...
        'LAX', ...
        };
      labels = {...
        dprintf('Right atrium'), ...
        dprintf('Left atrium'), ...
        dprintf('Long-axis'), ...
        };
      labelsdictionnary = containers.Map(keys,labels);
    end

    %----------------------------------------------
    function colorsdictionnary = getcolorsdictionary
      %----------------------------------------------
      %TODO: handle translation of colours name
      colorname = [...
        "Tangerine Yellow", ...
        "Olive", ...
        "Mango Tango", ...
        ..."Saddle Brown", ...
        "Myrtle", ...
        "Grey", ...
        "Purple", ...
        "Teal", ...
        "Maroon", ...
        "Salmon", ...
        "Medium Violet Red", ...
        "Dark Sea Green", ...
        "Pale Violet Red", ...
        "Light Slate Blue", ...
        ];
      colorcode = [...
        "#FFCC00", ...
        "#808000", ...
        "#F07800", ...
        ..."#804000", ...
        "#004000", ...
        "#808080", ...
        '#A752FF', ...
        "#008080", ...
        "#800000", ...
        "#FA8072", ...
        "#C71585", ...
        "#8FBC8F", ...
        "#DB7093", ...
        "#8080FF", ...
        ];
      colorsdictionnary = containers.Map(colorname,colorcode);
    end

  end %end protected static methods

  methods (Static, Sealed = true)
    %----------------------------------------------
    function label = getdefaultlabel
      %----------------------------------------------
      %Return a default label for a general pen object
      label = dprintf("New object");
    end

    %----------------------------------------------
    function numobjects = getnumobjects(no)
      %----------------------------------------------
      %Return the number of general pen objects in the current image stack
      global SET NO
      if nargin == 0
        no = NO;
      end
      numobjects = length(SET(no).GeneralPenObjects);
    end

    %----------------------------------------------
    function objectind = getcurrentobject
      %----------------------------------------------
      %Return the index of the currently selected general pen object in the
      %current image stack
      global SET NO
      objectind = [];
      if ~isempty(SET(NO).GeneralPenObjects)
        objectind = find(cat(1,SET(NO).GeneralPenObjects(:).Selected));
      end
    end

    %----------------------------------------------
    function setcurrentobject(objectind)
      %----------------------------------------------
      %Set the general pen object corresponding to the provided index as 
      %the currently selected object in the current image stack
      global SET NO

      %Unselect all objects
      for loop = 1:length(SET(NO).GeneralPenObjects)
        select(SET(NO).GeneralPenObjects(loop),false);
      end

      %Select the desired object
      if ~isempty(objectind)
        select(SET(NO).GeneralPenObjects(objectind));
      end
    end
  end %end sealed static methods
end