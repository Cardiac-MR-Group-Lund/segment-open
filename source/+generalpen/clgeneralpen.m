classdef clgeneralpen < generalpen.clpen
  % This class implements a general segmentation pen. 
  % One instance of the class corresponds to one segmentated object.
  %
  % Justine Le Douaron, Medviso, 2023

  properties (SetAccess = protected, GetAccess = public)
    Icon myicon %Icon to be displayed in the objects' toolbar
  end

  methods (Access = public, Sealed = true)
    %----------------------------------------------
    function self = clgeneralpen(label,color)
      %----------------------------------------------
      %Constructor

      %Object instantiation
      self = self@generalpen.clpen; %call superclass' constructor

      %Assign default properties
      self.setlabel(label);
      self.setcolor(color);
      self.setnumpoints(160);

      %Create the icon associated with this object and add it to the
      %objects' toolbar
      self.createicon;
      generalpen.generalpenfunctions('addicontotoolbar',self.Icon);
    end

    %----------------------------------------------
    function createicon(self)
      %----------------------------------------------
      %Create a myicon object to be used in the toolbar that represents
      %General Pen objects
      global DATA

      %Set flags to initialise myicon object
      fdacleared = true;
      defaulttype = 1; %indented, but only one in the group
      toolbargroup = 2; %group
      defaultindent = []; %automatically generates an indented icon

      %Get all necessary arguments for creating the icon
      plaincolorcdata = self.getplaincolorcdata(self.getcolorRGB);
      tooltip = self.getlabel;
      iconname = sprintf('generalpenobject%d',randi([100 999]));
      callback = @() generalpen.generalpenfunctions('setcurrentobject_Callback');

      %Create the icon and store it in the Icon property
      icon{1,1} = myicon(iconname,DATA.Handles.configiconholder,plaincolorcdata,tooltip,callback,defaulttype,toolbargroup,defaultindent,fdacleared);
      self.Icon = icon{1,1};
    end

    %----------------------------------------------
    function updateicon(self)
      %----------------------------------------------
      %Update the properties of the object's icon
      self.Icon.cdata = self.getplaincolorcdata(self.getcolorRGB);
      self.Icon.cdataDisplay = self.Icon.cdata;
      self.Icon.generateclickeddisabledandindent(self.Icon.cdata);
      self.Icon.mouseovertext = self.getlabel;
    end

  end %end sealed public methods

  methods (Access = protected, Sealed = true)
  end %end sealed protected methods

  methods (Access = public)
  end %end public methods

  methods (Static, Access = protected)
  end %end protected static methods

  methods (Static, Access = private, Sealed = true)
    %----------------------------------------------
    function plaincolorcdata = getplaincolorcdata(colorRGB)
      %----------------------------------------------
      %Method that return a RGB CData that is used for myicon objects. The
      %CData represents a square filled with a plain color. 

      %Define the width of the icon's frame
      framewidth = 20;

      %Create a plain color CData
      plaincolorcdata = uint8(240*ones(128,128,3)); %background color
      plaincolorcdata(:,:,1) = colorRGB(1); % Set red channel
      plaincolorcdata(:,:,2) = colorRGB(2); % Set green channel
      plaincolorcdata(:,:,3) = colorRGB(3); % Set blue channel

      %Add a frame
      plaincolorcdata(1:framewidth,:,:) = 240;
      plaincolorcdata(end-framewidth:end,:,:) = 240;
      plaincolorcdata(:,1:framewidth,:) = 240;
      plaincolorcdata(:,end-framewidth:end,:) = 240;
    end
  end %end sealed private static methods
end