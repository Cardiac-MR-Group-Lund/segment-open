classdef timebarfunctions
  % Collection of functions for handling time/timebar
  
  %#ok<*GVMIS>

  methods (Static)
    %------------------------------------
    function tvec = gettimevectorinms(no)
      %----------------------------------
      % Function to get time vector in milliseconds for all time frames
      global SET
      tvec = 1000*SET(no).TimeVector;
    end

    %------------------------------------
    function [xlimstart,xlimend] = getxlimsinms(no)
      %----------------------------------
      % Function to get x.axis limits in milliseconds
      global SET
      xlimstart = 1000*SET(no).TimeVector(1);
      xlimend = 1000*SET(no).TimeVector(end);
    end

    %---------------------------------------
    function t = getcurrenttimeframeinms(no)
      %-------------------------------------
      % Function to get the time in milliseconds for the current time frame
      global SET
      tind = SET(no).CurrentTimeFrame;
      t = 1000*SET(no).TimeVector(tind);
    end

    %---------------------------------------
    function t = gettimeframeinms(no,tind)
      %-------------------------------------
      % Function to get the time in milliseconds for provided time frame index
      global SET
      t = 1000*SET(no).TimeVector(tind);
    end

    %---------------------------------------
    function t = getcurrenttimeframestring(no)
      %-------------------------------------
      % Function to get current time as string
      global SET
      t = num2str(SET(no).CurrentTimeFrame);
    end

    %---------------------------------------
    function str = getcombinedtimeframestring(no)
      %-------------------------------------
      % Function to get current time as string
      str1 = timebarfunctions.getcurrenttimeframestring(no);
      str2 = timebarfunctions.getendtimeframestring(no);
      str = sprintf('%s/%s',str1,str2);
    end

    %---------------------------------------
    function t = getendtimeframestring(no)
      %-------------------------------------
      % Function to get end time as string
      global SET
      t = num2str(SET(no).TSize);
    end

    %-----------------------------------
    function tf = getclickedtime(obj,no)
      %---------------------------------
      % Function to get the time frame based on the clicked position on a plot
      global SET
      % Check if clicked object is an axes handle, take the parent
      % otherwise
      if ~isa(obj, 'matlab.graphics.axis.Axes')
        obj = obj.Parent;
      end
      % Get the x-coordinate of the clicked position
      x = mygetcurrentpoint(obj);
      % Find the closest time frame index based on the clicked position
      [~,tf] = min(abs(1000*SET(no).TimeVector-x));
      % Ensure found time frame is within valid bounds
      tf = max(min(tf,SET(no).TSize),1);
    end

    %------------------------------
    function setclickedtime(obj,no)
      %----------------------------
      % Function to set the current time frame based on the clicked position
      global SET
      % Get the time frame based on the clicked position
      tf = timebarfunctions.getclickedtime(obj,no);
      % Set the current time frame
      SET(no).CurrentTimeFrame = tf;
    end
  end
end