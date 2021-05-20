function papillaryslider(arg)
%Adjust slider for papillary muscle detection threshold

global DATA SET NO

if nargin == 0
  arg = 'init';
end

switch arg
  case 'init'
    %init
    DATA.GUI.PapillarySlider = mygui('papillaryslider.fig','blocking'); %Create object and store in global the variable DATA
  case 'slider'
    try
      gui = DATA.GUI.PapillarySlider;
      
      %This solution does not seem to work
      %warning off
      %jFig = get(DATA.fig,'JavaFrame');
      %jFig.requestFocus;
      %warning on

      v = get(gui.handles.papillarythresholdslider,'value');
      if isapproxequal(abs(SET(NO).PapillaryThreshold-v),0.0082)
        v = SET(NO).PapillaryThreshold;
        set(gui.handles.papillarythresholdslider,'value',v);
      else
        SET(NO).PapillaryThreshold = v;
      end
      
      if ~isempty(SET(NO).EndoX)
        lvpeter('segmentestimatepapilaryvolume_Callback');
        figure(gui.fig);
      end
    catch me
      mydispexception(me);
    end
  case 'close'
    try
    DATA.GUI.PapillarySlider = close(DATA.GUI.PapillarySlider);
    catch
      delete(gcbf);
      DATA.GUI.PapillarySlider = [];
    end
end

