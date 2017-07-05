function papillaryslider(arg)
%Adjust slider for papillary muscle detection threshold

global DATA SET NO

if nargin == 0
  arg = 'init';
end;

switch arg
  case 'init'
    %init
    DATA.GUI.PapillarySlider = mygui('papillaryslider.fig','blocking'); %Create object and store in global the variable DATA
  case 'slider'
    try
      gui = DATA.GUI.PapillarySlider;
      
      v = get(gui.handles.papillarythresholdslider,'value');
      SET(NO).PapillaryThreshold = v;
           
      if ~isempty(SET(NO).EndoX)
        lvpeter('segmentestimatepapilaryvolume_Callback');
        figure(gui.fig);
      end;
    catch me
      mydispexception(me);
    end;
  case 'close'
    try
    DATA.GUI.PapillarySlider = close(DATA.GUI.PapillarySlider);
    catch
      delete(gcbf);
      DATA.GUI.PapillarySlider = [];
    end
end;

