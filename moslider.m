function moslider(arg)
%Adjust slider for MO threshold

global DATA SET NO

if nargin == 0
  arg = 'init';
end;

switch arg
  case 'init'
    %init
    DATA.GUI.MOSlider = mygui('moslider.fig','blocking'); %Create object and store in global the variable DATA
  case 'slider'
    try
      gui = DATA.GUI.MOSlider;
      
      v = get(gui.handles.mothresholdslider,'value');
      if ~isempty(SET(NO).Scar)
        SET(NO).Scar.MOThreshold = v;
        viability('viabilitycalc');
        drawfunctions('drawimagepanel',DATA.CurrentPanel);
        figure(gui.fig);
      end;
    catch me
      mydispexception(me);
    end;
  case 'close'
    try
      DATA.GUI.MOSlider = close(DATA.GUI.MOSlider);
    catch
      delete(gcbf);
      DATA.GUI.MOSlider = [];
    end
end;

