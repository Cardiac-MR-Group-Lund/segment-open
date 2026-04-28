function papillaryslider(arg)
%Adjust slider for papillary muscle detection threshold

global DATA SET NO %#ok<GVMIS>

if nargin == 0
  arg = 'init';
end

switch arg
  case 'init'
    %init
    DATA.GUI.PapillarySlider = mygui('papillaryslider.fig','blocking'); %Create object and store in global the variable DATA
    gui = DATA.GUI.PapillarySlider;
    % set maingui's keypressfcn for current figure as well
    set(gui.fig,'keypressfcn',DATA.fig.KeyPressFcn)
    if ~isempty(SET(NO).PapillaryThreshold)
      threshold = SET(NO).PapillaryThreshold;
    else
      threshold = 0;
    end
    set(gui.handles.papillarythresholdslider,'value',threshold);
  case 'slider'
    try
      gui = DATA.GUI.PapillarySlider;

      % Request focus for the whole papillary slider figure,
      % resulting that key arrows cannot be used to move slider
      myrequestfocus(gui.fig); %Request focus to avoid problem with sliders
      
      v = get(gui.handles.papillarythresholdslider,'value');

      SET(NO).PapillaryThreshold = v;
      if ~isempty(SET(NO).EndoX)
        lvrvtools('segmentestimatepapilaryvolume_Callback');
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

