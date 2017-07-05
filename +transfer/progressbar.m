function progressbar(f,msg)
global DATA

try
  gui = DATA.GUI.Segment;

  set(gui.handles.overalltext,'string',msg);
  set(gui.handles.overallpatch,'xdata',[0 0 f f]);
  drawnow;
catch %#ok<CTCH>
end;