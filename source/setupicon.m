function setupicon(fig)
%Set up software icon for mygui figures. 
%Use Segment's icons instead of Mathworks icon.
%
%Note: this function uses the undocumented JavaFrame property which
%may become obsoleted in a future Matlab release (still supported in
%Matlab 2021a). This produces a warning message in the Command Window,
%which is disabled.

if ~strcmp(version('-release'),'2022a')
  logdisp('Attempted to change icon, ignored for now in future versions.')
  return
end

global DATA %#ok<*GVMIS> 

%--- Get icon according to software in use
%No icon if no software is running, e.g. software testing or translating
if isempty(DATA) %|| ~isfield(DATA,'ProgramName')
  return;
end
try 
  prgmname = DATA.ProgramName;
catch
  return
end

if isdeployed
  file = 'icon.png';
else 
  if strcmp(prgmname, 'Segment')
  %--- EiH: 2022-03-07. Rewriting to easier code and ensure always return one icon.
  file = fullfile('segment_resources','icon.png');
  end
  if strcmp(prgmname, 'Segment CMR')
    file = fullfile('segmentcmr_resources','icon.png');
  end

  if strcmp(prgmname, 'Segment CT')
    file = fullfile('segmentct_resources','icon.png');
  end

  if strcmp(prgmname, 'Segment 3DPrint')
    file = fullfile('segment3dp_resources','icon.png');
  end
end

%--- Set up icon
if isa(fig,'matlab.apps.AppBase')
  fig.figure1.Icon = file;
else
  %Disable warning message about JavaFrame
  warning('off','MATLAB:ui:javaframe:PropertyToBeRemoved');
  warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame'); %MATLAB 2019
  %Access the underlying Java object
  h = get(handle(fig),'JavaFrame');
  %Set new icon
  icon = javax.swing.ImageIcon(file);
  h.setFigureIcon(icon);
  %Enable warning message about JavaFrame
  warning('on','MATLAB:ui:javaframe:PropertyToBeRemoved');
  warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
end