function setupicon(fig)
%Set up software icon for mygui figures. 
%Use Segment's icons instead of Mathworks icon.
%
%Note: this function uses the undocumented JavaFrame property which
%may become obsoleted in a future Matlab release (still supported in
%Matlab 2021a). This produces a warning message in the Command Window,
%which is disabled.

global DATA

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

if strcmp(prgmname, 'Segment')
  if isdeployed
    file = 'icon_segment.png';
  else
    file = fullfile([DATA.SegmentFolder filesep 'segment_resources'],'icon_segment.png');
  end  
  
elseif strcmp(prgmname, 'Segment CMR') 
  if isdeployed
    file = 'icon_cmr.png';
  else
    file = fullfile([DATA.SegmentFolder filesep 'segmentcmr_resources'],'icon_cmr.png');
  end  
  
elseif strcmp(prgmname, 'Segment CT')
  if isdeployed
    file = 'icon_ct.png';
  else
    file = fullfile([DATA.SegmentFolder filesep 'segmentct_resources'],'icon_ct.png');
  end 
  
elseif strcmp(prgmname, 'Segment 3DPrint')
  if isdeployed
    file = 'icon_3dp.png';
  else
    file = fullfile([DATA.SegmentFolder filesep 'segment3dp_resources'],'icon_3dp.png');
  end
end

%--- Set up icon
%Disable warning message about JavaFrame
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
%Access the underlying Java object
h = get(handle(fig),'JavaFrame');
%Set new icon
icon = javax.swing.ImageIcon(file);
h.setFigureIcon(icon);
%Enable warning message about JavaFrame
warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');