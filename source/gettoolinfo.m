function [tooltip,fullinfo,hotkey] = gettoolinfo(toolid)
% function to get tooltip, full description and hotkey based on tool ID,
% that can be icon name or hotkey

%#ok<*GVMIS>
global DATA
tooltip = '';
fullinfo = '';
hotkey = '';
shortinfo = '';

if length(toolid) == 1
  % just one letter, then take upper case as hotkey
  hotkey = upper(toolid);
end

switch lower(toolid)
  case '1'
    shortinfo = dprintf('Shift mode in panel between one and all frames mode');
  case 'a'
    shortinfo = dprintf('Select Analysis tab');
  case 'angio'
    shortinfo = dprintf('Angio view');
  case 'alt-a'
    hotkey = sprintf('%s-A',getaltkey);
    shortinfo = dprintf('Translate contours left (selected slices)');
  case 'alt-d'
    hotkey = sprintf('%s-D',getaltkey);
    shortinfo = dprintf('Set end diastole at current time frame');
  case 'alt-f'
    hotkey = sprintf('%s-F',getaltkey);
    switch DATA.ProgramName
      case {'Segment CMR', 'Segment'}
        shortinfo = dprintf('Propagate Flow ROI forward and refine');
%       case 'Segment 3DPrint'
%         shortinfo = dprintf('View from front');
%       case 'Segment CT'
    end
  case 'alt-s'
    hotkey = sprintf('%s-S',getaltkey);
    shortinfo = dprintf('Set end systole at current time frame');
  case 'alt-t'
    hotkey = sprintf('%s-T',getaltkey);
    switch DATA.ProgramName
      case {'Segment CMR', 'Segment'}
        shortinfo = dprintf('Track vessel over time');
%       case 'Segment 3DPrint'
%         shortinfo = dprintf('View from top');
%       case 'Segment CT'
    end
  case 'alt-up'
    hotkey = sprintf('%s-%s',getaltkey,getuparrow);
    switch DATA.ProgramName
      case {'Segment CMR', 'Segment'}
        shortinfo = dprintf('Even out myocardium wall');
%       case 'Segment 3DPrint'
%         shortinfo = dprintf('Rotate 3D view up');
%       case 'Segment CT'
    end
  case {'alt-right','faster'}
      hotkey = sprintf('%s-%s',getaltkey,getrightarrow);
    switch DATA.ProgramName
      case {'Segment CMR', 'Segment'}
        shortinfo = dprintf('Faster frame rate');
    end
  case {'alt-left','slower'}
      hotkey = sprintf('%s-%s',getaltkey,getleftarrow);
    switch DATA.ProgramName
      case {'Segment CMR', 'Segment'}
        shortinfo = dprintf('Slower frame rate');
    end
  case 'alt-w'
    hotkey = sprintf('%s-W',getaltkey);
    shortinfo = dprintf('Translate contours up (selected slices)');
  case 'alt-x'
    hotkey = sprintf('%s-X',getaltkey);
    shortinfo = dprintf('Translate contours right (selected slices)');
  case 'alt-z'
    hotkey = sprintf('%s-Z',getaltkey);
    shortinfo = dprintf('Translate contours down (selected slices)');
  case 'autolax'
    shortinfo = dprintf('AI-based automatic LAX segmentation for ED (all chambers)');
  case 'autolv'
    switch DATA.ProgramName
      case {'Segment CMR', 'Segment'}
        return % for now
%       case 'Segment 3DPrint'
%         hotkey = sprintf('%s-L',getctrlkey);
%         shortinfo = dprintf('Automatic CT LV segmentation');
      case 'Segment CT'
        hotkey = sprintf('%s-L',getctrlkey);
        shortinfo = 'Automatic LV segmentation';
    end
  case 'autolvedes'
    hotkey = sprintf('%s-L',getctrlkey);
    shortinfo = dprintf('AI-based automatic LV SAX segmentation for ED/ES');
  case 'autolvlaxwand'
    shortinfo = dprintf('AI-based automatic LV LAX segmentation for ED (all chambers)');
  case 'autolvlaxsemi'
    shortinfo = dprintf('AI-based semi-automatic LV LAX segmentation for ED (all chambers)');
  case 'autolvone'
    shortinfo = sprintf('%s (%s)',dprintf('AI-based automatic LV SAX segmentation'),dprintf('current timeframe'));
  case 'autolvwand'
    hotkey = sprintf('%s-L',getctrlkey);
    shortinfo = dprintf('AI-based automatic LV SAX segmentation');
  case 'autolvlge'
    shortinfo = dprintf('AI-based LV segmentation in LGE images');
  case 'autorvendo'
    hotkey = sprintf('%s-%s-M',getctrlkey,getaltkey);
    shortinfo = dprintf('Semi-automatic RV Endo segmentation');
  case 'autorvwand'
    hotkey = sprintf('%s-%s-M',getctrlkey,getaltkey);
    shortinfo = dprintf('AI-based automatic RV Endo segmentation for ED/ES');
%   case 'autoteeth'
%     shortinfo = sprintf('%s - %s (%s)', dprintf('AI teeth'), dprintf('CT'), dprintf('Research use only'));
%   case 'autovesselbone'
%     shortinfo = sprintf('%s - %s (%s)', dprintf('AI bone vessel'), dprintf('CT'), dprintf('Research use only'));
  case 'autozoom'
    hotkey = sprintf('%s-Z',getshiftkey);
    shortinfo = dprintf('Auto zoom');
  case {'bullseye','alt-b'}
    hotkey = sprintf('%s-B',getctrlkey);
    shortinfo = dprintf('Bullseye analysis');
  case 'c'
    switch DATA.ProgramName
      case {'Segment CMR', 'Segment'}
        shortinfo = dprintf('Play cine thumbnail');
%       case 'Segment 3DPrint'
%         shortinfo = dprintf('Set coronal view as current view');
    end
  case {'closeall','ctrl-shift-w'}
    hotkey = sprintf('%s-%s-W',getctrlkey,getshiftkey);
    shortinfo = dprintf('Close all image stacks');
  case {'contract'}
    shortinfo = dprintf('Contract');
  case {'contractendo','ctrl-k'}
    hotkey = sprintf('%s-K',getctrlkey);
    shortinfo = dprintf('Contract LV Endo');
  case {'contractepi','ctrl-alt-k'}
    hotkey = sprintf('%s-%s-K',getctrlkey,getaltkey);
    shortinfo = dprintf('Contract LV Epi');
  case {'copylvdown','ctrl-d'}
    hotkey = sprintf('%s-D',getctrlkey);
    shortinfo = dprintf('Copy LV downwards and refine');
  case {'copylvup','ctrl-u'}
    hotkey = sprintf('%s-U',getctrlkey);
    shortinfo = dprintf('Copy LV upwards and refine');
  case {'copyrvdown','ctrl-alt-d'}
    hotkey = sprintf('%s-%s-D',getctrlkey,getaltkey);
    shortinfo = dprintf('Copy RV Endo downwards and refine');
  case {'copyrvup','ctrl-alt-u'}
    hotkey = sprintf('%s-%s-U',getctrlkey,getaltkey);
    shortinfo = dprintf('Copy RV Endo upwards and refine');
  case 'ctrl-0'
    hotkey = sprintf('%s-0 (%s)',getctrlkey,dprintf('zero'));
    shortinfo = dprintf('Reset GUIs position');
  case 'ctrl-1'
    hotkey = sprintf('%s-1',getctrlkey);
    switch DATA.ProgramName
      case {'Segment CMR', 'Segment'}
        shortinfo = dprintf('One slice view');
%       case 'Segment 3DPrint'
%         shortinfo = dprintf('2D Split out to one object');
    end
  case 'ctrl-a'
    hotkey = sprintf('%s-A',getctrlkey);
    shortinfo = dprintf('Select all slices');
  case 'ctrl-c'
    hotkey = sprintf('%s-C',getctrlkey);
    shortinfo = dprintf('Apply auto contrast to all image stacks');
  case 'ctrl-alt-m'
    hotkey = sprintf('%s-%s-M',getctrlkey,getaltkey);
    if DATA.isGPUCompatible
      shortinfo = dprintf('AI-based automatic RV Endo segmentation for ED/ES');
    else
      shortinfo = dprintf('Semi-automatic RV Endo segmentation');
    end
  case 'ctrl-l'
    hotkey = sprintf('%s-L',getctrlkey);
    switch DATA.ProgramName
      case {'Segment CMR', 'Segment'}
        shortinfo = dprintf('Perform automatic LV segmentation');
    end
  case 'ctrl-alt-l'
    hotkey = sprintf('%s-%s-L',getctrlkey,getaltkey);
    switch DATA.ProgramName
      case 'Segment'
        shortinfo = sprintf('Perform alternative automatic LV segmentation in selected slices');
    end
  case 'ctrl-n'
    hotkey = sprintf('%s-N',getctrlkey);
    shortinfo = dprintf('Open next .mat file');
  case {'ctrl-o','openfromdisc'}
    hotkey = sprintf('%s-O',getctrlkey);
    shortinfo = dprintf('Open from disc');
  case 'ctrl-q'
    hotkey = sprintf('%s-Q',getctrlkey);
    shortinfo = dprintf('Quit program');
  case 'ctrl-s'
    hotkey = sprintf('%s-S',getctrlkey);
    switch DATA.ProgramName
      case {'Segment CMR', 'Segment 3DPrint','Segment CT'}
        shortinfo = dprintf('Save to patient database');
      case 'Segment'
        shortinfo = dprintf('Save to disc');
    end
  case 'ctrl-mousewheel'
    hotkey = sprintf('%s-%s',getctrlkey,getmousewheel);
    shortinfo = dprintf('Zoom');
  case 'ctrl-shift-mousewheel'
    hotkey = sprintf('%s-%s-%s',getctrlkey,getshiftkey,getmousewheel);
    shortinfo = dprintf('Windowing');
  case 'ctrl-w'
    hotkey = sprintf('%s-W',getctrlkey);
    shortinfo = dprintf('Close current image stack');
  case 'd'
    shortinfo = dprintf('Go to end diastole');
  case {'database','ctrl-p'}
    hotkey = sprintf('%s-P',getctrlkey);
    shortinfo = dprintf('Open patient database');
  case 'databaseadd'
    hotkey = sprintf('%s-S',getctrlkey);
    shortinfo = dprintf('Save to patient database');
  case 'databasepacsadd'
    shortinfo = dprintf('Save to patient database and PACS');
  case 'drawingguide'
    shortinfo = dprintf('Drawing guidance for strain analysis');
  case 'down'
    hotkey = getdownarrow;
    shortinfo = dprintf('Next slice');
    fullinfo = dprintf('View next slice in apical direction');
  case 'evenoutwall'
    hotkey = sprintf('%s-%s',getaltkey,getuparrow);
    shortinfo = dprintf('Even out myocardium wall');
    case 'expand'   
    shortinfo = dprintf('Expand');
  case {'expandendo','ctrl-e'}
    hotkey = sprintf('%s-E',getctrlkey);
    shortinfo = dprintf('Expand LV Endo');
  case {'expandepi','ctrl-alt-e'}
    hotkey = sprintf('%s-%s-E',getctrlkey,getaltkey);
    shortinfo = dprintf('Expand LV Epi');
  case 'f'
    shortinfo = dprintf('Select ROI/Flow tab');
%   case 'fill'
%     shortinfo = dprintf('Fill in 2D [F]/3D [%s-F]',getshiftkey);
  %case 'g'
   % shortinfo = dprintf('Select General Pen tab');
  case {'h'}
    hotkey = 'H';
    shortinfo = dprintf('Hide/show all overlays (segmentation, point, text, etc.), remember previous hide setting');
  case 'i'
    shortinfo = dprintf('Select Image tab');
%   case 'keep'
%     shortinfo = dprintf('Select isolated object and create new object');
  case 'l'
    shortinfo = dprintf('Select Function tab');
  case 'm'
    shortinfo = dprintf('Go to next panel');
  case 'mousewheel'
    hotkey = getmousewheel;
    shortinfo = dprintf('Scroll through slices');
  case 'mousewheelbutton'
    hotkey = dprintf('Mouse wheel button');
    shortinfo = dprintf('Pan');
  case 'move'
    shortinfo = dprintf('Translate contour');
  case 'n'
    shortinfo = dprintf('Go to previous panel');
  case {'next','right'}
    hotkey = getrightarrow;
    shortinfo = dprintf('Next frame');
  case 'p'
    switch DATA.ProgramName
      case {'Segment CMR', 'Segment'}
        shortinfo = dprintf('Play movie');
%       case 'Segment 3DPrint'
%         shortinfo = dprintf('Set point to current tool');
    end
  case 'panel1'
    hotkey = sprintf('%s-1',getshiftkey);
    shortinfo = dprintf('View one image panel');
  case {'panel2','shift-2'}
    hotkey = sprintf('%s-2',getshiftkey);
    shortinfo = dprintf('View two image panels');
  case {'panel2x1','alt-2'}
    hotkey = sprintf('%s-2',getaltkey);
    shortinfo = dprintf('View two image panels');
    fullinfo = dprintf('View two image panels as rows');
  case {'panel3','shift-3'}
    hotkey = sprintf('%s-3',getshiftkey);
    shortinfo = dprintf('View three image panels');
  case {'panel3x1','alt-3'}
    hotkey = sprintf('%s-3',getaltkey);
    shortinfo = dprintf('View three image panels');
    fullinfo = dprintf('View three image panels as rows');
  case 'panel4'
    hotkey = sprintf('%s-4',getshiftkey);
    shortinfo = dprintf('View four image panels');
  case {'panel6','shift-6'}
    hotkey = sprintf('%s-6',getshiftkey);
    shortinfo = dprintf('View six image panels');
  case 'play'
    hotkey = 'P';
    shortinfo = dprintf('Play movie');
  case {'plotflow','ctrl-t'}
    hotkey = sprintf('%s-T',getctrlkey);
    shortinfo = dprintf('Plot flow');
  case 'point'
    shortinfo = dprintf('Place annotation point');
    switch DATA.ProgramName
      case 'Segment 3DPrint'
        hotkey = 'P';
    end
  case {'prev','left'}
    hotkey = getleftarrow;
    shortinfo = dprintf('Previous frame');
  case {'propagateendo','ctrl-f'}
    hotkey = sprintf('%s-F',getctrlkey);
    shortinfo = dprintf('Propagate LV Endo forward and refine');
  case {'propagateepi','ctrl-shift-f'}
    hotkey = sprintf('%s-%s-F',getctrlkey,getshiftkey);
    shortinfo = dprintf('Propagate LV Epi forward and refine');
  case 'r'
    shortinfo = dprintf('Select Review tab');
  case {'refineendo','ctrl-r'}
    hotkey = sprintf('%s-R',getctrlkey);
    shortinfo = dprintf('Refine LV Endo');
  case {'refineepi','ctrl-shift-r'}
    hotkey = sprintf('%s-%s-R',getctrlkey,getshiftkey);
    shortinfo = dprintf('Refine LV Epi');
  case {'refineroi','alt-r'}
    hotkey = sprintf('%s-R',getaltkey);
    shortinfo = dprintf('Refine Flow ROI');
  case 'refineroinext'
    hotkey = sprintf('%s-F',getaltkey);
    shortinfo = dprintf('Propagate Flow ROI forward and refine');
  case {'refinervendo','ctrl-alt-r'}
    hotkey = sprintf('%s-%s-R',getctrlkey,getaltkey);
    shortinfo = dprintf('Refine RV Endo');
  case 's'
    shortinfo = dprintf('Go to end systole');
  case 'scale'
    shortinfo = dprintf('Scale object');
  case 'select'
    switch DATA.ProgramName
      case {'Segment CMR', 'Segment','Segment CT'}
        shortinfo = dprintf('Select image stack or object');
%       case 'Segment 3DPrint'
%         if strcmpi(DATA.CurrentTheme,'lv')
%           shortinfo = dprintf('Select image stack or object');
%         else
%           hotkey = 'Y';
%           shortinfo = dprintf('Cursor tool, set cross-hair position');
%         end
    end
%   case 'select2d3d'
%     shortinfo = dprintf('Select 2D [2] or 3D [3]');
  case 'selectoneall'
    hotkey = '1';
    shortinfo = dprintf('Select one or all frames mode');
%   case 'set3dviewotherside'
%     hotkey = sprintf('%s-V',getaltkey);
%     shortinfo = dprintf('View 3D from other side');
  case 'shift-1'
    hotkey = sprintf('%s-1',getshiftkey);
    switch  DATA.ProgramName
      case {'Segment CMR', 'Segment','Segment CT'}
        shortinfo = dprintf('View one image panel');
%       case 'Segment 3DPrint'
%         shortinfo = dprintf('View 3D view only');
    end
  case 'shift-4'
    hotkey = sprintf('%s-4',getshiftkey);
    switch  DATA.ProgramName
      case {'Segment CMR', 'Segment','Segment CT'}
        shortinfo = dprintf('View four image panels');
      case 'Segment 3DPrint'
        shortinfo = dprintf('View four image panels');
    end
  case 'shift-a'
    hotkey = sprintf('%s-A',getshiftkey);
    shortinfo = dprintf('View all image stacks');
  case 'shift-alt-a'
    hotkey = sprintf('%s-%s-A',getshiftkey,getaltkey);
    shortinfo = dprintf('Translate contours and image left (selected slices)');
  case 'shift-alt-r'
    hotkey = sprintf('%s-%s-R',getshiftkey,getaltkey);
    shortinfo = dprintf('Refine LV Endo contour for Alternative LV segmentation method');
  case 'shift-alt-w'
    hotkey = sprintf('%s-%s-W',getshiftkey,getaltkey);
    shortinfo = dprintf('Translate contours and image up (selected slices)');
  case 'shift-alt-x'
    hotkey = sprintf('%s-%s-X',getshiftkey,getaltkey);
    shortinfo = dprintf('Translate contours and image right (selected slices)');
  case 'shift-alt-z'
    hotkey = sprintf('%s-%s-Z',getshiftkey,getaltkey);
    shortinfo = dprintf('Translate contours and image down (selected slices)');
  case 'shift-d'
    hotkey = sprintf('%s-D',getshiftkey);
    shortinfo = dprintf('Go to SAX end diastole in all visible image stacks');
  case {'shift-h','hideall'}
    hotkey = sprintf('%s-H',getshiftkey);
    shortinfo = dprintf('Hide/show all overlays');
    fullinfo = dprintf('Hide/show all overlays (segmentation, point, text, etc.)');
  case 'shift-mousewheel'
    hotkey = sprintf('%s-%s',getshiftkey,getmousewheel);
    shortinfo = dprintf('Scroll through time frames');
  case 'shift-s'
    hotkey = sprintf('%s-S',getshiftkey);
    shortinfo = dprintf('Go to SAX end systole in all visible image stacks');
  case 'shift-u'
    hotkey = sprintf('%s-U',getshiftkey);
    shortinfo = dprintf('Unselect all slices');
  case 'shift-z'
    hotkey = sprintf('%s-Z',getshiftkey);
    shortinfo = dprintf('Auto zoom');
  case 'space'
    hotkey = dprintf('Space');
    shortinfo = dprintf('Toggle tool in toolbar menu (epi/endo, RA/LA)');
  case 'shift-space'
    hotkey = sprintf('%s-%s',getshiftkey,dprintf('Space'));
    shortinfo = dprintf('Toggle tool in toolbar menu (all modes)');
  case 'ctrl-space'
    hotkey = sprintf('%s-%s',getctrlkey,dprintf('Space'));
    shortinfo = dprintf('Toggle between thumbnail modes "All" and "Mode"');
  case {'smooth','o'}
    hotkey = 'O';
    shortinfo = dprintf('Smooth current contour');
  case 't'
    shortinfo = dprintf('Select T1/T2/T2* tab');
  case 'trackingvessel'
    hotkey = sprintf('%s-T',getaltkey);
    shortinfo = dprintf('Track vessel over time');
  case  {'undo','ctrl-z'}
    hotkey = sprintf('%s-Z',getctrlkey);
    shortinfo = dprintf('Undo last operation');
  case 'up'
    hotkey = getuparrow;
    shortinfo = dprintf('Previous slice');
    fullinfo = dprintf('View next slice in basal direction');
  case 'v'
    switch DATA.ProgramName
      case {'Segment CMR', 'Segment'}
        fullinfo = dprintf('Shift mode in panel between montage and one slice');
%       case 'Segment 3DPrint'
%         shortinfo = dprintf('Hide/show segmentation or 3D');
%         fullinfo = dprintf('Hide/show segmentation or 3D view');
    end
%   case 'view3d'
%     hotkey = 'V';
%     shortinfo = dprintf('Hide/show 3D model');
  case 'view4'
    hotkey = sprintf('%s-4',getshiftkey);
    shortinfo = dprintf('View four image panels');
  case 'view2'
    hotkey = sprintf('%s-2',getshiftkey);
    shortinfo = dprintf('View two image panels');     
  case 'view1'
    hotkey = sprintf('%s-1',getshiftkey);
    shortinfo = dprintf('View one image panel');
  case {'viewall','ctrl-3'}
    hotkey = sprintf('%s-3',getctrlkey);
    shortinfo = dprintf('View all slices');
    fullinfo = dprintf('Montage view');
  case 'viewcoronal'
    hotkey = sprintf('%s-C',getshiftkey);
    shortinfo = dprintf('View coronal');
  case 'viewone'
    hotkey = sprintf('%s-1',getctrlkey);
    shortinfo = dprintf('View one slice');
  case {'viewrow','ctrl-4'}
    hotkey = sprintf('%s-4',getctrlkey);
    shortinfo = dprintf('View all slices in 2 rows');
    fullinfo = dprintf('Montage row view');
  case 'viewsagittal'
    hotkey = sprintf('%s-S',getshiftkey);
    shortinfo = dprintf('View sagittal');
  case 'viewtransversal'
    hotkey = sprintf('%s-T',getshiftkey);
    shortinfo = dprintf('View transversal');
  case 'w'
    shortinfo = dprintf('Select Viability tab');
%   case 'wand'
%     hotkey = 'W';
%     shortinfo = dprintf('Magic wand');
%   case 'y'
%     shortinfo = dprintf('Select the Selection tool');
  case 'z'
    shortinfo = dprintf('Select Strain tab');
  case {'zoomin','ctrl-plus'}
    hotkey = sprintf('%s-%s',getctrlkey,getpluskey);
    shortinfo = dprintf('Zoom in');
  case {'zoomout', 'ctrl-minus'}
    hotkey = sprintf('%s-%s',getctrlkey,getminuskey);
    shortinfo = dprintf('Zoom out');
  case 'ribbonlv'
    hotkey = 'L';
    shortinfo = dprintf('LV');
  otherwise
    return

end
if ~isempty(hotkey) && ~(strcmp(DATA.ProgramName,'Segment 3DPrint') && strcmpi(DATA.CurrentTheme,'lv'))
  tooltip = sprintf('%s [%s]',shortinfo,hotkey);
else
  tooltip = shortinfo;
end
if isempty(fullinfo)
  fullinfo = shortinfo;
end

function outputkey = getaltkey
outputkey =  dprintf('Alt');

function outputkey = getshiftkey
outputkey = dprintf('Shift');

function outputkey = getctrlkey
outputkey = dprintf('Ctrl');

function outputkey = getminuskey
outputkey =  dprintf('minus');

function outputkey = getpluskey
outputkey =  dprintf('plus');

function outputkey = getuparrow
outputkey =  dprintf('Up');

function outputkey = getdownarrow
outputkey =  dprintf('Down');

function outputkey = getrightarrow
outputkey =  dprintf('Right');

function outputkey = getleftarrow
outputkey =  dprintf('Left');

function outputkey = getmousewheel
outputkey = dprintf('Mouse wheel');