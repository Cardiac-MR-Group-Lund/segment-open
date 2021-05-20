function versionhello(fcn, varargin)
%GUI for displaying message when starting a Segment session
feval(fcn,varargin{:})

%------------------------
function init(stri,title) %#ok<DEFNU>
%------------------------
global DATA

if isfield(DATA.Pref,'NoPopUp')
  if DATA.Pref.NoPopUp
    return
  end
end

gui = openfig('versionhello.fig');
setupicon(gui);
handles = guihandles(gui);
set(handles.figure1,'Name',title);
set(handles.titletext,'String',stri{1});
if length(stri) == 9
  set(handles.infotext,'String',stri{2});
  set(handles.refertext,'enable','off');
  set(handles.clinicaltext,'enable','off'); 
  set(handles.howtoreferpushbutton,'enable','off');
  set(handles.getsegmentcmrpushbuttont,'visible','off');
  set(handles.getsegmentctpushbutton,'visible','off');
  set(handles.getsegment3dppushbutton,'visible','off');
  set(handles.getstrainpushbutton,'visible','off');
  set(handles.segmentcmrtext,'visible','off');
  set(handles.segmentcttext,'visible','off');
  set(handles.segment3dptext,'visible','off');
  set(handles.refertext,'visible','off');
  set(handles.clinicaltext,'visible','off'); 
  set(handles.howtoreferpushbutton,'visible','off');
  set(handles.getsegmentcmrpushbuttont,'enable','off');
  set(handles.getsegmentctpushbutton,'enable','off');
  set(handles.getsegment3dppushbutton,'enable','off');
  set(handles.getstrainpushbutton,'enable','off');
else
  set(handles.refertext,'String',stri{2});
  set(handles.clinicaltext,'String',stri{3});
  set(handles.infotext,'enable','off'); 
  set(handles.infotext,'visible','off');
  
  try
    im = imread('segmentcmrimage.png');
    image(im,'Parent',handles.segmentcmraxes);
    axis(handles.segmentcmraxes,'equal');
    axis(handles.segmentcmraxes,'off');
    
    im = imread('segmentstrainimage.png');
    image(im,'Parent',handles.strainaxes);
    axis(handles.strainaxes,'equal');
    axis(handles.strainaxes,'off');
    
    im = imread('segmentctimage.png');
    image(im,'Parent',handles.segmentctaxes);
    axis(handles.segmentctaxes,'equal');
    axis(handles.segmentctaxes,'off');
    
    im = imread('segment3dpimage.png');
    image(im,'Parent',handles.segment3dpaxes);
    axis(handles.segment3dpaxes,'equal');
    axis(handles.segment3dpaxes,'off');
  catch %#ok<CTCH>
  end
end

if length(stri{2}) > 200
  set(handles.pushbutton5,'visible','off');
end

%logo
try
  im = imread('medviso.png');
  image(im,'Parent',handles.medvisoaxes);
  axis(handles.medvisoaxes,'equal');
  axis(handles.medvisoaxes,'off');
catch %#ok<CTCH>
end
set(gui,'visible','on');

%-------------------
function ok_Callback %#ok<DEFNU>
%-------------------
closereq;

%---------------------------
function howtorefer_Callback %#ok<DEFNU>
%---------------------------
mybrowser('http://medviso.com/research/how-to-refer/');
closereq;

%------------------------------
function getsegmentcmr_Callback %#ok<DEFNU>
%------------------------------
mybrowser('http://medviso.com/cmr/');
closereq;

%------------------------------
function getsegment3dp_Callback %#ok<DEFNU>
%------------------------------
mybrowser('http://medviso.com/segment-3dprint/');
closereq;

%------------------------------
function getstrain_Callback %#ok<DEFNU>
%------------------------------
mybrowser('http://medviso.com/strain/');
closereq;

%------------------------------
function getsegmentct_Callback %#ok<DEFNU>
%------------------------------
mybrowser('http://medviso.com/ct/');
closereq;
