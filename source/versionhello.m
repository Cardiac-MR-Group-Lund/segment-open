function versionhello(fcn, varargin)
%GUI for displaying message when starting a Segment session
feval(fcn,varargin{:})

%------------------------
function init(stri,title)
%------------------------
gui = openfig('versionhello.fig');
handles = guihandles(gui);
set(handles.figure1,'Name',title);
set(handles.titletext,'String',stri{1});
if length(stri) == 9
  set(handles.infotext,'String',stri{2});
  set(handles.refertext,'enable','off');
  set(handles.clinicaltext,'enable','off'); 
  set(handles.howtoreferpushbutton,'enable','off');
  set(handles.getsegmentcmrpushbutton,'enable','off');
  set(handles.getsegmentctpushbutton,'enable','off');
  set(handles.refertext,'visible','off');
  set(handles.clinicaltext,'visible','off'); 
  set(handles.howtoreferpushbutton,'visible','off');
  set(handles.getsegmentcmrpushbutton,'visible','off');
  set(handles.getsegmentctpushbutton,'visible','off'); 
else
  set(handles.refertext,'String',stri{2});
  set(handles.clinicaltext,'String',stri{3});
  set(handles.infotext,'enable','off'); 
  set(handles.infotext,'visible','off'); 
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

['We are deeply committed to continue to improve Segment, and keep ' ...
  'Segment freely available for research and educational purposes. ' ...
  'Please support us in this commitment by recommending the commercial ' ...
  'version Segment CMR to your clinical colleagues. Furthermore, there ' ...
  'is no other CMR software platform aside from Segment that allows ' ...
  'users to implement customized additions to clinical post-processing ' ...
  'or reporting software. Notably, in future procurement of MRI ' ...
  'scanners, you may want to consider adding Segment post-processing ' ...
  'software as a part of your purchasing requirement.'];

%-------------------
function ok_Callback
%-------------------
closereq;

%---------------------------
function howtorefer_Callback
%---------------------------
mybrowser('http://medviso.com/research/how-to-refer/');
closereq;

%------------------------------
function getsegmentcmr_Callback
%------------------------------
mybrowser('http://medviso.com/products/cmr/');
closereq;

%------------------------------
function getsegmentct_Callback
%------------------------------
mybrowser('http://medviso.com/products/ct/');
closereq;
