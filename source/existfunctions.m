function varargout = existfunctions(varargin)
% EXISTFUNCTIONS
% Functions for checking existence of segmentation

% Moved out from segment_main by Nisse Lundahl

%#ok<*GVMIS> 

%Invoke subfunction
if (nargout)
  [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
else
  feval(varargin{:}); % FEVAL switchyard
end

%-------------------------
function y = existendo(no) 
%-------------------------
%True if endocardium exist in some slices
global SET NO 
if nargin<1
  no = NO;
end

if isempty(SET(no).EndoX) || all(isnan(SET(no).EndoX(:)))
  y = false;
else
  %Tested enough, make true
  y = true;
end

%------------------------
function y = existepi(no) 
%------------------------
%True if epicardium exist in some slices
global SET NO

if nargin<1
  no = NO;
end

if isempty(SET(no).EpiX) || all(isnan(SET(no).EpiX(:)))
  y = false;
else
  %Tested enough, make true
  y = true;
end

%-------------------------
function y = existrvendo(no) 
%-------------------------
%True if RV endocardium exist in some slices
global SET NO

if nargin < 1
  no = NO;
end

if isempty(SET(no).RVEndoX) || all(isnan(SET(no).RVEndoX(:)))
  y = false;
else
  %Tested enough, make true
  y = true;
end

%------------------------
function y = existrvepi(no) 
%------------------------
%True if RV epicardium exist in some slices
global SET NO

if nargin < 1
  no = NO;
end

if isempty(SET(no).RVEpiX) || all(isnan(SET(no).RVEpiX(:)))
  y = false;
else
  %Tested enough, make true
  y = true;
end

%------------------------
function bool = existatrium(no,heartpart)
%------------------------
%True if LA or RA contour exists in no
switch heartpart
  case 'LA'
    bool = existla(no);
  case 'RA'
    bool = existra(no);
  case 'any'
    bool = existla(no) || existra(no);
end

%------------------------
function bool = existlv(no)
%------------------------
%True if LV endo contour exists in no
bool = existventriculecontour(no,'LV');

%------------------------
function bool = existrv(no)
%------------------------
%True if RV endo contour exists in no
bool = existventriculecontour(no,'RV');

%------------------------
function bool = existlvedt(no)
%------------------------
%True if LV contour exists in no in EDT
bool = existventriculeedt(no,'LV');

%------------------------
function bool = existrvedt(no)
%------------------------
%True if RV contour exists in no in EDT
bool = existventriculeedt(no,'RV');

%------------------------
function bool = existla(no)
%------------------------
%True if LA contour exists in no
bool = existatriumcontour(no,'LA');

%------------------------
function bool = existra(no)
%------------------------
%True if LA contour exists in no
bool = existatriumcontour(no,'RA');

%------------------------
function bool = existlaedt(no)
%------------------------
%True if LA contour exists in no in EDT
bool = existatriumedt(no,'LA');

%------------------------
function bool = existraedt(no)
%------------------------
%True if RA contour exists in no in EDT
bool = existatriumedt(no,'RA');

%------------------------
function bool = existatriumedt(no,atrium)
%------------------------
%True if specified atrium contour exists in no
arguments
  no
  atrium {mustBeMember(atrium,{'LA','RA'})}
end
global SET

edt = SET(no).EDT;
bool = existatriumcontour(no,atrium) && ...
  ~all(isnan(SET(no).(atrium).X(:,edt)));

%--------------------------------
function bool = existraest(no)
%------------------------
%True if RA contour exists in no in EST
bool = existatriumest(no,'RA');

%--------------------------------
function bool = existlaest(no)
%------------------------
%True if RA contour exists in no in EST
bool = existatriumest(no,'LA');

%------------------------
function bool = existatriumest(no,atrium)
%------------------------
%True if specified atrium contour exists in no in EST
arguments
  no
  atrium {mustBeMember(atrium,{'LA','RA'})}
end
global SET

est = SET(no).EST;
bool = existatriumcontour(no,atrium) && ...
  ~all(isnan(SET(no).(atrium).X(:,est)));

%------------------------
function bool = existventriculeedt(no,ventricule)
%------------------------
%True if specified ventricule contour exists in no
arguments
  no
  ventricule {mustBeMember(ventricule,{'LV','RV'})}
end
global SET

switch ventricule
  case 'LV'
    xcontourfield = 'EndoX';
  case 'RV'
    xcontourfield = 'RVEndoX';
end
bool = existventriculecontour(no,ventricule);
if bool
  edt = SET(no).EDT;
  edtcontour = squeeze(SET(no).(xcontourfield)(:,edt,:));
  booledt = ~all(isnan(edtcontour(:)));
  bool = bool & booledt;
end

%------------------------
function bool = existatriumcontour(no,atrium)
%------------------------
%True if specified atrium contour exists in no
arguments
  no
  atrium {mustBeMember(atrium,{'LA','RA'})}
end
global SET
bool = ~isempty(SET(no).(atrium)) && ~isempty(SET(no).(atrium).X);

%------------------------
function bool = existventriculecontour(no,ventricule)
%------------------------
%True if specified atrium contour exists in no
arguments
  no
  ventricule {mustBeMember(ventricule,{'LV','RV'})}
end
global SET

switch ventricule
  case 'LV'
    xcontourfield = 'EndoX';
  case 'RV'
    xcontourfield = 'RVEndoX';
end

bool = ~isempty(SET(no).(xcontourfield)) && ...
  ~all(isnan(SET(no).(xcontourfield)(:)));

%--------------------------------------
function y = existendoinselected(no,t,s)
%--------------------------------------
%True if endocardium exists in selected slices.
%If argument is specified, only look in timeframe t.
global DATA SET NO

if nargin<1
  no = NO;
end

if isempty(SET(no).EndoX) || all(isnan(SET(no).EndoX(:)))
  y = false;
  return;
end

if nargin == 3
  y = not(anyall(isnan(SET(no).EndoX(...
    1,t,s))));
  return
end

if nargin == 2
  y = not(anyall(isnan(SET(no).EndoX(...
    1,t,SET(no).StartSlice:SET(no).EndSlice))));
  return
end

% if strcmp(DATA.ProgramName,'Segment CMR')
%   DATA.ThisFrameOnly = false;
% end

if DATA.ThisFrameOnly
  y = not(anyall(isnan(SET(no).EndoX(...
    1,SET(no).CurrentTimeFrame,SET(no).StartSlice:SET(no).EndSlice))));
else
  y = not(anyall(isnan(SET(no).EndoX(...
    1,:,SET(no).StartSlice:SET(no).EndSlice))));
end

%------------------------------------
function y = existepiinselected(no,t,s)
%------------------------------------
%True if epicardium exist in selected slices.
%If argument is specified, only look in timeframe t.
global DATA SET NO

if nargin<1
  no = NO;
end

if isempty(SET(no).EpiX)
  y = false;
  return;
end

if nargin == 2
  y = not(anyall(isnan(SET(no).EpiX(...
    1,t,SET(no).StartSlice:SET(no).EndSlice))));
  return
end

if nargin == 3
  y = not(anyall(isnan(SET(no).EpiX(...
    1,t,s))));
  return
end

% if strcmp(DATA.ProgramName,'Segment CMR')
%   DATA.ThisFrameOnly = false;
% end

if DATA.ThisFrameOnly
  y = not(anyall(isnan(SET(no).EpiX(...
    1,SET(no).CurrentTimeFrame,SET(no).StartSlice:SET(no).EndSlice))));
else
  y = not(anyall(isnan(SET(no).EpiX(...
    1,:,SET(no).StartSlice:SET(no).EndSlice))));
end

%-----------------------------------
function y = existendoinslices(no,t)
%-----------------------------------
%True for the slices where endocardium exists.
%If argument is specified, only look in timeframe t, oterwise look if there
%is endocardium in any time frame
global SET NO

if nargin<1
  no = NO;
end

if isempty(SET(no).EndoX) || all(isnan(SET(no).EndoX(:)))
  y = false(SET(no).ZSize,1);
  return;
end

if nargin == 2
  y = squeeze(~isnan(SET(no).EndoX(1,t,:)));
  return
end

y = squeeze(any(~isnan(SET(no).EndoX(1,:,:))));


%-----------------------------------
function y = existepiinslices(no,t)
%-----------------------------------
%True for the slices where endocardium exists.
%If argument is specified, only look in timeframe t, oterwise look if there
%is endocardium in any time frame
global SET NO

if nargin<1
  no = NO;
end

if isempty(SET(no).EpiX) || all(isnan(SET(no).EpiX(:)))
  y = false(SET(no).ZSize,1);
  return;
end

if nargin == 2
  y = squeeze(~isnan(SET(no).EpiX(1,t,:)));
  return
end

y = squeeze(any(~isnan(SET(no).EpiX(1,:,:))));


%-------------------------------
function y = existrvendoinselected(no)
%-------------------------------
%True if RV exist in selected slices.
global DATA SET NO

if nargin<1
  no = NO;
end

if isempty(SET(no).RVEndoX)
  y = false;
  return;
end

% if strcmp(DATA.ProgramName,'Segment CMR')
%   DATA.ThisFrameOnly = false;
% end

if DATA.ThisFrameOnly
  y = not(anyall(isnan(SET(no).RVEndoX(...
    1,SET(no).CurrentTimeFrame,SET(no).StartSlice:SET(no).EndSlice))));
else
  y = not(anyall(isnan(SET(no).RVEndoX(...
    1,:,SET(no).StartSlice:SET(no).EndSlice))));
end

%-------------------------------
function y = existendoonlyinedes(no)
%-------------------------------
%True if endocardium exists in ED and ES slices.
global SET NO

if nargin<1
  no = NO;
end

if isempty(SET(no).EndoX) || all(isnan(SET(no).EndoX(:)))
  y = false;
  return;
end

slices=SET(no).StartSlice:SET(no).EndSlice;
esframe=SET(no).EST;
edframe=SET(no).EDT;

notedesframes=true(1,SET(no).TSize);
notedesframes(edframe)=0;
notedesframes(esframe)=0;

y =anyall(isnan(SET(no).EndoX(1,notedesframes,slices))) && ...
  not(anyall(isnan(SET(no).EndoX(1,esframe,slices)))) &&...
  not(anyall(isnan(SET(no).EndoX(1,edframe,slices))));

%----------------------------------
function y = existepionlyinedes(no)
%----------------------------------
%True if epicardium exists in ED and ES slices.
global SET NO

if nargin<1
  no = NO;
end

if isempty(SET(no).EpiX) || all(isnan(SET(no).EpiX(:)))
  y = false;
  return;
end

slices=SET(no).StartSlice:SET(no).EndSlice;
esframe=SET(no).EST;
edframe=SET(no).EDT;

notedesframes=true(1,SET(no).TSize);
notedesframes(edframe)=0;
notedesframes(esframe)=0;

y =anyall(isnan(SET(no).EpiX(1,notedesframes,slices))) && ...
  not(anyall(isnan(SET(no).EpiX(1,esframe,slices)))) &&...
  not(anyall(isnan(SET(no).EpiX(1,edframe,slices))));

%-------------------------------------
function y = existrvendoonlyinedes(no)
%-------------------------------------
%True if RV endocardium exists in all slices in ED and ES time frames.
global SET NO

if nargin < 1
  no = NO;
end

if isempty(SET(no).RVEndoX) || all(isnan(SET(no).RVEndoX(:)))
  y = false;
  return;
end

slices=SET(no).StartSlice:SET(no).EndSlice;
esframe = SET(no).EST;
edframe = SET(no).EDT;

notedesframes = true(1,SET(no).TSize);
notedesframes(edframe) = 0;
notedesframes(esframe) = 0;

y = anyall(isnan(SET(no).RVEndoX(1,notedesframes,slices))) && ...
  not(anyall(isnan(SET(no).RVEndoX(1,esframe,slices)))) &&...
  not(anyall(isnan(SET(no).RVEndoX(1,edframe,slices))));


%-------------------------------------
function boolvalue = existrvendoinedes(no)
%-------------------------------------
%True if any RV endocardium exists in ED and ES time frames.
global SET NO

if nargin < 1
  no = NO;
end

if isempty(SET(no).RVEndoX) || all(isnan(SET(no).RVEndoX(:)))
  boolvalue = false;
  return;
end

esframe = SET(no).EST;
edframe = SET(no).EDT;

nslice = size(SET(no).RVEndoX,3);
slices = 1:nslice;

boolvalue = any(~isnan(SET(no).RVEndoX(1,esframe,slices))) && ...
  any(~isnan(SET(no).RVEndoX(1,edframe,slices)));

%------------------------------------
function y = existrvepionlyinedes(no)
%------------------------------------
%True if RV epicardium exists in ED and ES slices.
global SET NO

if nargin<1
  no = NO;
end

if isempty(SET(no).RVEpiX) || all(isnan(SET(no).RVEpiX(:)))
  y = false;
  return;
end

slices=SET(no).StartSlice:SET(no).EndSlice;
esframe=SET(no).EST;
edframe=SET(no).EDT;

notedesframes=true(1,SET(no).TSize);
notedesframes(edframe)=0;
notedesframes(esframe)=0;

y =anyall(isnan(SET(no).RVEpiX(1,notedesframes,slices))) && ...
  not(anyall(isnan(SET(no).RVEpiX(1,esframe,slices)))) &&...
  not(anyall(isnan(SET(no).RVEpiX(1,edframe,slices))));

%---------------------
function z = anyall(a)
%---------------------
%Equivalent to z = any(a(:)); 
z = any(a(:));

