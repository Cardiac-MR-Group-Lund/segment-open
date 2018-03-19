function varargout = existfunctions(varargin)
% EXISTFUNCTIONS
% Functions for checking existence of segmentation

% Moved out from segment_main by Nisse Lundahl

%Invoke subfunction
macro_helper(varargin{:}); %future macro recording use
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
end;

if isempty(SET(no).EndoX) || all(isnan(SET(no).EndoX(:)))
  y = false;
else
  %Tested enough, make true
  y = true;
end;

%------------------------
function y = existepi(no) 
%------------------------
%True if epicardium exist in some slices
global SET NO

if nargin<1
  no = NO;
end;

if isempty(SET(no).EpiX) || all(isnan(SET(no).EpiX(:)))
  y = false;
else
  %Tested enough, make true
  y = true;
end;


%--------------------------------------
function y = existendoinselected(no,t,s) %#ok<DEFNU>
%--------------------------------------
%True if endocardium exists in selected slices.
%If argument is specified, only look in timeframe t.
global DATA SET NO

if nargin<1
  no = NO;
end;

if isempty(SET(no).EndoX) || all(isnan(SET(no).EndoX(:)))
  y = false;
  return;
end;

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
end;

%------------------------------------
function y = existepiinselected(no,t,s) %#ok<DEFNU>
%------------------------------------
%True if epicardium exist in selected slices.
%If argument is specified, only look in timeframe t.
global DATA SET NO

if nargin<1
  no = NO;
end;

if isempty(SET(no).EpiX)
  y = false;
  return;
end;

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
end;

%-----------------------------------
function y = existendoinslices(no,t) %#ok<DEFNU>
%-----------------------------------
%True for the slices where endocardium exists.
%If argument is specified, only look in timeframe t, oterwise look if there
%is endocardium in any time frame
global SET NO

if nargin<1
  no = NO;
end;

if isempty(SET(no).EndoX) || all(isnan(SET(no).EndoX(:)))
  y = false(SET(no).ZSize,1);
  return;
end;

if nargin == 2
  y = squeeze(~isnan(SET(no).EndoX(1,t,:)));
  return
end

y = squeeze(any(~isnan(SET(no).EndoX(1,:,:))));


%-----------------------------------
function y = existepiinslices(no,t) %#ok<DEFNU>
%-----------------------------------
%True for the slices where endocardium exists.
%If argument is specified, only look in timeframe t, oterwise look if there
%is endocardium in any time frame
global SET NO

if nargin<1
  no = NO;
end;

if isempty(SET(no).EpiX) || all(isnan(SET(no).EpiX(:)))
  y = false(SET(no).ZSize,1);
  return;
end;

if nargin == 2
  y = squeeze(~isnan(SET(no).EpiX(1,t,:)));
  return
end

y = squeeze(any(~isnan(SET(no).EpiX(1,:,:))));


%-------------------------------
function y = existrvendoinselected(no) %#ok<DEFNU>
%-------------------------------
%True if RV exist in selected slices.
global DATA SET NO

if nargin<1
  no = NO;
end;

if isempty(SET(no).RVEndoX)
  y = false;
  return;
end;

% if strcmp(DATA.ProgramName,'Segment CMR')
%   DATA.ThisFrameOnly = false;
% end

if DATA.ThisFrameOnly
  y = not(anyall(isnan(SET(no).RVEndoX(...
    1,SET(no).CurrentTimeFrame,SET(no).StartSlice:SET(no).EndSlice))));
else
  y = not(anyall(isnan(SET(no).RVEndoX(...
    1,:,SET(no).StartSlice:SET(no).EndSlice))));
end;

%-------------------------------
function y = existendoonlyinedes(no) %#ok<DEFNU>
%-------------------------------
%True if endocardium exists in ED and ES slices.
global SET NO

if nargin<1
  no = NO;
end;

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
function y = existepionlyinedes(no) %#ok<DEFNU>
%----------------------------------
%True if epicardium exists in ED and ES slices.
global SET NO

if nargin<1
  no = NO;
end;

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
function y = existrvendoonlyinedes(no) %#ok<DEFNU>
%-------------------------------------
%True if RV endocardium exists in ED and ES slices.
global SET NO

if nargin<1
  no = NO;
end;

if isempty(SET(no).RVEndoX) || all(isnan(SET(no).RVEndoX(:)))
  y = false;
  return;
end

slices=SET(no).StartSlice:SET(no).EndSlice;
esframe=SET(no).EST;
edframe=SET(no).EDT;

notedesframes=true(1,SET(no).TSize);
notedesframes(edframe)=0;
notedesframes(esframe)=0;

y =anyall(isnan(SET(no).RVEndoX(1,notedesframes,slices))) && ...
  not(anyall(isnan(SET(no).RVEndoX(1,esframe,slices)))) &&...
  not(anyall(isnan(SET(no).RVEndoX(1,edframe,slices))));

%------------------------------------
function y = existrvepionlyinedes(no) %#ok<DEFNU>
%------------------------------------
%True if RV epicardium exists in ED and ES slices.
global SET NO

if nargin<1
  no = NO;
end;

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

