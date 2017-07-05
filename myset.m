function myset(handlevec,varargin)
%MYSET Same as SET, but removes NaN's in the handles.

%Einar Heiberg

handlevec = handlevec(ishandle(handlevec));
if ~isempty(handlevec)
  set(handlevec,varargin{:});
end;