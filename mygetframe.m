function z = mygetframe(varargin)
%Same as getframe, but with error checking for multiple screens
%Also makes sure that the frame is visible on screen
%Einar Heiberg

if nargin > 0
  h = varargin{1};
  set(h,'Units','normalized');
  pos = get(h,'Position');
  if any(pos(1:2) < 0) || any(pos(1:2)+pos(3:4) > 1)
    set(h,'Position',[0.1 0.1 0.8 0.7]);
  end
end
try
  z = getframe(varargin{:});
catch %#ok<CTCH>
  myfailed('Could not perform screen capture. Could be depending on two monitors. Please move main window to primary monitor and retry.');
end;