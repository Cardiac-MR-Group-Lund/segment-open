function mymaximize(fig)
%MYMAXIMIZE(HANDLE) Maximizes a figure or a mygui object 
%and aligns it to correct monitor.

%Einar Heiberg

units=get(fig,'units');
set(fig,'units','normalized');
pos=[0.05 0.05 0.95 0.85];
set(fig,'position',pos);
set(fig,'units',units);

% switch alignment
%   case 'left'
%     set(fig,'units','pixels','position',[...
%       DATA.PrimaryMonitorStart+5 ... %Left
%       40 ... %Bottom
%       DATA.PrimaryMonitorWidth-15 ... Width
%       DATA.PrimaryMonitorHeight-150]);
%   case 'right'
%     set(fig,'units','pixels','position',[...
%       DATA.SecondaryMonitorStart+5 ... %Left 
%       40 ... %Bottom
%       DATA.SecondaryMonitorWidth-15 ... Width
%       DATA.SecondaryMonitorHeight-150]);    
%   otherwise
%     myfailed('Unknown alignment specification.');    
% end;