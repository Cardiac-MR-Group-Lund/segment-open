function pacstransferloop(~,varargin)
%-----------------------------
%This functions loops and updates the PACS transfer GUI until called upon
%with no input arguments. Helper function to PACS module.

%Einar Heiberg
global DATA

%Get handle
pacstransfergui = DATA.GUI.PacsTransfer;

if nargin>1  
  %Many arguments if called from receive => stop this and exits. This is called by the send
  %command to the server. When a response is received then
  pacstransfergui.TakeNext = true;
  
  %Add here varargin{3} to pacslog with addtopacslog in pacs.m EH:
  pacs('addtopacslog',['SLAVE: ' varargin{3}]);
else
  pacstransfergui.TakeNext = false;  
end
  
try
  while ~pacstransfergui.TakeNext;
    
    %Update the number of files transfered in this series graphically.
    f = dir(pacstransfergui.CurrentFolder);
    nfiles = sum(~cat(1,f(:).isdir));
    set(pacstransfergui.handles.filesinthisseriestext,'String',...
      sprintf('Files Transferred in this Serie:%04d',nfiles));
    pause(0.1);    
  end;
catch me
  errormsg = 'Error when updating PACS Transfer progress.';
  disp(errormsg);
  pacs('addtopacslog',['SEGMENT: ' errormsg]);
  mydispexception(me);
  pacstransfergui.TakeNext = false;
  pacstransfergui.Stop = true;
end;