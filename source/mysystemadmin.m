function [status,result] = mysystemadmin(commandstri)
%Runs a command as administrator.
%Syntax just as system commmand. Note only works on Windows Platform.

%Example 1:
%----------
%[status,result] = mysystemadmin('dir');
%
%Example 2:
%----------
%command = 'nssm.exe';
%stri = sprintf('"%s%s%s" set SegmentSorterServer Description "This sorts incoming files"',DATA.SegmentFolder,filesep,command);
%[status,result] = mysystemadmin(stri);
%
%See also MYMOVETOSEGMENTFOLDER, MYDELFROMSEGMENTFOLDER

%Einar Heiberg

%stri = sprintf('Elevate64.exe -wait4exit %s',commandstri);
if isfile(commandstri)
    stri = sprintf('powershell start-process ''%s'' -verb runas',commandstri); %Powershell syntax: '[filepath]'
else
    stri = sprintf('powershell start-process %s -workingdirectory ''%s'' -verb runas',commandstri,pwd);
end
disp(sprintf('mysystemadmin: %s',stri)) %#ok<DSPS>

[status,result] = system(stri);



% start-process .\nssm.exe 'install SegmentSync "C:\Medviso Development\Source\branches\develop\segmentsync.bat"' -verb runas
