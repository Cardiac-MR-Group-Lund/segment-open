function deltree(fileLocation)
% Works like deltree /y in good old Windows. 
% You can specify files, folders, or wildcards.
%
% Example: deltree('C:\Temp\Administrator\*.*');

if( ispc )
    executeString = sprintf('del /s /q %s', fileLocation);
    system(executeString);
    rmdir(fileLocation, 's');
else
    executeString = sprintf('rm -Rf %s', fileLocation);
    system(executeString);
end
