function ok = mymovetosegmentfolder(tempfilename,newname,docopy)
%ok = MYMOVETOSEGMENTFOLDER(tempfilename,newname,[docopy])
%
%Moves the file tempfilename to the Segment path (where segment is
%installed) and gives it a new name. It checks first for write access 
%and if not then it ask for admin rights. If docopy is specified then
%a copy operation is performed instead of a move operation.
%
%tempfilename should be full path and filename to the file.
%
%newname should not contain any path just filename (can be relative path)
%
%See also MYSYSTEMADMIN, MYDELFROMSEGMENTFOLDER.

%Einar Heiberg

global DATA %#ok<GVMIS> 

%Changed from move to copy; as copy function assigns the same permissions
%as for SegmentFolder to the copied file; while move function moves the
%file togheter with permissions that are set for AppData folder, which
%might make the file unreadable for another users that have proper
%permission assigned to the SegmentFolder

if nargin<3
  docopy = false;
end

newfilename = [DATA.SegmentFolder filesep newname];

%just copy file 
commandstri = sprintf('copy "%s" "%s"',tempfilename,newfilename);
[failed, ~] = system(commandstri);

if ~failed
 ok = true;

 if ~docopy
   %delete temp file
   commandstri = sprintf('del "%s" ',tempfilename);
   [status,msg] = system(commandstri);
   if ~isequal(status,0)
     sprintf('Cannot delete temporary files. Error message: %s',msg);
   end
 end

 return
end

%ok that did not work, need to work harder...
tempbatfilename = [getpreferencespath filesep 'tempmove.bat'];
fid = fopen(tempbatfilename,'wt');

fprintf(fid,'copy "%s" "%s"\n',tempfilename,newfilename);  
fclose(fid);

%stri = sprintf('Elevate64.exe -wait4exit "%s"',tempbatfilename);
%stri = sprintf('myelevate.bat "%s"',tempbatfilename);
[status, result] = mysystemadmin(tempbatfilename);
%[status,result] = system(stri);
pause(1) % pause before checking if the file exists
if isequal(status,0)
  for loop = 1: 3
    %Ok check that it is really there
    if exist(newfilename,'file')
      ok = true;
      break
    else
      ok = false;
      if loop < 3
         pause(1) % give time for file to arrive
      end
    end
  end
else
  ok = false;
end

if ~ok
  disp(sprintf('Could not move file %s to %s, write protect? Message: %s',tempfilename,newfilename,result)); %#ok<DSPS>
end

%delete temp file(s)
if docopy
  commandstri = sprintf('del "%s"',tempbatfilename);
else
  commandstri = sprintf('del "%s" "%s"',tempbatfilename,tempfilename);
end

[status,msg] = system(commandstri);
if ~isequal(status,0)
  sprintf('Cannot delete temperory files. Error message: %s',msg);
end