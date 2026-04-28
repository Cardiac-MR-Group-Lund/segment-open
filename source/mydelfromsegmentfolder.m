function ok = mydelfromsegmentfolder(filename)
%ok = MYDELFROMSEGMENTFOLDER(filename)
%
%Deletes the file filename from Segment folder 
%(where segment is installed). It checks first 
%for write access  and if not then it ask for 
%admin rights.
%
%ok = mydelfromsegmentfolder(filename) where 
%filename should not contain any path just filename 
%(or relative path).
%
%See also MYMOVETOSEGMENTFOLDER, MYSYSTEMADMIN.

%Einar Heiberg

global DATA %#ok<GVMIS> 

fullfilename = [DATA.SegmentFolder filesep filename];

%Check if exist
if ~exist(fullfilename,'file')
  disp(sprintf('mydelfromsegmentfolder: %s does not exist, no work to do.',fullfilename)); %#ok<DSPS>
  ok = true;
  return
end

%try with just delete
try
  delete(fullfilename);
catch
end

if ~exist(fullfilename,'file')
  ok = true;
  return
end

%ok that did not work, need to work harder...
tempbatfilename = [getpreferencespath filesep 'tempdel.bat'];
fid = fopen(tempbatfilename,'wt');
fprintf(fid,'del "%s"\n',fullfilename);
fclose(fid);

%Call the new .bat-file
% stri = sprintf('Elevate64.exe -wait4exit "%s"',tempbatfilename);
% [~,result] = system(stri);
[~, result] = mysystemadmin(tempbatfilename);

%Ok check that it is really not there
if ~exist(fullfilename,'file')
  ok = true;
else
  ok = false;
end

if ~ok
  disp(sprintf('Could not delete file %s. Write protect? Message: %s',fullfilename,result)); %#ok<DSPS>
else
  disp(sprintf('Deleted %s.',fullfilename)); %#ok<DSPS>
end