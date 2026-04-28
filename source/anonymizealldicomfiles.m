function varargout = anonymizealldicomfiles(pathname, newpatientname,logfilename,actpatientnum,allpatientnum)
% Pseudonymizes all DICOM files found under pathname
% uses Medvisos anonymizedicomfile

% 2021-09 initial implementation by Jelena Bock

global DATA
	 
% Default arguments
if nargin == 0
  if isempty(DATA)
    startpath = pwd;
  else
    startpath = DATA.Pref.datapath;
  end
  pathname = myuigetdir(startpath,'Select directory to pseudonymize DICOM files in. Note: recursive behaviour.');
  if isequal(pathname,0)
    error('Segment Pseudonymization ERROR: Aborted by user');    
  end
end
if ~exist(pathname,'dir')
  % should this go into log file?
  return
end

if nargin < 2
  newpatientname = myinputdlg({'Enter new patient name'},'Patient Name',1,{'Hidden'});
  if isempty(newpatientname)
    error('Segment Pseudonymization ERROR: Invalid new name')    
  else
    newpatientname = newpatientname{1};
  end
end

% check if newpatientname is a string
if ~ischar(newpatientname)
  if isstring
    newpatientname = char(newpatientname);
  elseif isnumeric(newpatientname)
    newpatientname = num2str(newpatientname);
  else
    % log file?
    error('Segment Pseudonymization ERROR: Invalid new name')    
  end
end
newpatientname = regexprep(newpatientname, '[^\w'']', ''); %delete special characters
if length(newpatientname)>64 %can't be longer than 64
  newpatientname = newpatientname(1:64);
end

if nargin < 3
  filepath = getpreferencespath;
  logfilename = [filepath filesep sprintf('pseudonymization_errorlog_%s.log',datestr(now,'yyyymmddHHMMSS'))];
else
  [folderpath,~,~] = fileparts(logfilename);
  if ~exist(folderpath,'dir')
    % folder to write log file into does not exist
    error('Segment Pseudonymization ERROR: Invalid log file location %s',logfilename)    
  end
end

waitbarmsg = sprintf('Pseudonymizing patient. Please wait');
if nargin == 5
  if isnumeric(actpatientnum) && isnumeric(allpatientnum)
    waitbarmsg = sprintf('Pseudonymizing patient %d of %d. Please wait',actpatientnum,allpatientnum);
  end
end
% start to write into log file
writelog('-------------------',logfilename);

%Find files to pseudonymize
[filesinfolder,~,~] = createtree(pathname);

if all(cellfun('isempty',filesinfolder))
  myfailed('No files to pseudonymize');
  return;
end

disp('DICOM pseudonymization starting.');

numfiles = length(filesinfolder);


% loop over all found files
fileindex = 1;
createmode = 'create';
h = mywaitbarstart(numfiles+1,waitbarmsg,1,DATA.GUI.Anonymization);
d = [];

for loop = 1:numfiles
  h = mywaitbarupdate(h);
  currentfile = filesinfolder{loop};
  [~,filename,ext] = fileparts(currentfile);
  if isequal(filename,'.') || isequal(filename,'..')
     continue
  end
  
  

  if isequal(ext,'.cache')
    if ~isequal(filename,'folders') && ~isequal(filename,'thumbs')
      % delete cache, but not for folders.cache and thumbs.cache
      delete(currentfile);
    end
    continue
  end

  if isequal(ext,'.txt')
    continue
  end  

  if isempty(d)
    d = dicominfo(currentfile);
  end

  try
    successanon = anonymizedicomfile(currentfile, newpatientname,createmode);
  catch me
    successanon = false;
    writelog(me.message);
  end
  if strcmpi(createmode,'create') && ~successanon
    % try one more time with copy as crate mode
    createmode = 'copy';
    try
      successanon = anonymizedicomfile(currentfile, newpatientname,createmode);
    catch me
      writelog(sprintf('File %s failed with %s',currentfile,me.message))
    end
  end
  % compare patient details
  try
    verifiedbyname = comparepatientnames(currentfile,newpatientname);
  catch me
    writelog(sprintf('File %s failed with %s',currentfile,me.message))
  end

  if verifiedbyname
    if ~successanon
      % something went wrong at the pseudonymization, but it has new patient
      % name
      message = sprintf('File %s might failed in DICOM attribute writing',currentfile);
      writelog(message)
    end
    % rename file
    [wasrenamed,~] = renamefile(currentfile,fileindex);
    fileindex = fileindex+1;
    if ~wasrenamed
      % file was not renamed properly
      % write into error log
      message = sprintf('File %s NOT renamed',currentfile);
      writelog(message)
    end      
  else
    % not succesful pseudonymization
    % write into error log
    message = sprintf('File %s NOT pseudonymized',currentfile);
    writelog(message)
  end
end

name = '';
firstname = '';
date = '';
id = '';
if isfield(d,'PatientName')
  if isfield(d.PatientName, 'FamilyName')
  name = d.PatientName.FamilyName;
  end
end

if isfield(d,'PatientName')
  if isfield(d.PatientName, 'GivenName')
  firstname = d.PatientName.GivenName;
  end
end

if isfield(d,'PatientID')
  id = d.PatientID;
end

if isfield(d,'PatientBirthDate')
  date = d.PatientBirthDate;
end

output{1} = [name ' ' firstname] ;
output{2} = id;
output{3} = date;

mywaitbarclose(h);

if ~nargout==0
  varargout=output;
end

%----------------------------------------------------
function issuccess = comparepatientnames(filename,newpatientname)
%----------------------------------------------------
% compare if file after the pseudonymization has gotten the right patient name
issuccess = false;
info = fastdicominfo(filename);   
if ischar(info)
  return
end

if ~isequal(deblank(info.PatientName),deblank(newpatientname))
  return
else
  issuccess = true;
end

%-----------------------
function [success,msg] = renamefile(filepath,fileindex)
%-----------------------
% Rename file
success = false;
msg = 'rename failed';
if nargin < 2
  return
end
[pathstr,~,ext] = fileparts(filepath);
if all(isstrprop(ext(2:end),'digit'))
  % whole extension is a number, then skip extension
  ext = '.dcm';
end
newfilename = [datestr(now,'yyyymmddHHMMSS') sprintf('-%05d',fileindex),ext];

newfullname = fullfile(pathstr,newfilename);
try
  [success,msg] = movefile(filepath,newfullname,'f');
  if ~success
    [success,msg] = renamefile_helper(filepath,newfilename);
  end
catch
  % first attempt to rename failed, so try to rename using system command REN
  [success,msg] = renamefile_helper(filepath,newfilename); 
end

%-----------------------
function [success,msg] = renamefile_helper(filepath,newfilename)
%-----------------------
% fucntion to rename file using command line
success = false;
msg = 'rename failed';
try
  cmdstr = sprintf('REN "%s"="%s"',filepath,newfilename);
  [status,msg] = system(cmdstr);
  if status == 0
    success = true;
  else
    disp(msg);
  end
catch me
  mydispexception(me)
end