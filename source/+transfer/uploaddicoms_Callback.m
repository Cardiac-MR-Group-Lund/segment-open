function uploaddicoms_Callback
%This function upload DICOM files to central image repository.

%Written by Jonatan Wulcan. Commented and modified by Einar Heiberg

global DATA

%In the previous version it asked from where to upload data. Now assume and
%ask.

%Find folder of CD-ROM
pathname = myuigetdir('D:','Select the CD-ROM drive.');

if isequal(pathname,0)
  myfailed('Aborted.');
  return;
end;

%Check if folder containing DICOMS
f = dir(pathname);
numdicomfolders = 0;
dicomfolder = pathname;
for loop = 1:length(f)
  if f(loop).isdir
    if ~isempty(strfind(lower(f(loop).name),'dicom'))
      numdicomfolders = numdicomfolders+1;
      dicomfolder = f(loop).name;
    end;
  end;
end;

%If found only one then 
if numdicomfolders == 1
  if yesno('Found one folder named ''dicom''. Do you want to restrict to this folder? (recommended)');
    pathname = [pathname filesep dicomfolder];
  else
    disp('Taking whole CD.');
  end;
end;

try
  wb = @(n, msg) transfer.waitbar(n, msg);
  [~, user, pw, url, pattern] = transfer.readcredentials();
catch me
  disp('Could not read credential file.');
  mydispexception(me);
  transfer.failmessage;
  return; %Since fatal
end;

try
  %Create Server object
  s = transfer.server(user, pw, url);
catch me
  disp('Could not create server object.');
  mydispexception(me);
  transfer.failmessage;
  return; %Since fatal
end;

try
  %Prepare the study (sort and anonymize)
  study = transfer.study(pathname, pattern, wb);
catch me
  disp('Could not prepare the study (sort & anonymize)');
  mydispexception(me);
  transfer.failmessage;
  return; %Since fatal
end;

try
  %Send the study
  disp('Starting to send studies.');
  transfer.progressbar(0.3,'Transfering studies.');
  didfail = transfer.sendstudy(s, transfer.path2study(study.getpath), ...
    study.getname, wb);
catch me
  mydispexception(me);
  transfer.failmessage;
  return;
end;

try
  clear study  
catch me
  disp('Could not clear study.');
  mydispexception(me);
  transfer.failmessage;
end

% Display message
old_dna = DATA.Pref.DoNotAsk;
DATA.Pref.DoNotAsk = false; %Ensure that message box is displayed.
DATA.Buffer.KeyStroke = {}; %Ensure that message box is displayed.

if ~didfail
  mymsgbox('Upload completed successfully!');
else
  mymsgbox('Files uploaded. Errors were detected.!');
  transfer.failmessage;
end;
DATA.Pref.DoNotAsk = old_dna;

