% DOMATFILES Script that opens Segment Makematlab to load DICOM files and save 
% them to a .mat file

% Parse inputfile to get list of DICOMs
readok = true;
failmsg = '';
try
  fid = fopen('domatfilesinput.txt','r');
  if fid<0
    readok = false;
    failmsg = 'Could not open input file';
  else
    outputpath = fgetl(fid);
    transferid = fgetl(fid);
    nbrofseries = str2double(fgetl(fid));
    dicomcell = cell(nbrofseries,1);
    for i = 1:nbrofseries
      nbroffiles = str2double(fgetl(fid));
      subcell = cell(nbroffiles,1);
      for j = 1:nbroffiles
        subcell{j} = fgetl(fid);
      end
      dicomcell{i} = subcell;
    end
  end;
catch me
  readok = false;
  failmsg = 'Could not read input file';
end
try
  fclose(fid);
catch me
  readok = false;
  failmsg = 'Could close input file';
end;

% Load files in Segment CMR and save to a .mat file
global DATA SET
if readok
  try
    transfer.domatfilesgo; %Start Segment Makematlab
    DATA.GUISettings.AskToExitProgram = false;
    DATA.Pref.DoNotAsk = true;
    
    %Loop over stacks to load
    for i = 1:nbrofseries
      if ~isempty(dicomcell{i})
        openfile('loadfiles',dicomcell{i},false,[]);
      end
    end
    DATA.Silent = true;
    
    %autocropall; %Automatically crop it
    
    [pathname,filename] = fileparts(outputpath);
    [SET.transferid] = deal(transferid); %Update transfer information = ID stamp the source
    
    %Store it
    fail = filemenu('saveallas_helper',pathname,filename,false);
    if fail
      error('SEGMENT:ERROR','Could not save .mat file');
    end
  catch me
    mydispexception(me);
    failmsg = me.message;
    s = regexp(failmsg,'\n');
    for j = s
      failmsg(j) = ' ';
    end
  end
  filemenu('quit_Callback');
end

fid = fopen('domatfilesoutput.txt','w');
if isempty(failmsg)
  fprintf(fid,'SUCCESS\n');
else
  fprintf(fid,'FAILED\n');
  fprintf(fid,failmsg);
end
fclose(fid);