function dicomconvert(infile,outfile,dinfo)
%DICOMCONVERT(INFILE,OUTFILE,[SYNTAX]);
%
%Reads the infile and writes the outfile with a new transfersyntax (default
%is 1.2.840.10008.1.2.1 (Explicit VT little endian). If transfer syntax is specified then Matlab
%routines aree used, otherwise DCMTK routines are used (which is much faster).
%
%If called with only one file, then convert filename with the same
%outfilename.

%Einar Heiberg

outsyntax = '1.2.840.10008.1.2.1';

usedcmtk = true;  %true for a start

jpeglist=[{'1.2.840.10008.1.2.4.57'};...
{'1.2.840.10008.1.2.4.70'};...
{'1.2.840.10008.1.2.4.90'};...
{'1.2.840.10008.1.2.4.50'};...
{'1.2.840.10008.1.2.4.51'};...
{'1.2.840.10008.1.2.4.70'};...
{'1.2.840.10008.1.2.4.80'};...
{'1.2.840.10008.1.2.4.81'};...
{'1.2.840.10008.1.2.4.91'};...
{'1.2.840.10008.1.2.4.91'};...
{'1.2.840.10008.1.2.4.92'}];

warnId=[]; %#ok<NASGU>

if nargin<2
  outfile = infile;
end

if ~exist(infile,'file')
  error('Input file does not exist.');
end

if nargin<3
  try
    dinfo = fastdicominfo(infile);
  catch
    try
      dinfo = dicominfo(infile);
    catch me
      mydispexception(me);
      error('Could not read DICOM file %s.',infile); 
    end
  end
end

%if isequal(dinfo.TransferSyntaxUID,'jpeg') || isequal(dinfo.TransferSyntaxUID,'1.2.840.10008.1.2.4.90')|| isequal(dinfo.TransferSyntaxUID,'1.2.840.10008.1.2.4.70')
if isequal(dinfo.TransferSyntaxUID,'jpeg') || any(ismember(dinfo.TransferSyntaxUID,jpeglist))
  usedcmtk = false;
elseif isequal(dinfo.TransferSyntaxUID,'explicitbigendian') || isequal(dinfo.TransferSyntaxUID,'1.2.840.10008.1.2.2')
  usedcmtk = true;
end

if usedcmtk
  stri = sprintf('dcmconv%s +te "%s" "%s"',myexecutableext,infile,outfile); %+te write explicit VR little endian.
  [s,w] = system(stri);
  if ~isequal(s,0)
    error(w);
  end
else
  %Read the image
  try
    lastwarn(''); % Clear last warning message
    
    dinfo = dicominfo(infile); %re-read dinfo with stable loader to get all correctly read.
   
    [msg, warnId] =  lastwarn; %Check the warning status
    if ~isempty(warnId)
      warning('off',warnId);
      try
        deletedicomtag('2001,105f',infile) %Special fix to dr Paul Joy dataset; Philips JPEG
      catch
        sprintf('Unable to delete ''2001,105f'' tag. There might be another tag to remove : %s',msg)
      end
    end
    X = dicomread(dinfo);
  catch me
    mydispexception(me);
    error('Could not read DICOM image in file %s',infile);   
  end
  
  %Display information
  if isequal(infile,outfile)
    disp(sprintf('Converted encoding in file %s.',infile)); %#ok<DSPS>
  else
    disp(sprintf('Converted encoding in files %s => %s.',infile,outfile)); %#ok<DSPS>
  end
  
  %Set new transfersyntax
  dinfo.TransferSyntaxUID = outsyntax;
  
  %Write the dicom file
  try
    %dicomwrite(X,outfile,'WritePrivate',true,dinfo);
    dicomwrite(X,outfile,dinfo);
  catch me
    mydispexception(me);
    error('Could not write DICOM file to %s',outfile);   
  end
end

%-----------------------------------------------------------------------------------------------------------
function deletedicomtag(tag, file)
%-----------------------------------------------------------------------------------------------------------
%Function to delete tag from dicom file
global DATA
source=[DATA.SegmentFolder filesep 'dcmodify'];

stri = sprintf('"%s" -nb -e (%s) "%s"', source, tag, file);
[status,~] = system(stri);
if ~isequal(status,0)
    fprintf('Warning, error in deleting tag (%s) in %s.\n',tag, file);
end