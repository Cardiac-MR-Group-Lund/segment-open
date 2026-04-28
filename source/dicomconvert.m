function dicomconvert(infile,outfile,dinfo,usematlabforjpeg,usematlabforall)
%Reads the infile and writes the outfile with a new transfersyntax (Explicit VT little endian). 
%
%If called with only one file, then convert filename with the same
%outfilename.

%Einar Heiberg

%#ok<*GVMIS>

outsyntax = '1.2.840.10008.1.2.1';

usedcmtk = true;  %true for a start

jpeglist=[
{'1.2.840.10008.1.2.4.57'};...
{'1.2.840.10008.1.2.4.70'};...
{'1.2.840.10008.1.2.4.90'};...
{'1.2.840.10008.1.2.4.50'};...
{'1.2.840.10008.1.2.4.51'};...
{'1.2.840.10008.1.2.4.70'};...
{'1.2.840.10008.1.2.4.91'};...
{'1.2.840.10008.1.2.4.91'};...
{'1.2.840.10008.1.2.4.92'}];

jpegdcmtk = [ ...
  {'1.2.840.10008.1.2.4.50'};
  {'1.2.840.10008.1.2.4.51'};
  {'1.2.840.10008.1.2.4.53'};
  {'1.2.840.10008.1.2.4.55'};
  {'1.2.840.10008.1.2.4.57'};
  {'1.2.840.10008.1.2.4.70'};
  ];

jpegexcluded = [...            %These are JPEG-LS, not supported
{'1.2.840.10008.1.2.4.80'};... %JPEG-LS Lossless
{'1.2.840.10008.1.2.4.81'}];    %JPEG-LS Near-Lossless (Lossy)



warnId=[]; %#ok<NASGU>

if nargin<2
  outfile = infile;
end

if ~exist(infile,'file')
  error('Input file does not exist.');
end

if (nargin<3) || isempty(dinfo)
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

if (nargin < 4) || isempty(usematlabforjpeg)
  usematlabforjpeg = true;
end

%if isequal(dinfo.TransferSyntaxUID,'jpeg') || isequal(dinfo.TransferSyntaxUID,'1.2.840.10008.1.2.4.90')|| isequal(dinfo.TransferSyntaxUID,'1.2.840.10008.1.2.4.70')
if isequal(dinfo.TransferSyntaxUID,'jpeg') || any(ismember(dinfo.TransferSyntaxUID,jpeglist))
  usedcmtk = false;
elseif isequal(dinfo.TransferSyntaxUID,'explicitbigendian') || isequal(dinfo.TransferSyntaxUID,'1.2.840.10008.1.2.2')
  usedcmtk = true;
end

%Check if we should force it to Matlab
if (nargin > 4)
  usedcmtk = ~usematlabforall;
end

if usedcmtk
  %Use DCMTK to convert
  stri = sprintf('dcmconv%s +te "%s" "%s"',myexecutableext,infile,outfile); %+te write explicit VR little endian.
  [s,w] = system(stri);
  if ~isequal(s,0)
    error(w);
  end
else
  %Use Matlab code to convert.

  %---Read the image
  try
    lastwarn(''); % Clear last warning message
    
    dinfo = dicominfo(infile); %re-read dinfo with stable loader to get all correctly read.
    if ~usematlabforjpeg
      % check if supported transfer syntax in dcmtk
      isdcmdjpegsupported = contains(dinfo.TransferSyntaxUID,jpegdcmtk);
    else
      isdcmdjpegsupported = false;
    end
   
    [msg, warnId] =  lastwarn; %Check the warning status
    if ~isempty(warnId)
      warning('off',warnId);
      try
        deletedicomtag('2001,105f',infile) %Special fix to dr Paul Joy dataset; Philips JPEG
      catch
        sprintf('Unable to delete ''2001,105f'' tag. There might be another tag to remove : %s',msg)
      end
    end
    if usematlabforjpeg || ~isdcmdjpegsupported
      X = dicomread(dinfo);
    end
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
    MediaStorageSOPClass = deblank(string(dinfo.MediaStorageSOPClassUID));
    if isequal(MediaStorageSOPClass,'1.2.840.10008.5.1.4.1.1.4.1') %Enhanced DICOM
      if usematlabforjpeg || ~isdcmdjpegsupported
        dicomwrite(X,outfile, dinfo, 'CreateMode', 'Copy', 'MultiframeSingleFile', 'true');
      else
        stri = sprintf('dcmdjpeg%s "%s" "%s"',myexecutableext,infile,outfile);
        [s,w] = system(stri);
        if ~isequal(s,0)
          error(w);
        end
      end
    else
      if usematlabforjpeg 
        dicomwrite(X,outfile,dinfo, 'CreateMode', 'Copy'); %to assure that all main tags are copied (except Private Tags)
      elseif ~isdcmdjpegsupported
        %to assure that all tags together with Private tags are copied); not always stable
        dicomwrite(X,outfile,dinfo, 'CreateMode', 'Copy','WritePrivate', true);
      else
        stri = sprintf('dcmdjpeg%s "%s" "%s"',myexecutableext,infile,outfile);
        [s,w] = system(stri);
        if ~isequal(s,0)
          error(w);
        end
      end
    end
  catch me
    mydispexception(me);
    error('Could not write DICOM file to %s',outfile);   
  end

end

%---------------------------------
function deletedicomtag(tag, file)
%---------------------------------
%Function to delete tag from dicom file

global DATA

source = [DATA.SegmentFolder filesep 'dcmodify'];

stri = sprintf('"%s" -nb -e (%s) "%s"', source, tag, file);
[status,~] = system(stri);
if ~isequal(status,0)
  fprintf('Warning, error in deleting tag (%s) in %s.\n',tag, file);
end