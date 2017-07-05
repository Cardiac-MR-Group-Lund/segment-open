function dicomconvert(infile,outfile,syntax)
%DICOMCONVERT(INFILE,OUTFILE,[SYNTAX]);
%
%Reads the infile and writes the outfile with a new transfersyntax (default
%is 1.2.840.10008.1.2.1 (Explicit VT little endian). If transfer syntax is specified then Matlab 
%routines aree used, otherwise DCMTK routines are used (which is much faster).  
% 
%If called with only one file, then convert filename with the same 
%outfilename.

%Einar Heiberg

usedcmtk = false;
if nargin<3
  syntax = '1.2.840.10008.1.2.1';
end;
usedcmtk = true; 

if nargin<2
  outfile = infile;
end;

if ~exist(infile,'file')
  error('Input file does not exist.');
end;

%Read info from the dicom file
try
  dinfo = dicominfo(infile);
catch me
  mydispexception(me);
  error(sprintf('Could not read DICOM file %s.',infile)); %#ok<SPERR>
end;

%Read the image
% try
%   X = dicomread(dinfo);
% catch me
%   mydispexception(me);
%   error(sprintf('Could not read DICOM image in file %s',infile));   %#ok<SPERR>
% end;
if usedcmtk
  stri = sprintf('dcmconv%s +te "%s" "%s"',myexecutableext,infile,outfile); %+te write explicit VR little endian.
  [s,w] = system(stri);
  if ~isequal(s,0)
    error(w);
  end;
else
  %Read info from the dicom file
  try
    dinfo = dicominfo(infile);
  catch me
    mydispexception(me);
    error(sprintf('Could not read DICOM file %s.',infile)); %#ok<SPERR>
  end;
  
%   %Read the image
%   try
%     X = dicomread(dinfo);
%   catch me
%     mydispexception(me);
%     error(sprintf('Could not read DICOM image in file %s',infile));   %#ok<SPERR>
%   end;
%   
%Display information 
if isequal(infile,outfile)
  disp(sprintf('Converted encoding in file %s.',infile)); %#ok<DSPS>
else
  disp(sprintf('Converted encoding in files %s => %s.',infile,outfile)); %#ok<DSPS>
end

%Set new transfersyntax
dinfo.TransferSyntaxUID = syntax;
  
  %Write the dicom file
  try
    %dicomwrite(X,outfile,'WritePrivate',true,dinfo);
    dicomwrite(X,outfile,dinfo);
  catch me
    mydispexception(me);
    error(sprintf('Could not write DICOM file to %s',outfile));   %#ok<SPERR>
  end;
end;

% %Set new transfersyntax
% dinfo.TransferSyntaxUID = syntax;
% 
% %Write the dicom file
% try
%   %dicomwrite(X,outfile,'WritePrivate',true,dinfo);
%   dicomwrite(X,outfile,dinfo);
% catch me
%   mydispexception(me);
%   error(sprintf('Could not write DICOM file to %s',outfile));   %#ok<SPERR>
% end;


