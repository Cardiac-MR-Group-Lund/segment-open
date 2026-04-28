function success = anonymizedicomfile(inputfilename, newpatientname,createmode)
% Anonymizes DICOM files
% 2021-09 initial implementation by Jelena Bock
% implemetation is based on Matlab's dicomanon.m

% embed this function in try/catch to get error messages

if nargin == 0
  error('Segment Pseudonymization ERROR: no input arguments provided')
end
% check if the file exist
if ~exist(inputfilename,'file')
  error('Segment Pseudonymization ERROR: File %s does not exist.',inputfilename)
end
% check if newpatientname was provided
if nargin < 2 || isempty(newpatientname)
  newpatientname = 'Hidden';
end

if nargin < 3 || ~any(strcmpi(createmode,{'create','copy'}))
  createmode = 'create'; % can be 'copy', however this mode does not verify parameters
end
createmode = lower(createmode);

% check if it is a string
if ~ischar(newpatientname)
  if isstring
    newpatientname = char(newpatientname);
  elseif isnumeric(newpatientname)
    newpatientname = num2str(newpatientname);
  else
    error('Segment Pseudonymization ERROR: Invalid new name')
  end
end

success = false;

% fixed settings (at the moment of initial implemetation)
writeprivateflag = true;

% Get original  info data
%  'UseVRHeuristic' is true by default? Read noncompliant DICOM files that switch VR modes incorrectly
try
  metadata = dicominfo(inputfilename);
catch me
  error('Segment Pseudonymization ERROR: dicominfo failed in %s',inputfilename)  
end

[X, map] = dicomread(metadata);

if isempty(X) %no pixel data, must call dcmodify
  global DATA %#ok<TLEV,GVMIS> #only call global variable if needed
  
  dcmstr = sprintf('%s\\dcmodify.exe -imt -nb -ie -m "(0010,0020)=%s" -m "(0010,0010)=%s" -m "(0008,0080)= " -m "(0010,1040)= " -m "(0008,0081)= " -m "(0008,0081)= "-m "(0010,0030)= " -m "(0010,1000)= " -m "(0010,2160)= " -m "(0010,2180)= " -m "(0010,21B0)= " -m "(0010,4000)= " "%s"',  DATA.SegmentFolder, newpatientname, newpatientname, inputfilename);
  system(dcmstr)
  success = true;
  return
end
dictionary = dicomdict('get_current');

metadata = changeAttr(metadata, '0008', '0080', '',  dictionary); %InstitutionName
metadata = removeAttr(metadata, '0008', '0081',  dictionary); %InstitutionAddress
metadata = changeAttr(metadata, '0010', '0010', newpatientname,  dictionary); %Patient Name
metadata = changeAttr(metadata, '0010', '0020', newpatientname,  dictionary); %PatientID
metadata = changeAttr(metadata, '0010', '0030', '',  dictionary); %Date of Birth
metadata = changeAttr(metadata, '0010', '1040', '',  dictionary); %Patient's Address
metadata = removeAttr(metadata, '0010', '1000',  dictionary); %OtherPatientID's
metadata = removeAttr(metadata, '0010', '2160',  dictionary); %Ethnic Group
metadata = removeAttr(metadata, '0010', '2180',  dictionary); %Occupation
metadata = removeAttr(metadata, '0010', '21B0',  dictionary); %AdditionalPatientHistory
metadata = removeAttr(metadata, '0010', '4000',  dictionary); %PatientComments

% Write the new data file.
try
  if (~isempty(map))
    sts = dicomwrite(X, map, inputfilename, metadata,'createmode', createmode ,'dictionary', dictionary, 'WritePrivate', writeprivateflag, 'UseMetadataBitDepths', true);
  else
    sts = dicomwrite(X, inputfilename, metadata, 'createmode', createmode, 'dictionary', dictionary, 'WritePrivate', writeprivateflag, 'UseMetadataBitDepths', true);
  end
catch me
  mydispexception(me) 
  error('Segment Pseudonymization ERROR: dicomwrite failed in %s',inputfilename)  
end
% check status after dicom write
if isempty(sts) || (isempty(sts.BadAttribute) && isempty(sts.MissingCondition) && isempty(sts.MissingData))
  success = true;
else
  disp(sts); 
  %   'BadAttribute'
  % The attribute's internal description is bad. It might be missing from the data dictionary or have incorrect data in its description.
  % 'MissingCondition'
  % The attribute is conditional but no condition has been provided for when to use it.
  % 'MissingData'
  % No data was provided for an attribute that must appear in the file.
  % 'SuspectAttribute'
  % Data in the attribute does not match a list of enumerated values in the DICOM specification.
end

%------------
function metadata = changeAttr(metadata, group, element, newValue,  dictionary)
%------------
%CHANGEATTR  Update an attribute's value.

name = getattrname(group, element, dictionary);


if ((~isempty(name)) && (isfield(metadata, name)))
    metadata.(name) = newValue;
end


%------------
function metadata = removeAttr(metadata, group, element, dictionary)
%------------
%REMOVEATTR  Remove an attribute.

name = getattrname(group, element, dictionary);


if ((~isempty(name)) && (isfield(metadata, name)))
    metadata = rmfield(metadata, name);
end

%------------
function name = getattrname(groupStr, elementStr, dictionary)
%------------
attr = getattrfromdicom(groupStr, elementStr, dictionary);

if (isempty(attr))
    name = '';
else
    name = attr.Name;
end

%------------
function attr = getattrfromdicom(group, element, dictionary)
%------------


MAX_GROUP = 65535;   % 0xFFFF
MAX_ELEMENT = 65535;  % 0xFFFF

%
% Load the data dictionary.
%

persistent tags values prev_dictionary;
mlock;

% Load dictionary for the first time or if dictionary has changed.
if ((isempty(values)) || (~isequal(prev_dictionary, dictionary)))
    
    [tags, values] = images.internal.dicom.loadDictionary(dictionary);
    prev_dictionary = dictionary;
    
end

%
% Convert hex strings to decimals.
%

if (ischar(group))
    group = sscanf(group, '%x');
end

if (ischar(element))
    element = sscanf(element, '%x');
end

if (group > MAX_GROUP)
    error(message('images:dicom_dict_lookup:groupOutOfRange', sprintf( '%x', group ), sprintf( '(%x,%04x)', group, element )))
end


if (element > MAX_ELEMENT)
    error(message('images:dicom_dict_lookup:elementOutOfRange', sprintf( '%x', element ), sprintf( '(%04x,%x)', group, element )))
end

%
% Look up the attribute.
%

% Group and Element values in the published data dictionary are 0-based.
index = tags((group + 1), (element + 1));

if (index == 0)
    attr = struct([]);
else
    attr = values(index);
end
