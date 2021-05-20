 function varargout = anonymization(varargin)
%---------------------------------------------
% ANONYMIZATION GUI for pseudonymization of dicom-files or mat-files as
%extensive, partial or corrective pseudonymization, for a single file/folder or for
%a whole working list.

% Written by Jane Tufvesson

%Todo 
%     :code to make partial/normal anonymization of dicom files
%     
%     : code to make full anonymization of dicom files
%
%     : code to make corrective anonymization of dicom files by removing
%     the added tags
%    
%    :rewrite find patient details to be compatible with import to work
%    list

macro_helper(varargin{:});
[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%-----------------------------
function initgui(filetype)
%-----------------------------
global SET DATA WL

if ~isempty(SET)
  myfailed('Please save and close your current analysis in order to pseudonymize multiple files.');
  return
end

%open figure
if isequal(filetype,'mat')
	DATA.GUI.Anonymization = mygui('anonymizationmat.fig');
	WL=[];
elseif isequal(filetype,'dcm')
	DATA.GUI.Anonymization = mygui('anonymizationdicom.fig');
	WL=[];
end
fig = DATA.GUI.Anonymization.fig;

set(fig,'DeleteFcn','anonymization(''deletegui_Callback'')');

updatelist;

%--------------------------
function deletegui_Callback
%--------------------------
global DATA WL

try
  DATA.GUI.Anonymization = close(DATA.GUI.Anonymization);
catch %#ok<CTCH>
  delete(gcf);
  if isa(DATA,'maingui')
    DATA.GUI.Anonymization = [];
  end
end

WL=[];

%--------------------------
function updategui_Callback
%--------------------------
updatelist;

%-----------------------------
function addsinglemat_Callback
%-----------------------------
global DATA

%select file
[filename{1},filepath{1}] = myuigetfile('*.mat','Select .mat-file.');
if isequal(filepath{1},0)||isequal(filename,0)
  myfailed('Function aborted. No file selected.',DATA.GUI.Segment);
  return;
end

[~,filename{1}] = fileparts(filename{1});

%add to worklist
addtolistmat(filepath,filename);

%-----------------------------
function addmatfolder_Callback
%-----------------------------
global DATA

%select folder
pathname = DATA.Pref.datapath;
pathname = myuigetdir(pathname,'Select a folder (including subfolders) with .mat files');

if isequal(pathname,0)
	myfailed('Aborted.');
	return;
end

ext = '.mat';

%find files in folder
filelist = createtree(pathname,ext);
numfiles = length(filelist);

for loop = 1:numfiles
  [filepath{loop},filename{loop}] = fileparts(filelist{loop});
end

%add to work list
addtolistmat(filepath,filename);


%-------------------------------------
function importlist_Callback(filetype)
%-------------------------------------
global DATA WL

temp = pwd;
cd(DATA.Pref.datapath);
[filename,pathname] = myuigetfile(...
  '*.*','Select Excel file to load worklist from');
cd(temp);

fullfilename=[pathname filesep filename];

%Load Excel file
[~, ~,textfromfile] = xlsread(fullfilename);
%remove header and extract correct amount of columns
WLfromfile=textfromfile(2:end,1:4);

%format variable anonymize in WL from file 

if isempty(WL)
	WL=WLfromfile;
else
	WL=[WL; WLfromfile]; %add to current WL
end

updatelist;

%-------------------------------------
function clearlist_Callback(filetype)
%-------------------------------------
global WL

WL=[];
updatelist;

%-----------------------------------
function savelist_Callback
%-----------------------------------
global WL DATA

getlistfromgui;

gui=DATA.GUI.Anonymization;
titles=get(gui.handles.worklisttable,'ColumnName');
WLtofile=[titles';WL];%add header row to file

segment('cell2clipboard',WLtofile,true);%true in order to writetofile

%-------------------------------
function adddicomfolder_Callback
%-------------------------------
global DATA

if isempty(DATA)
	pathname = myuigetdir(pwd,'Select directory to pseudonymize DICOM files in. Note: recursive behaviour.');
else
	pathname = myuigetdir(DATA.Pref.datapath,'Select directory to pseudonymize DICOM files in. Note: recursive behaviour.');
end
if isequal(pathname,0)
	myfailed('Aborted.',DATA.GUI.Segment);
	return;
end

filesepindex=strfind(pathname,filesep);
lastfilesep=find(filesepindex<length(pathname),1,'last');
foldername=pathname(filesepindex(lastfilesep)+1:end);

%reformat pathname and foldername to cell and add to list
pathnametolist{1}=pathname;
foldernametolist{1}=foldername;
addtolistdicom(pathnametolist,foldernametolist);
updatelist;

%-----------------------------
function anonymizemat_Callback
%-----------------------------
%Anonymize .mat-files according to anonymization type selcted by
%radiobuttons

global DATA SET WL

if isempty(WL)
  myfailed('Worklist is empty. Please add files to pseudonymize.');
  return;
end

%get type of anonymization from radiobutton
gui=DATA.GUI.Anonymization;
if get(gui.handles.fullradiobutton,'value')
  anontype='full';
elseif get(gui.handles.partialradiobutton,'value')
  anontype='partial';
elseif get(gui.handles.correctiveradiobutton,'value')
  anontype='corrective';
else
  myfailed('Select which type of pseudonymization to perform.');
  return;
end

if isequal(anontype,'full')
  stri = 'This feature will remove all subject identity and info (patient ID, name, sex, age, length, weight, acquisition date and filename) on the selected .mat files. Are you sure?';
elseif isequal(anontype,'partial')
  stri = 'This feature will remove subject identity (patient ID, name and filename) on the selected .mat files. Are you sure?';
else
  stri = 'This feature will only remove original pathname and filename on the selected .mat files. Are you sure?';
end

if ~yesno(stri)
  mymsgbox('Aborted.','');
  return;
end

saveasnewname=get(gui.handles.renamecheckbox,'value');
savekeyfile=get(gui.handles.savelistcheckbox,'value');

%update list from gui
getlistfromgui;

%find which files to anonymize from worklist
anonymize=zeros(size(WL,1),1);
[anonymize(:)]=[WL{:,1}];
indextoanonymize=find(anonymize);
numfiles=length(indextoanonymize);

%doublecheck with user if files shall not be renamed
if not(saveasnewname)
	if isequal(anontype,'full')
		if ~yesno('You have not chosen to save .mat files as new name. For extensive pseudonymization, it is recommended to rename the files. Are you sure that you want to continue without saving .mat-files as new name?')
			return;
		end
	elseif isequal(anontype, 'partial')
		if ~yesno('You have not chosen to save .mat files as new name. The new name will therefore only be changed inside the file. Are you sure that you want to continue without saving .mat-files as new name?')
			return;
		end
	end
end

%doublecheck with the user if key file shall be saved
if not(savekeyfile)
	if ~yesno('You have not chosen to save the key as Excel file. This will be your only key to identify the subjects. Are you sure that you want to continue without saving the key file?');
		return;
	end
end

%initialise output
if savekeyfile
	output{1,1}='New Name';
	output{1,2}='New FileName';
	output{1,3}='File Name';
	output{1,4}='File Path';
	output{1,5}='Patient Name';
	output{1,6}='Patient ID';
	output{1,7}='Birth Date';
	if isequal(anontype,'full')
		output{1,8}='Sex';
		output{1,9}='Age';
		output{1,10}='Acquisition Date';
		output{1,11}='Length';
		output{1,12}='Weight';
		output{1,13}='BSA';
		output{1,14}='Institution';
	end
end

try
	
	%initialize for anonymization
	setstruct = [];
	info=[]; %#ok<NASGU>
	h = mywaitbarstart(numfiles,'Pseudonymizing files. Please wait.',1);
	corruptedfiles='';
	segment('filecloseall_Callback',true);
	for fileloop=1:numfiles
		%--- Load file
		DATA.Silent = true; %Turn on "silent" mode to avoid to much update on screen when loading etc.
		listindex=indextoanonymize(fileloop);
		newname=WL{listindex,2};
		filename=WL{listindex,3};
		pathname=WL{listindex,4};
		ext='.mat';
		
		disp(dprintf('Loading %s.',filename));
		%Make sure a fresh start
		
		%---- try
		%Load
		SET=[];
		load([pathname filesep filename ext],'-mat','setstruct','info');
		
		%Assign
		if not(isempty(setstruct))
			SET = setstruct;
			clear setstruct;
			
			%Call to intialize all variables correcly after loaded data.
			openfile('setupstacksfrommat',1);
			segment('renderstacksfrommat');
			
			if savekeyfile
				output{fileloop+1,1}=newname;
				if saveasnewname
					output{fileloop+1,2}=newname;
				else
					output{fileloop+1,2}=filename;
				end
				output{fileloop+1,3}=filename;
				output{fileloop+1,4}=pathname;
				output{fileloop+1,5}=SET(1).PatientInfo.Name;
				output{fileloop+1,6}=SET(1).PatientInfo.ID;
				output{fileloop+1,7}=SET(1).PatientInfo.BirthDate;
				if isequal(anontype,'full')
					output{fileloop+1,8}=SET(1).PatientInfo.Sex;
					output{fileloop+1,9}=SET(1).PatientInfo.Age;
					output{fileloop+1,10}=SET(1).PatientInfo.AcquisitionDate;
					output{fileloop+1,11}=SET(1).PatientInfo.Length;
					output{fileloop+1,12}=SET(1).PatientInfo.Weight;
					output{fileloop+1,13}=SET(1).PatientInfo.BSA;
					output{fileloop+1,14}=SET(1).PatientInfo.Institution;
				end
			end
			
			%anonymize dataset
			if isequal(anontype,'full')
				tools('anonymoustotal_Callback',true,newname);
			elseif isequal (anontype, 'partial')
				tools('anonymous_Callback',true,newname);
			else
				tools('anonymouscorrective_Callback',true,newname);
			end
			
			%Create thumbnails before storing.
			calcfunctions('calcdatasetpreview');
			
			% Set view settings
			DATA.ViewPanels = 1;
			DATA.ViewPanelsType = {'one'};
			DATA.ViewMatrix = [1 1];
			DATA.ThisFrameOnly = 0;
			DATA.CurrentPanel = 1;
			DATA.CurrentTheme = 'lv';
			DATA.CurrentTool = 'select';
			
			%Save the file.
			disp('Saving...');
			%segment('filesaveall_Callback');
			corrupted=segment('checkcorrupteddataforautomaticsave');
			if corrupted
				corruptedfiles=sprintf('%s, %s',corruptedfiles,filename);
				mywarning(dprintf('Image file %s seems to be corrupted from last save. Please load and pseudonymize manually to ensure that the image is not corrupted before saving.',filename));
			else
				if saveasnewname && not(isequal(filename,newname))
					[sucess,m] = mymovefile([pathname filesep filename ext],[pathname filesep newname ext],'f'); %rename file
					if not(sucess)
						myfailed(dprintf('Failed to rename file %s',filename));
					else
						filemenu('saveallas_helper',pathname,newname);
					end
				else
					filemenu('saveallas_helper',pathname,filename);
				end
			end
			segment('filecloseall_Callback',true);
		end
		setstruct=[];
		info=[]; %#ok<NASGU>
		h = mywaitbarupdate(h);
  end
	mywaitbarclose(h);
	
	if isempty(corruptedfiles)
		if isequal(anontype,'full')
			donestri = 'Patient identity and info removed in the selected .mat files.';
		elseif isequal(anontype,'partial')
			donestri = 'Patient identity removed in the selected .mat files.';
		else
			donestri = 'Filename and pathname removed in the selected .mat files.';
		end
		mymsgbox(donestri,'Done!',DATA.GUI.Segment);
	else
		myfailed(dprintf('The following files have not been pseudonymized: %s',corruptedfiles),DATA.GUI.Segment);
	end
	
	segment('filecloseall_Callback',true);
catch me
	mydispexception(me);
	myfailed('Error in pseudonymization process. The files are not correclty pseudonymized. Please review that all filepaths are correct.');
end

if savekeyfile
	segment('cell2clipboard',output,true);%true in order to writetofile
end

DATA.Silent=0;

%-------------------------------
function anonymizedicom_Callback
%-------------------------------
%Anonymize dicom-files according to anonymization type selcted by
%radiobuttons

global DATA WL

if isempty(WL)
  myfailed('Worklist is empty. Please add files to pseudonymize.');
  return;
end

%get type of anonymization from radiobutton
gui=DATA.GUI.Anonymization;
if get(gui.handles.fullradiobutton,'value')
	anontype='full';
elseif get(gui.handles.partialradiobutton,'value')
	anontype='partial';
elseif get(gui.handles.correctiveradiobutton,'value')
	anontype='corrective';
else
	myfailed('Select which type of pseudonymization to perform.');
	return;
end
if isequal(anontype,'full')||isequal(anontype,'corrective')
	myfailed('Only partial pseudonymization is implemented');
	return;
end

if isequal(anontype,'full')
  stri = 'This feature will remove all subject identity and info on the selected DICOM files. Are you sure?';
elseif isequal(anontype,'partial')
  stri = 'This feature will remove subject identity (patient ID, name and birth date) on the selected DICOM files. Are you sure?';
else
	stri = 'This feature will only remove obscure tags which might contain patient name/ID/birthdate without editing the original patient ID/name/birthdate on the selected DICOM files. Are you sure?';
end

if ~yesno(stri)
  mymsgbox('Aborted.','');
  return;
end

savekeyfile=get(gui.handles.savelistcheckbox,'value');

%update list from gui
getlistfromgui;

%find which files to anonymize from worklist
anonymize=zeros(size(WL,1),1);
[anonymize(:)]=[WL{:,1}];
indextoanonymize=find(anonymize);
numfiles=length(indextoanonymize);

%doublecheck with the user if worklist shall be saved
if not(savekeyfile)
	if ~yesno('You have chosen to not save the work list. This will be your only key to identify the subjects. Are you sure that you want to continue without saving the work list?')
		return;
	end
end
try
	silent=true;
	showwaitbar=true;
	h = mywaitbarstart(numfiles,'Pseudonymizing files. Please wait.',1);
	corruptedfiles='';
	for fileloop=1:numfiles
		%--- Load file
		DATA.Silent = true; %Turn on "silent" mode to avoid to much update on screen when loading etc.
		listindex=indextoanonymize(fileloop);
		newname=WL{listindex,2};
		foldername=WL{listindex,3};
		pathname=WL{listindex,4};
		
		disp(dprintf('Loading %s.',foldername));
		
		%recursive anonymisation for current path to the current new name
		dicomanonymize(pathname,newname,silent,showwaitbar);
		h = mywaitbarupdate(h);
  end
	mywaitbarclose(h);
	
	if isempty(corruptedfiles)
		if isequal(anontype,'full')
			donestri = 'Patient identity and info removed in the selected DICOM files.';
		elseif isequal(anontype,'partial')
			donestri = 'Patient identity removed in the selected DICOM files.';
		else
			donestri = 'Filename and pathname removed in the selected DICOM files.';
		end
		mymsgbox(donestri,'Done!',DATA.GUI.Segment);
	else
		myfailed(dprintf('The following files have not been pseudonymized: %s',corruptedfiles),DATA.GUI.Segment);
	end

	segment('filecloseall_Callback',true);
	
catch me
	mydispexception(me);
	mydisp('Error in pseudonymization process. The files are not correclty pseudonymized.');
end

if savekeyfile
	savelist_Callback;
end

DATA.Silent=0;


%---------------------------------------------
function addtolistmat(filepath,filename,newname)
%--------------------------------------------
global WL

if nargin<3
	newname=filename;
end

%initialize WL
if isempty(WL)
	numlines=0;
	WL=cell(1,4);
else
	numlines=size(WL,1);
	getlistfromgui;
end

%add data in WL
numnewlines=length(filepath);
for loop=1:numnewlines
	WL{numlines+loop,1}=true;%anonymize
	WL{numlines+loop,2}=newname{loop};%default:new name same as filename
	WL{numlines+loop,3}=filename{loop};
	WL{numlines+loop,4}=filepath{loop};%filepath in last column
end

updatelist;

%---------------------------------------------
function subdirbatch_Callback %#ok<DEFNU>
%---------------------------------------------

global DATA

%select folder
pathname = DATA.Pref.datapath;
pathname = myuigetdir(pathname,'Select a folder (including subfolders) with DICOM files');

if isequal(pathname,0)
	myfailed('Aborted.');
	return
end

subdirbatch(pathname)


%---------------------------------------------
function subdirbatch(basefolder)
%--------------------------------------------

global WL

%initialize WL
if isempty(WL)
  numlines=0;
  WL = cell(1,4);
else
  numlines = size(WL,1);
  getlistfromgui;
end

flist = dir(basefolder);

addtoline = numlines;

for floop = 1:length(flist)
  f = flist(floop);
  
  if strcmp(f.name, '.') || strcmp(f.name, '..')
    continue
  end
  
  if ~f.isdir
    continue
  end
  
  % OK, we have a subdir, add it to the list
  addtoline = addtoline + 1;
  
  WL{addtoline,1} = true;
  WL{addtoline,2} = f.name;
  WL{addtoline,3} = f.name;
  WL{addtoline,4} = fullfile(f.folder, f.name);
end

updatelist;


%---------------------------------------------
function addtolistdicom(filepath,foldername, newname)
%--------------------------------------------
global WL

if nargin<3
	newname=cell(size(filepath));
	newname{:}='Hidden';
end

%initialize WL
if isempty(WL)
	numlines=0;
	WL=cell(1,4);
else
	numlines=size(WL,1);
	getlistfromgui;
end

%add data in WL
numnewlines=length(filepath);
for loop=1:numnewlines
	WL{numlines+loop,1}=true;%anonymize
	WL{numlines+loop,2}=newname{loop};%default:new name same as filename
	WL{numlines+loop,3}=foldername{loop};
	WL{numlines+loop,4}=filepath{loop};%filepath in last column
end

updatelist;

%------------------
function updatelist
%------------------
global DATA WL

gui=DATA.GUI.Anonymization;
set(gui.handles.worklisttable,'Data',WL);
DATA.GUI.Anonymization=gui;

%-----------------------
function getlistfromgui
%----------------------
global DATA WL

gui=DATA.GUI.Anonymization;
WL=get(gui.handles.worklisttable,'Data');


%-----------------------------
function initguisubject
%-----------------------------
%initiate the GUI for pseydonymization for current .mat-file
global DATA SET

if isempty(SET)
  myfailed('This feature applies to open files. Please open the file to pseudonymize, or use the batch pseudonymization features.');
  return
end

%open figure
DATA.GUI.Anonymization = mygui('anonymizationsubject.fig');

fig = DATA.GUI.Anonymization.fig;

set(fig,'DeleteFcn','anonymization(''deleteguisubject_Callback'')');


%--------------------------
function deleteguisubject_Callback
%--------------------------
%close the GUI for pseydonymization for current .mat-file
global DATA

try
  DATA.GUI.Anonymization = close(DATA.GUI.Anonymization);
catch %#ok<CTCH>
  delete(gcf);
  if isa(DATA,'maingui')
    DATA.GUI.Anonymization = [];
  end
end

%-----------------------------
function anonymizesubject_Callback
%-----------------------------
%Anonymize current .mat-file according to anonymization type selcted by
%radiobuttons

global DATA

handles = DATA.GUI.Anonymization.handles;

if get(handles.partialradiobutton,'value')
  tools('anonymous_Callback'); %partial
elseif get(handles.extensiveradiobutton,'value')
  tools('anonymoustotal_Callback'); %extensive
end

deleteguisubject_Callback;

