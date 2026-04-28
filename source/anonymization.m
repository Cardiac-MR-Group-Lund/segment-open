 function varargout = anonymization(varargin)
%---------------------------------------------
% ANONYMIZATION GUI for pseudonymization of dicom-files or mat-files as
%extensive, partial or corrective pseudonymization, for a single file/folder or for
%a whole working list.

% Written by Jane Tufvesson, addition by Fanny Månefjord
%#ok<*GVMIS> 

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

[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%-----------------------------
function initgui(filetype)
%-----------------------------
%Initialize GUI
%Input filetype can be mat or dcm

global SET DATA WL 

if ~isempty(SET)
  myfailed('Please save and close your current analysis in order to pseudonymize multiple files.');
  return
end

%open figure
if isequal(filetype,'mat')
	gui = mygui('anonymizationmat.fig');
	WL=[];
  gui.ext = '.mat';
elseif isequal(filetype,'dcm')
	gui = mygui('anonymizationdicom.fig');
	WL=[];
  gui.ext = '.dcm';
end
fig = gui.fig;

set(fig,'DeleteFcn','anonymization(''deletegui_Callback'')');
DATA.GUI.Anonymization = gui;
updatelist;

%--------------------------
function copy = copieschecked
%--------------------------
%Function to see if checkbox for creating copies is checked
%Returns 1 if copies should be creates, otherwise 0
global DATA 

gui = DATA.GUI.Anonymization;
copy = get(gui.handles.copiescheckbox,'value'); %should we create copies (1) or not (0)

%--------------------------
function rename = renamechecked
%--------------------------
%Function to see if rename checkbox is checked or not
%Returns 1 if the files should be renamed, otherwise 0
global DATA 

gui = DATA.GUI.Anonymization;
rename = get(gui.handles.renamecheckbox,'value'); %should we create new names (1) or not (0)

%--------------------------
function rename_Callback
%--------------------------
%update worklist in GUI when rename is clicked

updateoutput;
updatelist;

%--------------------------
function copiescheckbox_Callback
%--------------------------
%Function to enable textbox and browse button when copies checkbox is
% checked.
global DATA 

gui = DATA.GUI.Anonymization;
saveascopies = copieschecked;
bgrcolor = [1 1 1]; % white
frgcolor = [0 0 0]; % black
if (saveascopies)
  enablestatus = 'on'; %Browse button enable
else
  enablestatus = 'off'; 
end

set(gui.handles.outputpathtext, 'BackgroundColor',bgrcolor,'ForegroundColor',frgcolor,'Enable',enablestatus);
set(gui.handles.browsepushbutton,'Enable',enablestatus);
updateoutput;

%--------------------------
function updateoutput
%--------------------------
%Update output dir in WL and update in GUI
global WL DATA

gui = DATA.GUI.Anonymization;
ext = gui.ext;
saveascopies = copieschecked;
renameindex=2;

if (strcmp(ext,'.mat')) %let's check if the files should be renamed in the mat case
  rename = renamechecked;
  if ~rename
    renameindex=3;
  end
end

outputdir = gui.handles.outputpathtext.String;
if (strcmp(ext, '.mat'))
  if (saveascopies  && ~isempty(outputdir))
    outputdir = gui.handles.outputpathtext.String;
    outputdir = strcat(outputdir, '\');
    for loop = 1:size(WL,1)
      outputpath = strcat(outputdir, WL{loop,renameindex}, ext);
      WL{loop,5}=outputpath;
    end
  else
    for loop = 1:size(WL,1)
      filepath = fileparts(WL{loop,4});
      outputpath = strcat(filepath, '\', WL{loop,renameindex}, ext);
      WL{loop,5}=outputpath;
    end
  end

else %dcm
  if (saveascopies  && ~isempty(outputdir))
    outputdir = gui.handles.outputpathtext.String;
    %outputdir = strcat(outputdir, '\');
    for loop = 1:size(WL,1)
      %outputpath = strcat(outputdir, WL{loop,2}, ext);
      WL{loop,5}=outputdir;
    end
  else
    for loop = 1:size(WL,1)
      WL{loop,5}=WL{loop,4};
    end
  end
end
updatelist;

%--------------------------
function browsebutton_Callback
%--------------------------
%User selects an output path

global DATA

%select folder
pathname = DATA.Pref.datapath;
pathname = myuigetdir(pathname,'Select a folder to put output files');

if isequal(pathname,0)
  return;
end
gui = DATA.GUI.Anonymization;
gui.handles.outputpathtext.String = pathname;

updateoutput;

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

%Add \ at the end of the path if it's not there
if (pathname(length(pathname))~='\')
  pathname = strcat(pathname, '\');
end

ext = '.mat';

%find files in folder
filelist = createtree(pathname,ext);
numfiles = length(filelist);
if numfiles == 0
  return
end
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

updateoutput;
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
WLtofile=[titles';WL];  %add header row to file

segment('cell2clipboard',WLtofile,true);  %true in order to writetofile

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

gui=DATA.GUI.Anonymization;

if isempty(WL)
  myfailed('Worklist is empty. Please add files to pseudonymize.');
  return;
end

s = gui.handles.outputpathtext.String;

if copieschecked
  if (isempty(s))
    myfailed('Please choose a directory to put the files in if you wish to save the files as copies.');
    return;
  end
  dup = findduplicates;
  
  if (sum(dup)>0 && strcmp(gui.ext, '.mat')) %There are duplicate names
    
    questn = sprintf('%s\n%s\n%s\n%s',...
      dprintf('You have some files with the same name.'),...
      dprintf('If you proceed the files will be overwritten.'),...
      dprintf('You can cancel and change the duplicate names.'),...
      dprintf('Do you wish to proceed?')...
      );

    answcancel = dprintf('Cancel');
    answproceed = dprintf('Proceed');
    answ = myquestdlg(questn,{answcancel,answproceed},answcancel);
    
    if strcmp(answ,answcancel)
      return;
    end
  end
end


%---get type of anonymization from radiobutton

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

%---find which files to anonymize from worklist
anonymize=zeros(size(WL,1),1);
[anonymize(:)]=[WL{:,1}];
indextoanonymize = find(anonymize == 1);
numfiles = length(indextoanonymize);
if numfiles == 0
  myfailed('No valid data to pseudononymize')
  return
end
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

%---initialise output
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
silentstate = DATA.Silent;

try
	
	%initialize for anonymization
	setstruct = [];
	info=[]; %#ok<NASGU>
	h = mywaitbarstart(numfiles,'Pseudonymizing files. Please wait.');
	corruptedfiles = strings(numfiles,1);
	segment('filecloseall_Callback',true);
  
  %If the copy box is checked, the files need to be copied to the output
  %folder
  if(copieschecked)
    copytooutput(indextoanonymize);
  end

	for fileloop = 1:numfiles
		%--- Load file
		DATA.Silent = true; %Turn on "silent" mode to avoid to much update on screen when loading etc.
		listindex=indextoanonymize(fileloop);
		newname = WL{listindex,2};
    if isnumeric(newname)
      newname = num2str(newname);
    end
		filename = WL{listindex,3};
		pathname = WL{listindex,5};
    newdir = fileparts(pathname);
    if (newdir(length(newdir))~='\')
      newdir = strcat(newdir, '\');
    end

		ext ='.mat';
    if ~ischar(newname) || ~ischar(filename) || ~ischar(pathname)
      continue
    end
		fprintf('Loading %s\n.',filename);
		%Make sure a fresh start
		
		%---- try
		%Load
		SET=[];
    currentfile = [newdir filename ext]; 
    try
      load(currentfile,'-mat','setstruct','info');      
    catch
      corruptedfiles(fileloop) = string(filename);      
      continue
    end
		
		%Assign
		if exist('setstruct','var') && not(isempty(setstruct))
			SET = setstruct;
			clear setstruct;
			
			%---Call to intialize all variables correcly after loaded data.
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
			corrupted = segment('checkcorrupteddataforautomaticsave');
      if corrupted
        wrnstr = dprintf('Image file %s seems to be corrupted. Please check manually if successfully pseudonymized.',filename);
        mywarning(wrnstr);
      end
      try
        if saveasnewname && not(isequal(filename,newname))
          [sucess,~] = mymovefile(currentfile, pathname, 'f');  %mymovefile([pathname filesep filename ext],[pathname filesep newname ext],'f'); %rename file

          if not(sucess)
            myfailed(dprintf('Failed to rename file %s',filename));
          else
            filemenu('saveallas_helper',newdir,newname);
          end
        else
          filemenu('saveallas_helper',newdir,filename);
        end
      catch me
        mydispexception(me)
        corruptedfiles = sprintf('%s, %s',corruptedfiles,filename);
      end
      segment('filecloseall_Callback',true);
    else
      corruptedfiles(fileloop) = string(filename);      
		end
		setstruct = [];
		info = []; %#ok<NASGU>
		h = mywaitbarupdate(h);
  end
	mywaitbarclose(h);
	indcorfiles = (corruptedfiles == '');
	
	if all(indcorfiles) % all strings are empty
		if isequal(anontype,'full')
			donestri = 'Patient identity and info removed in the selected .mat files.';
		elseif isequal(anontype,'partial')
			donestri = 'Patient identity removed in the selected .mat files.';
		else
			donestri = 'Filename and pathname removed in the selected .mat files.';
		end
		mymsgbox(donestri,'Done!',DATA.GUI.Segment);
  else
    DATA.Silent = false;
    errstr = dprintf('The following files have not been pseudonymized:');
		myfailed([errstr,corruptedfiles(~indcorfiles)'],DATA.GUI.Segment);
	end
	
	segment('filecloseall_Callback',true);
catch me
	mydispexception(me);
	myfailed('Error in pseudonymization process. The files are not correctly pseudonymized. Please review that all filepaths are correct.');
end

if savekeyfile
	segment('cell2clipboard',output,true);%true in order to writetofile
end

DATA.Silent = silentstate;

%-------------------------------
function copytooutput(indextoanonymize)
%-------------------------------
%Function to copy files to output directory

global WL DATA

gui = DATA.GUI.Anonymization;
ext = gui.ext;
outputpath = DATA.GUI.Anonymization.handles.outputpathtext.String;

if (strcmp(ext, '.dcm'))
  for loop=1:length(indextoanonymize)
    inputpath = WL{indextoanonymize(loop),4};
    output = strcat(outputpath, '\', WL{indextoanonymize(loop),2});
    copyfile(inputpath, output);
  end
else %mat file

  if ~isfolder(outputpath)
    status = mkdir(outputpath);
    if (status~=1)
      str = dprintf('Could not create directory %s',outputpath);
      myfailed(str);
      return;
    end
  end


  for loop=1:length(indextoanonymize)
    path = WL{indextoanonymize(loop),4};
    copyfile(path, outputpath);
  end
end

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

gui=DATA.GUI.Anonymization;

s = gui.handles.outputpathtext.String;
if copieschecked
  if (isempty(s))
    myfailed('Please choose a directory to put the files in if you wish to save the files as copies.');
    return;
  end
end
%---get type of anonymization from radiobutton
gui = DATA.GUI.Anonymization;
if get(gui.handles.fullradiobutton,'value')
  anontype = 'full';
elseif get(gui.handles.partialradiobutton,'value')
  anontype = 'partial';
elseif get(gui.handles.correctiveradiobutton,'value')
  anontype = 'corrective';
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

savekeyfile = get(gui.handles.savelistcheckbox,'value');

%update list from gui
getlistfromgui;

%find which files to anonymize from worklist
anonymize = zeros(size(WL,1),1);
[anonymize(:)] = [WL{:,1}];
indextoanonymize = find(anonymize);
numfiles = length(indextoanonymize);

%If the copy box is checked, the files need to be copied to the output
%folder
copy=copieschecked;
if(copy)
  copytooutput(indextoanonymize);
end

%doublecheck with the user if worklist shall be saved
if not(savekeyfile)
  if ~yesno('You have chosen to not save the work list. This will be your only key to identify the subjects. Are you sure that you want to continue without saving the work list?')
    return;
  end
end

if savekeyfile
  output{1,1}='New Name';
  output{1,2}='New FolderName';
  output{1,3}='Folder Name';
  output{1,4}='Folder Path';
  output{1,5}='Patient Name';
  output{1,6}='Patient ID';
  output{1,7}='Birth Date';
end

try
  % start writing log file
  logfilename = [getpreferencespath filesep sprintf('pseudonymization_errorlog_%s.log',datestr(now,'yyyymmddHHMMSS'))];
  msg = sprintf('Pseudonymization started at %s',datestr(now,'yyyy-mm-dd HH:MM:SS'));
  writelog(msg,logfilename);

  corruptedfiles = strings(numfiles,1);
  silentstate = DATA.Silent;
  DATA.Silent = true; %Turn on "silent" mode to avoid to much update on screen when loading etc.

  for fileloop = 1:numfiles
    %--- Load file
    writelog('-------------------');
    listindex = indextoanonymize(fileloop);
    newname = WL{listindex,2};
    foldername = WL{listindex,3};
    pathname = WL{listindex,5};
    if (copy)
      pathname = strcat(WL{listindex,5}, '\', newname);
    end

    fprintf('Loading %s.\n',foldername);

    %recursive anonymisation for current path to the current new name
    try
      [name, id, date] = anonymizealldicomfiles(pathname, newname,logfilename,fileloop,numfiles);
    catch me
      mydispexception(me);
      msg = sprintf('Failed for %s', pathname);
      writelog(msg);
      writelog(me.message);
      corruptedfiles(fileloop) = string(pathname);
    end

    if savekeyfile %write the patient info in the excel file
      output{fileloop+1,1}=newname;
      output{fileloop+1,2}=newname;
      output{fileloop+1,3}=foldername;
      output{fileloop+1,4}=pathname;
      output{fileloop+1,5}=name;
      output{fileloop+1,6}=id;
      output{fileloop+1,7}=date;
    end


  end

  indcorfiles = (corruptedfiles == '');
  donestri = '';
  if all(indcorfiles) % all strings are empty
    if isequal(anontype,'full')
      donestri = 'Patient identity and info removed in the selected DICOM files.';
    elseif isequal(anontype,'partial')
      donestri = 'Patient identity removed in the selected DICOM files.';
    else
      donestri = 'Filename and pathname removed in the selected DICOM files.';
    end
  end

  segment('filecloseall_Callback',true);

catch me
  mydispexception(me);
  disp('Error in pseudonymization process. The files are not correctly pseudonymized.');
end

isopen = openerrorlogfile(logfilename);
if ~isopen
  % error log file was not opened
  mymsgbox(donestri,'Done!',DATA.GUI.Segment);
else
  myfailed('Pseudonymization failed on one or several files. Please check the log-file for details',DATA.GUI.Segment)
end

if savekeyfile
  %savelist_Callback;
  segment('cell2clipboard',output,true); %true in order to writetofile
end

DATA.Silent = silentstate;

%------------------------------------
function shouldopen = openerrorlogfile(logfilename)
%------------------------------------
%Open log file for this session in browser.
global DATA
DATA.Silent = false;
if ispc
  if exist(logfilename,'file')
    filetext = fileread(logfilename);
    shouldopen = contains(filetext,["error","file"],'IgnoreCase',true);
    if shouldopen
      dos(sprintf('notepad.exe "%s" &',logfilename));
    end
  else  
    myfailed('Error log file does not exist.')
  end
else
  if exist(logfilename,'file')
    filetext = fileread(logfilename);
    errorexists = contains(filetext,["error","file"],'IgnoreCase',true);
    if errorexists
      msgstr = sprintf('Warning: There was an error. Error log file for this run is %s',logfilename);
    else
      msgstr = sprintf('Error log file for this run is %s',logfilename); 
    end
    mymsgbox(msgstr,'',DATA.GUI.Segment); 
  else
    myfailed('Error log file does not exist.')
  end
end


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
	WL=cell(1,5);
else
	numlines=size(WL,1);
	getlistfromgui;
end

%add data in WL
numnewlines=length(filepath);
for loop=1:numnewlines
  %add a \ at the end of the path
  if (filepath{loop}(length(filepath{loop}))~='\')
    filepath{loop} = strcat(filepath{loop}, '\');
  end
  ipath = strcat(filepath{loop}, newname{loop}, '.mat'); 

	WL{numlines+loop,1}=true;%anonymize
	WL{numlines+loop,2}=newname{loop};%default:new name same as filename
	WL{numlines+loop,3}= filename{loop}; 
	WL{numlines+loop,4}=ipath; %filepath in last column
end
updateoutput;

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
  WL = cell(1,5);
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
updateoutput;
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
	WL=cell(1,5);
else
	numlines=size(WL,1);
	getlistfromgui;
end

%add data in WL
numnewlines=length(filepath);
for loop=1:numnewlines
	WL{numlines+loop,1}=true;%anonymize
	WL{numlines+loop,2}=newname{loop};
	WL{numlines+loop,3}=foldername{loop};
	WL{numlines+loop,4}=filepath{loop};%filepath 
end

updateoutput;
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

%-----------------------------------------------
function dup = findduplicates
%-----------------------------------------------
%Function to see if they are duplicate new file names
global WL

%find which files to anonymize from worklist
anonymize = zeros(size(WL,1),1);
[anonymize(:)] = [WL{:,1}];
indextoanonymize = find(anonymize);
numfiles = length(indextoanonymize);
names = strings(numfiles,1);

for loop = 1 : numfiles
  names(loop) = WL{indextoanonymize(loop),2};
end

% Unique values
[~,idxu,idxc] = unique(names);
% count unique values
[count, ~, idxcount] = histcounts(idxc,numel(idxu));
% Where is greater than 1 occurence
dup = count(idxcount)>1;

%-----------------------------------------------
function worklist_Callback
%-----------------------------------------------
%Update WL when table is updated

getlistfromgui;
updateoutput;

%-----------------------------------------------
function generate_Callback
%-----------------------------------------------

global WL

if isempty(WL)
  myfailed('Worklist is empty. Please add files to pseudonymize.');
  return;
end

%find which files to anonymize from worklist
anonymize = zeros(size(WL,1),1);
[anonymize(:)] = [WL{:,1}];
indextoanonymize = find(anonymize == 1);
numfiles = length(indextoanonymize);
if numfiles == 0
  myfailed('No files chosen.')
  return
end

getlistfromgui;
generatefilenames('init');
updatelist;

