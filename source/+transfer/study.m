classdef study < handle
  properties(SetAccess = private, GetAccess = private)
    temppath
    studypath
    waitbar
    pattern
    newname
  end
  
  methods(Access = public)
    
    function self = study(path, pattern, waitbar)
      
      global DATA
      
      gui = DATA.GUI.Segment;
      
      self.waitbar = waitbar;
      self.pattern = pattern;
      
      %--- Ask for new name for the study
      mymsgbox(dprintf([...
        'You need to give the study a new. Please follow the pattern according to the CRF. '...
        'In this trial the pattern is %s where the characters means: \n\n' ...
        'n any numeric digit (0,1,2,3,4,5,6,7,8,9)\n' ...
        'x any character, (A,B,C....,a,b,c...)\n' ...
        'X (upper case letter). Must be exact this character\n' ...        
        '0 numerical digit. Must be exactly this digit\n' ...
        '- dash\n' ...
        ': colon\n' ...
        '. point\n'],pattern));
      
      docheck = true;
      if isempty(pattern)
        docheck = false;
        pattern = 'Pnnnn';
      end;
      
      valid = false;
      
      %Keep asking until a valid name has been found.
      while ~valid
        s = inputdlg({['Enter new patient name according to the pattern ' pattern]},'Patient Name',1,{pattern});
        if isempty(s)
          error('SEGMENT:ERROR','User aborted on naming study.');
        else
          self.newname = s{1};
        end;
      
        %Ensure no uncessesary spaces
        self.newname = strtrim(self.newname);
      
        %--- Check validity of pattern
        if docheck
          [valid,msg] = transfer.checkpattern(pattern,self.newname); %Throws an error if not matching
          if ~valid
            mymsgbox(dprintf('Invalid name. Pattern is %s. %s',pattern,msg));
          end;
        else
          valid = true;
        end;
      end;
      
      %Enable the progress bar
      set(gui.handles.overallpatch,'visible','on');
      set(gui.handles.overallline,'visible','on');
      set(gui.handles.overalltext,'visible','on');
      
      self.temppath = tempname;
      mkdir(self.temppath);
     
      %--- Sort the files (1).
      transfer.progressbar(0.1,'Sort and copy files.');
      logargs.recurse = true;
      logargs.anonym = false;
      logargs.stable = false; % but slow
      logargs.keepname = false;
      logargs.makethumbs = false;
      logargs.usedesc = true; % sort by series description
      dicomsorter('sortit', path, 1, self.temppath, [], logargs);
      delete([self.temppath filesep 'log.txt']);
      
      % Check that we only have one study and get the patientpath name
      tt = dir(self.temppath);
      numfolders = sum(cat(1,tt(:).isdir));
      
      if isequal(numfolders,2)
        error('SEGMENT:ERROR', 'Didn''t find any valid DICOM files in choosen folder');
      end
      
      %Make a list of studies
      studies = [];
      studies.patient = [];      
      studies.studydate = [];      
      studies.folder = [];
      studies.series = [];
      
      numstudies = 0;
      for loop = 1:length(tt)
        if tt(loop).isdir && (~isequal(tt(loop).name,'.')) && (~isequal(tt(loop).name,'..'))
          %Folder => check for studies
          
          %get files & folders
          f = dir([self.temppath filesep tt(loop).name]);
          
          for sloop = 1:length(f)
            if f(sloop).isdir && (~isequal(f(sloop).name,'.')) && (~isequal(f(sloop).name,'..'))
              %Add to list of studies
              if isempty(strfind(f(sloop).name,'20XX-XX-XX'))
                %Ignored study without valid dates
                numstudies = numstudies+1;
                studies(numstudies).patient = tt(loop).name; %#ok<AGROW>
                studies(numstudies).studydate = f(sloop).name;                                 %#ok<AGROW>
                studies(numstudies).folder = [self.temppath filesep tt(loop).name filesep f(sloop).name]; %#ok<AGROW>

                %Find number of series
                fs = dir(studies(numstudies).folder);
                studies(numstudies).series = sum(cat(1,fs(:).isdir))-2; %#ok<AGROW>
                
                if studies(numstudies).series<0
                  studies(numstudies).series = 0; %#ok<AGROW>
                end;
                
              end;
            end;
            
          end; %Loop over studies
        end; %Valid patient
      end; %Loop over patients

      if numstudies == 0
        error('ERROR:SEGMENT', 'Found no studies.');
      end;
      
      if numstudies>1
        mywarning('Detected multiple studies. This should normally not occur.');
        
        %Create a menu structure
        studycell = cell(1,length(studies));
        for loop = 1:length(studycell)
          studycell{loop} = sprintf('Patient:%s Date:%s Series:%d',studies(loop).patient,studies(loop).studydate,studies(loop).series);
        end;

        selstudy = mymenu('Select which study to take',studycell{:});
        if isequal(selstudy,0)
          error('SEGMENT:ERROR', 'Aborted when chosing which study to take.');
        end;              
        
        if ~yesno(dprintf('Selected study ''%s''. Are you sure?',studycell{selstudy}))
          error('SEGMENT:ERROR', 'Aborted when confirming which patient to upload.');
        end;                          

      else
        selstudy = 1; %There was only one...
      end;            

      self.studypath = studies(selstudy).folder;
        
      %--- Anonymize (2)
      transfer.progressbar(0.2,'Anonymizing data.');      
      
      anonsilent = true;
      anonwaitbar = true;
      dicomanonymize(self.studypath, self.newname, anonsilent, anonwaitbar); %Need only to sort the study path
      disp(sprintf('Sorted folder was:%s', self.studypath));
      
      %Rename the folder
      oldfolder = [self.temppath filesep studies(selstudy).patient];
      newfolder = [self.temppath filesep self.newname ];

      disp(sprintf('Renaming the folder to %s',newfolder));
      
      [ok,msg] = movefile(oldfolder, newfolder);
      if ~ok && ~strcmp(oldfolder,newfolder)
        disp(msg);
        error('SEGMENT:ERROR','Could not rename folder.');
      end;
      
      %Assign the newly renamed folder.
      self.studypath = [self.temppath filesep self.newname filesep studies(selstudy).studydate];
    end
    
    function p = getpath(self)
      p = self.studypath;
    end
    
    function n = getname(self)
      n = self.newname;
    end
    
    function delete(self)
      w = self.waitbar(1, 'Clearing up temporary files...'); %#ok<NASGU>
      transfer.deltree(self.temppath);
      clear w;
    end
  end
end