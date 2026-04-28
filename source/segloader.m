classdef segloader < handle  
  %class for loading data into the software
  
  properties(Access = private)
    type
    dicoms
    set
    setpreview; %A preview of a .mat file if existing, otherwise [].
    setinfo; % Info about a .mat file if existing, otherwise [].
    previewmode;
  end
  
  methods
    %-----------------------
    function self = segloader()
    %-----------------------
    % Constuctor. Initiate all properties.
    self.type = 'empty';
    self.dicoms = [];
    self.set = [];
    self.setpreview = [];
    self.setinfo = [];
    self.previewmode = false;
    end
    
    %-----------------------
    function [instancenumber] = getcanonfields(self)
    %-----------------------
    if ~isempty(self.dicoms)
      instancenumber = self.dicoms(1).getinstancenumber;
    else
      instancenumber = [];
    end
    end
    
    %-----------------------
    function setpreviewmode(self)
    %-----------------------
    %Set preview mode to on
    self.previewmode = true;
    end
        
    %-----------------------
    function addmatfile(self, filename)
    %-----------------------
    % Adds a mat file to the object
    
    % Check that the type is empty
    if not(strcmp(self.type, 'empty'))
      error('SEGMENT:ERROR', ...
        'Loader must be empty when loading multidata files');
    end
    
    % Try to load the file
    try
      % append  '\\?\' to avoid MAX_PATH_LENGTH problem
      %first check that it is Windows, long filepath and that the
      %files is local (does not start with \\)
      if ispc && length(filename) > 250 && ~strcmp(filename(1:2), '\\')
        filename = append('\\?\',filename);
      end
      if self.previewmode
        s = load(filename, '-mat', 'preview', 'info');
      end
      if not(self.previewmode) || not(isfield(s,'preview')) || not(isfield(s,'info')) || isempty(s.preview)
         s = load(filename, '-mat', 'setstruct');
      end
    catch e
      error('SEGMENT:ERROR', 'Invalid mat file error was: %s', ...
        e.message);
    end
    
    self.type = 'set';
    
    if self.previewmode && isfield(s,'preview') && isfield(s,'info')
      self.setpreview = s.preview;
      self.setinfo = s.info;
      fnames = fieldnames(self.setinfo);
      for i=1:numel(fnames)
        if isnumeric(self.setinfo.(fnames{i}))
          self.setinfo.(fnames{i}) = sprintf('%g', self.setinfo.(fnames{i}));
        end
      end
    else
      self.set = s.setstruct;
    end
    
    segloaderprogressbar('update', struct('name', 'readmat'));
    end
      
    %-----------------------
    function adddicomfiles(self, filenames)
    %-----------------------
    % Adds DICOM files to the object
    
    % Check that the type isn't 'set'
    if strcmp(self.type, 'set')
      error('SEGMENT:ERROR', ...
        'Can''t continue loading after multidatafile');
    end
    
    % Read the dicoms   
    tic
    dicomdata = segdicomtags.readfiles(filenames);
    t = toc;
    logdisp(sprintf('Loading DICOMs took %0.5g s',t));
    
    if numel(dicomdata) == 0
      error('SEGMENT:ERROR', 'No data to load');
    end

    % Check if we loaded a segdicom file
%     allHasSegmentData = 1;
%     noneHasSegmentData = 1;
  
    hassegmentdata = true(numel(dicomdata),1);
    for i=1:numel(dicomdata)
      hassegmentdata(i) = dicomdata(i).hassegmentdata(); 
%       if dicomdata(1).hassegmentdata()
%         noneHasSegmentData = 0;
%       else
%         allHasSegmentData = 0;
%       end
    end
    allHasSegmentData = all(hassegmentdata);
    noneHasSegmentData = all(not(hassegmentdata));
    if not(allHasSegmentData) && not(noneHasSegmentData)
      error('SEGMENT:ERROR', 'Either all files or none should be segdicom files');      
    end
    
    if allHasSegmentData && strcmp(self.type, 'empty')
      self.type = 'set';
      self.set = [];
      for i=1:numel(dicomdata)
        % get data
        newdata = dicomdata(i).getsegmentdata();
        if ~isrow(newdata)
          newdata = newdata';
        end
        self.set = [self.set , newdata];
      end
      return
    end
    segloaderprogressbar('update', struct('name', 'hassegmentdata1'));
    
    % if we didn't load a segdicom file SegmentData should not be present
    for i=1:length(dicomdata)
      if dicomdata(i).hassegmentdata()
        error('SEGMENT:ERROR', 'Unexpected SegDicom file detected');
      end
    end
    segloaderprogressbar('update', struct('name', 'hassegmentdata2'));
    
    % Check if the dicoms want to be ignored
    i = 1;
    while i <= numel(dicomdata)
      if dicomdata(i).ignoreme()
        dicomdata(i) = [];
      else
        i = i + 1;
      end
    end
    segloaderprogressbar('update', struct('name', 'ignoreme'));
    
    if numel(dicomdata) > 0
      self.type = 'stacks';
      self.dicoms = [self.dicoms dicomdata];
    end
    segloaderprogressbar('update', struct('name', 'adddicomsdone'));
    end
        
    %-----------------------
    function [type, r, ignoreddata] = render(self, datapath, cropbox,donotignore)
    %-----------------------
    % Renders the images in the loader object to a
    % Preview or a set struct suitable for passing on
    % to segment.m.
    ignoreddata=[];
    if nargin < 4
      donotignore = false; %used in case stacks can be loaded anyways, for example for research use in Segment
    end
    switch self.type
      case 'empty'
        error('SEGMENT:ERROR', 'Can''t render empty loader')
      case 'stacks'
        type = 'stacks';
        [r,ignoreddata] = renderstacks(self,datapath,cropbox,donotignore);        
      case 'set'
        type = 'set';
        r = self.set;
        segloaderprogressbar('update', struct('name', 'rendermat'));
      otherwise
        error('SEGMENT:PANIC', 'Unknown segloader type');
    end
    end
    
    %-----------------------
    function [im, cmap, photometricinterpretation, desc, filetype, resolutionx, resolutiony, cancrop] = renderpreview(self) 
    %-----------------------    
    % Renders a preview of the files in the loader object.
    switch self.type
      case 'set'
        % Get im, a preview of the mat file.
        if isempty(self.setpreview) || isempty(self.setinfo)
          im = self.set(1).IM(:, :, 1, 1);
          cmap = [];
          photometricinterpretation = [];
          [acquisitiondate, acquisitiontime] = helperfunctions('getmostrecentstudydate',self.set);
          tfrac = acquisitiontime/(3600*24);
          acquisitiontimestr = sprintf('%02.0f:%02.0f:%02.0f',segloader.hour(tfrac),segloader.minute(tfrac),segloader.second(tfrac));
          [agedigit,ageunit] = calcfunctions('calcagewithunits',segloader.getpreviewdata(self.set(1).PatientInfo, 'Age'));
          desc = sprintf([...
            'Type:%s\n' ...
            'Seq:%s\n' ...
            'Ser:%s\n' ...
            'AcquisitionDate:\t%s\n',...
            'AcquisitionTime:\t%s\n',...
            'Name :\t\t%s\n' ...
            'ID :\t\t%s\n' ...
            'BirthDate :\t%s\n' ...
            'Sex :\t\t%s\n' ...
            'Age :\t\t%s %s\n' ...
            'HeartRate :\t%0.5g\n'], ...
            segloader.getpreviewdata(self.set(1), 'ImageType'), ...
            segloader.getpreviewdata(self.set(1), 'SequenceName'), ...
            segloader.getpreviewdata(self.set(1), 'SeriesDescription'), ...
            acquisitiondate, ...
            acquisitiontimestr, ...
            segloader.getpreviewdata(self.set(1).PatientInfo, 'Name'), ...
            segloader.getpreviewdata(self.set(1).PatientInfo, 'ID'), ...
            segloader.getpreviewdata(self.set(1).PatientInfo, 'BirthDate'), ...
            segloader.getpreviewdata(self.set(1).PatientInfo, 'Sex'), ...
            agedigit,ageunit, ...
            segloader.getpreviewdata(self.set(1).PatientInfo, 'HeartRate') );
        else
          im = self.setpreview;
          cmap = [];
          photometricinterpretation = [];
          [agedigit,ageunit] = calcfunctions('calcagewithunits',self.setinfo.Age);
          desc = sprintf([...
            'Type:%s\n' ...
            'AcquisitionDate:\t%s\n',...
            'Name :\t\t%s\n' ...
            'ID :\t\t%s\n' ...
            'BirthDate :\t%s\n' ...
            'Sex :\t\t%s\n' ...
            'Age :\t\t%s %s\n'], ...
            self.setinfo.ImageType, ...
            self.setinfo.AcquisitionDate, ...
            self.setinfo.Name, ...
            self.setinfo.ID, ...
            self.setinfo.BirthDate, ...
            self.setinfo.Sex, ...
            agedigit, ageunit);
        end
        
        filetype = self.type;
        
        resolutionx = 1;
        resolutiony = 1;
        cancrop = false;
      case 'stacks'
        images = self.dicoms(1).getimages();
        im = images(1).im;        
        cmap = self.dicoms(1).getcolormap();
        if isfield(images,'cmap')
          if ~isempty(images(1).cmap)
            cmap = images(1).cmap;
          end
        end
        if isfield(images,'photometricinterpretation')
          photometricinterpretation = images(1).photometricinterpretation;
        else
          photometricinterpretation = [];
        end
        patinfo = self.dicoms(1).getpatientinfo();
        acquisitiontime = self.dicoms(1).getacquisitiontime();
        if not(ischar(acquisitiontime))
          acquisitiontime = num2str(acquisitiontime);
        end
        secnbr = str2double(acquisitiontime);
        hours = floor(secnbr/3600);
        minutes = floor(mod(secnbr,3600)/60);
        seconds = floor(mod(secnbr,60));
        acquisitiontime = sprintf('%02.0f:%02.0f:%02.0f',hours,minutes,seconds);
        [agedigit,ageunit] = calcfunctions('calcagewithunits',patinfo.Age);
        desc = sprintf([...
          'Type:%s\n' ...
          'Seq:%s\n' ...
          'Ser:%s\n' ...
          'SlicePosition:\t%0.5g\n',...
          'AcquisitionDate:\t%s\n',...
          'AcquisitionTime:\t%s\n',...
          'Name :\t\t%s\n' ...
          'ID :\t\t%s\n' ...
          'BirthDate :\t%s\n' ...
          'Sex :\t\t%s\n' ...
          'Age :\t\t%d %s\n' ...
          'HeartRate :\t%0.5g\n'], ...
          self.dicoms(1).getimagetype(), ...
          self.dicoms(1).getsequencename(), ...
          self.dicoms(1).getseriesdescription(), ...
          self.dicoms(1).getsliceposition(),...
          patinfo.AcquisitionDate, ...
          acquisitiontime,...
          patinfo.Name, ...
          patinfo.ID, ...
          patinfo.BirthDate, ...
          patinfo.Sex, ...
          agedigit,ageunit,...
          ...patinfo.Age, ...
          self.dicoms(1).getheartrate() );
        
        filetype = self.type;
        
        resolutionx = self.dicoms(1).getresolutionx();
        resolutiony = self.dicoms(1).getresolutiony();
        cancrop = true;
      case 'empty'
        im = zeros(256);
        cmap = [];
        photometricinterpretation = [];
        desc = '';
        resolutionx = 1;
        resolutiony = 1;
        cancrop = false;
        filetype = 'empty';
    end
    end
end
  
methods (Access = private)

    %-----------------------
    function [r, ignoreddata] = renderstacks(self, datapath, cropbox, donotignore)
    %-----------------------
    % This rendering method is used when the files loaded are
    % dicom files. It's called from the render method.
    
    % Init rawstacks
    ignoreddata = [];
    if ~self.hasgantrydetectiontilt()
      %Normal case
      
      lines = self.uniquelines();
      [itis, imaxis] = self.isrotated;
      if itis
        r = self.renderrotstacks(datapath,cropbox,lines,imaxis);
        return
      end
      r = [];
      segloaderprogressbar('update', ...
        struct('name', 'renderstacksstart', 'numstacks', numel(lines)));
      for i = 1:numel(lines)
        curstack = segrawstack(lines(i));
        curdicoms = self.dicoms;
        j = 1;
        while j <= numel(curdicoms)
          if curstack.ismatch(curdicoms(j))
            j = j + 1;
          else
            curdicoms(j) = [];
          end
        end
        segloaderprogressbar('update', struct('name', 'sort'));
        [curdicoms, dcminfo] = self.removeduplicates(curdicoms);
        curstack.setdicoms(curdicoms,dcminfo);
        try
          [stack, excludecurrent] = curstack.render(datapath,cropbox,donotignore);
        catch me
          stack = [];
          excludecurrent.series = string(curdicoms(1).getfullseriesdesc);
          excludecurrent.msg = me.message;
          if strcmp(me.identifier, 'SEGMENT:ERROR')
            % display error and series description for stack to be excluded
            fprintf('%s\nIn series %s:\n%s\n',me.identifier,excludecurrent.series,me.message);
          elseif strcmp(me.identifier, 'SEGMENT:FATALERROR')
            rethrow(me)
          else
            mydispexception(me)
          end
        end
     
        if ~isempty(stack)
          r = [r stack]; %#ok<AGROW>
        else
          ignoreddata = [ignoreddata; excludecurrent]; %#ok<AGROW>
        end
        clear curstack;
        segloaderprogressbar('update', struct('name', 'renderstack'));
      end %loop over lines
    else
      %Gantry tilt case
      myworkoff;

      curdicoms = self.dicoms;
      r = [];

      if yesno([dprintf('Gantry tilt is not supported.') ' ' dprintf('Convert?')])
        gantrytilttool(); %Call without any input arguments uses files2load
        mymsgbox(dprintf('A new gantry tilt modified serie was created. Please use the new one instead when loading.'),'',[],true);
        ignoreddata = [];
      else
        failmsg = dprintf('Gantry tilt is not supported.');
        excludecurrent.series = string(curdicoms(1).getfullseriesdesc);
        excludecurrent.msg = failmsg;     
        ignoreddata = excludecurrent;
      end

      segloaderprogressbar('close'); %, struct('name', 'renderstack'));

    end    
    end
    
    %-----------------------
    function [r, imaxis] = isrotated(self)
    %-----------------------
    %Checks if the loaded files is a rotated image stack.
    
    imaxis = [];
    % Check that all types is 'mag'
    for i=1:numel(self.dicoms)
      if not(strcmp(self.dicoms(i).gettype(), 'mag'))
        r = false;
        return
      end
    end
    
    % Get normals and check that it's enough of them
    normals = self.uniquenormals();
    if size(normals, 1) < 6
      r = false;
      return
    end
    
    % Check that the number of images is divisible by the number
    % unique normals. There should be the same number of
    % timeframes in each normal.
    if mod(numel(self.dicoms), size(normals, 1)) ~= 0
      r = false;
      return
    end
    
    % Check that all normal is orthogonal to the same axis.
    [~,~,v] = svd(normals);
    imaxis = v(:,3)';
    if max(abs(acos(imaxis*normals')-pi/2)) > 10*pi/180
      r = false;
      return
    end
    
    r = true;
    
    end
    
    %-----------------------
    function r = renderrotstacks(self, datapath, cropbox, lines, imaxis)
    %-----------------------
    % This rendering method is used when the files loaded are
    % dicom files. It's called from the render method.
    
    % Init rawstacks
    segloaderprogressbar('update', ...
      struct('name', 'renderstacksstart', 'numstacks', numel(lines)));
    curstack = rotrawstack(lines, imaxis);
    curdicoms = self.dicoms;
    j = 1;
    while j <= numel(curdicoms)
      if curstack.ismatch(curdicoms(j))
        j = j + 1;
      else
        curdicoms(j) = [];
      end
    end
    segloaderprogressbar('update', struct('name', 'sort'));
    [curdicoms, dcminfo] = segloader.removeduplicates(curdicoms);
    curstack.setdicoms(curdicoms,dcminfo);
    r = curstack.render( datapath, cropbox );
    clear curstack;
    segloaderprogressbar('update', struct('name', 'renderstack'));
    end
        
    %-----------------------
    function r = uniquenormals(self)
    %-----------------------
    % Gets all the unique normal from the dicom files in
    % the loader object. Used in isrotaded method.
    
    % Put all the normals in r
    r = [];
    for i=1:numel(self.dicoms)
      r(end+1, :) = self.dicoms(i).getnormal();
    end
    
    % Sort it
    r = sortrows(r);
    
    % Remove dupes
    n = 1;
    while n < size(r, 1)
      %if isequal(r(n, :), r(n+1, :)) %Old code. Fix for #1443 on row below
      if norm(r(n,:)-r(n+1,:),2)<1e-6
        r(n+1, :) = [];
      else
        n = n+1;
      end
    end
    end
    
    %-----------------------
    function r = hasgantrydetectiontilt(self)
    %-----------------------
    %Returns true if has gantrydetectiontilt
    
    r = ~isequal(self.dicoms(1).getgantrydetectortilt(),0); 

    end
    
    %-----------------------
    function r = uniquelines(self)
    %-----------------------
    % Get all the unique lines from the dicom files in the loader object.
    % A line consist of the orientation vectors (vector parallel to
    % the x and y axis of the picture) and a projection of the
    % imageposition onto the subspaces that the orientation vectors
    % span. We need this projection to separate projection to
    % separate stacks with the same orientation but different positions.
    r = struct('orientation',{}, 'linepos',{}, 'seq',{},'type',{});
        
    segloaderprogressbar('update', struct('name', 'uniquelinesstart', ...
      'numdicoms', numel(self.dicoms)));
    for i=1:numel(self.dicoms)
      curline = self.dicoms(i).getline();
      found = false;
      for j=1:numel(r)
        found = found || segloader.linecmp(r(j), curline);
      end
      if not(found)
        r(end+1) = curline;
      end
      segloaderprogressbar('update', struct('name', 'uniqueline'));
    end
    end
  end
  
  methods(Static = true, Access = private)
    
    %-----------------------
    function [curdicoms,dcmimages] = removeduplicates(curdicoms)
    %-----------------------
    % Removes any duplicates in curdicoms
    numdicoms = numel(curdicoms);
    spacetimepos = zeros(5, numdicoms);
    dcmimages = cell(numdicoms,1);
    segloaderprogressbar('update', struct('name', 'removedupsstart', ...
      'numdicoms', numel(curdicoms)));
    for i = 1:numel(curdicoms)
      s = curdicoms(i).getimages();
      r = s(1).spacetimepos;
      r(4) = s(1).triggertime;
      spacetimepos(:, i) = r';
      dcmimages{i} = s;
      segloaderprogressbar('update', struct('name', 'removedup'));
    end
    tt = sortrows([spacetimepos; 1:size(spacetimepos, 2)]');
    
    i = 1;
    toremove = [];
    while i < size(tt, 1)
      if not(isequal(tt(i, 1:5), tt(i+1, 1:5)))
        i = i+1;
        continue;
      end
      
      if not(curdicoms(tt(i, 6)).isduplicate(curdicoms(tt(i+1, 6))))
        i = i+1;
        continue;
      end
      
      toremove = [toremove tt(i+1, 6)];
      tt(i+1, :) = [];
    end
    
    curdicoms(toremove) = [];
    dcmimages(toremove) = [];
    dcmimages = cell2mat(dcmimages);
    end

    %-----------------------
    function r = getpreviewdata(s, fname)
    %-----------------------
    % Gets an element from a struct if it exist and has type char.
    % else return ''. Used for generating preview.
    if isfield(s, fname) && isa(s.(fname), 'char')
      r = s.(fname);
    else
      r = '';
    end
    end
  end
    
  methods(Static = true, Access = public)
    
    %-----------------------
    function eq = linecmp(line1, line2)
    %-----------------------
    % Used to compare two lines generated by uniquelines.
    % (Accepting some errors).
    eq = false;
    
    or_dif1 = line1.orientation(1,:)*line2.orientation(1,:)';
    or_dif2 = line1.orientation(2,:)*line2.orientation(2,:)';
    tol = cos(2*pi/180); %5 degrees tolerance
    if or_dif1 < tol || or_dif2 < tol
      return
    end
    
    if norm(line1.linepos - line2.linepos) > 2 %changed from 10 NL
      return
    end
    
    phasestrings = {'\P\','\P','PCA/P','VELOCITY\NONE','\PCPHASE'};
    magstrings = {'\M\','PCA/M','M\FFE','T1\NONE','\M','\PCMAG'};
    psirstrings = {'\PSIR'};
    
    % exclude matchstring from type comparison
    matchstr = '\PRIMARY';
    type1 = erase(line1.type,matchstr);
    type2 = erase(line2.type,matchstr);
    
    if ~strcmp(line1.seq, line2.seq) %not equal sequence name
        if ~isapproxequal(norm(line1.linepos - line2.linepos),0,0.1) %20220803 EH&KG: changed tol from 0 to 0.1 to extend tolarance for position difference
            %different sequence name and there is position difference
            return
        end
        if strcmp(type1, type2)
          % exclusion for 4D flow phase stacks, since they are
          % same type but different sequence name
          %TODO check if \P last letter (UIH images)
          if cellfun(@(x,y) (contains(x,phasestrings)) && (contains(y,phasestrings)),{type1},{type2})
            
           if cellfun(@(x,y) (contains(x,psirstrings)) || (contains(y,psirstrings)),{type1},{type2})
             return
           end
           
          else
           return
          end
        end
    end

    if ~strcmp(type1, type2) %not equal image type
        % exclusion for 2D/4D flow phase and magnitude stacks
        %TODO check if \P last letter (UIH images)
         if cellfun(@(x,y) (contains(x,magstrings)) && (contains(y,phasestrings)) || (contains(y,magstrings)) && (contains(x,phasestrings)),{type1},{type2})  
           if cellfun(@(x,y) (contains(x,psirstrings)) || (contains(y,psirstrings)),{type1},{type2})  
             return
           end           
         else
           return
%         elseif cellfun(@(x,y)(contains(x,{'\M\IR','\CR\IR','\M\FFE',}) && contains(y,{'\M\IR','\CR\IR','\M\FFE'})),{line1.type},{line2.type}) %PSIR PHILIPS
%           return
        end
    end
    
    eq = true;
    end
    
    %-----------------------
    function h = hour(tfrac)
    %-----------------------
    h = floor(tfrac*24);
    end
    % Convert from tfrac to hours
    
    %-----------------------
    function m = minute(tfrac)
    %-----------------------
    % Convert from tfrac to minutes
    tfrac = tfrac-segloader.hour(tfrac)/24;
    m = floor(tfrac*24*60);
    end
    
    %-----------------------
    function s = second(tfrac)
    %-----------------------
    % Convert from tfrac to seconds
    tfrac = tfrac-segloader.hour(tfrac)/24;
    tfrac = tfrac-segloader.minute(tfrac)/(24*60);
    s = floor(tfrac*24*60*60);
    end
    
  end
end
