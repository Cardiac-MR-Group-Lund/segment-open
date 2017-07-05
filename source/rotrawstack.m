classdef rotrawstack < segrawstack
  
  properties (Access = private)
    lines %a rotrawstack has several lines
    angledim
    angles
  end
  
  methods
    %-----------------------
    function self = rotrawstack(lines, imaxis)
    %-----------------------
    %Constuctor. Initiate all properties.
    numlines = numel(lines);
    angs = zeros(1,numlines);
    breakall = false;
    %Loop over orientations to see which one is rotation axis
    for ori = 1:2
      lix = 1;
      loop = true;
      %Loop until one lines at the end of rotation is found
      while loop
        for i = 1:numlines
          angs(i) = acos(lines(lix).orientation(ori,:)*lines(i).orientation(ori,:)');
        end
        angdiff = diff(sort(angs));
        if all(angdiff >= 0.01)
          %Found proper rotation
          loop = false;
        else
          [~,lix] = max(angs);
        end
        if isempty(angdiff) || (min(angdiff) < 0.01 && std(angdiff) < 0.01)
          %No rotation found, proceed to next orientation
          break
        end
      end
      if ~loop
        break
      end
    end
    
    self = self@segrawstack(lines(lix));
    self.lines = lines;
    self.normal = [imaxis 0 0];
    self.angledim = ori;
    self.angles = angs;
    end
    
    %-----------------------
    function r = ismatch(self, dicom)
    %-----------------------
    % Check if a dicom files has one of the lines in this stack.
    r = false;
    for i = 1:numel(self.lines)
      if segloader.linecmp(self.lines(i), dicom.getline())
        r = true;
        return
      end
    end
    end
    
    %-----------------------
    function r = render(self, datapath, cropbox)
    %-----------------------
    % Renders the dicom files and images in this stack
    % to a preview and set struct suitable for
    % passing on to segment.m. If cropbox == []
    % load entire image.

    r = render@segrawstack(self, datapath, cropbox);
    cropbox = self.fixcropbox(cropbox);
    r.preview.Rotated = true;
    
    res = [r.preview.ResolutionY r.preview.ResolutionX];
    ori = self.angledim;
    res = res(ori);
    vec = (cropbox(3-ori,1):cropbox(3-ori,2))-1;
    x = zeros(length(self.lines),length(vec));
    y = x;
    z = x;
    for i = 1:length(self.lines)
      lin = self.lines(i);
      stpos = self.images(lin.images(1)).spacetimepos;
      x(i,:) = stpos(1)+vec*res*lin.orientation(ori,1);
      y(i,:) = stpos(2)+vec*res*lin.orientation(ori,2);
      z(i,:) = stpos(3)+vec*res*lin.orientation(ori,3);
    end
    d = diff(x).^2+diff(y).^2+diff(z).^2;
    [~,ix] = min(sum(d));
    
    %Set RotationCenter field (negative if image needs to be rotated)
    r.preview.RotationCenter = (3-2*self.angledim)*ix; %- cropbox(3)-1; %should be 86
    end
    
  end
    
  methods(Access = protected)
    
    %-----------------------
    function [base, offset] = makeimbase(self, dimsizes)
    %-----------------------
    % Return a base for the loaded image stack.
    
    % Change coordinate system of spacetimepos to timelinephase subspace
    spacetimepos = cat(2, self.images(:).spacetimepos);
    timelinephasepos = [0 0 0 1 0; self.normal; 0 0 0 0 1]*spacetimepos;
    angs = self.angles;
    
    % find the offset
    [~, first_time_ind] = min(timelinephasepos(1, :));
    [~, first_line_ind] = min(angs);
    offset = [self.images(first_line_ind).spacetimepos(1:3); ...
      self.images(first_time_ind).spacetimepos(4); 1];
    
    % find timeDist and sliceDist
    bounds = [...
      min(timelinephasepos(1, :)) max(timelinephasepos(1, :)); ...
      min(angs) max(angs)];
    
    if dimsizes(3) > 1
      timeDist = diff(bounds(1, :))/(dimsizes(3)-1);
    else
      timeDist = 0;
    end
    if dimsizes(4) > 1
      sliceDist = diff(bounds(2, :))/(dimsizes(4)-1);
    else
      sliceDist = 1;
    end
    
    % calc the base
    %base = [ [0 0 0 1 0]'*timeDist [implanepos(1:3,1)*sliceDist; implanepos(4:5,1)] [0 0 0 0 1]'];
    base = [ [0 0 0 1 0]'*timeDist self.normal'*sliceDist [0 0 0 0 1]'];
    end
    
    %-----------------------
    function [im,timepos] = makeimdata(self, imbase, offset, dimsizes, cropbox)
    %-----------------------
    % Combines all images in the stack to an im suitable for SET.im.
    self.memchunk = [];
    im = repmat(single(NaN), dimsizes);
    ori = self.angledim;
    vec1 = self.line.orientation(ori,:);
    angdist = norm(imbase(1:3,2))/norm(self.normal);
    imbase(1:3,2) = imbase(1:3,2)/angdist;
    
    %setup to enable finding the correct line
    nbrlines = numel(self.lines);
    posi = zeros(3,nbrlines);
    norma = zeros(3,nbrlines);
    for i = 1:numel(self.lines)
      lin = self.lines(i);
      posi(:,i) = lin.orientation(1,:)*lin.linepos(1) + lin.orientation(2,:)*lin.linepos(2);
      norma(:,i) = cross(lin.orientation(1,:),lin.orientation(2,:));
    end
    [self.lines.images] = deal([]);
    
    segloaderprogressbar('update', struct('name', 'makeimdatastart', ...
      'numimages', numel(self.images)));
    timepos = zeros(1,numel(self.images));
    for i=1:numel(self.images)
      
      % calc IMPos and make sure we get integers
      impos = pinv(imbase)*(self.images(i).spacetimepos - offset);
      stpos = self.images(i).spacetimepos(1:3);
      stproj = repmat(stpos,1,nbrlines) - repmat(stpos'*norma,3,1).*norma;
      [~,ix] = min(sum((stproj-posi).^2));
      self.lines(ix).images = [self.lines(ix).images i];
      impos(2) = acos(self.lines(ix).orientation(ori,:)*vec1')/angdist;
      
      if norm(round(impos)-impos, inf) > 0.2
        error('SEGMENT:ERROR', 'Image is out of place. This means that they are not forming a rectilinear stack.');
      end
      impos = round(impos);
      
      % Check that we haven't already used position and
      % store the image
      if not(isnan(im(1, 1, 1+impos(1), 1+impos(2), 1+impos(3))))
        error('SEGMENT:ERROR', 'Two images in the same spacetime pos');
      end
      
      %Store
      timepos(i) = 1+impos(1);
      im(:, :, timepos(i), 1+impos(2), 1+impos(3)) = ...
        self.images(i).im( ...
        cropbox(1, 1):cropbox(1, 2), ...
        cropbox(2, 1):cropbox(2, 2));
      
      %Update progress bar
      segloaderprogressbar('update', struct('name', 'makeimdata'));
    end
    
    % Check that there is no NaN:s left
    if not(isequal(im, im))
      error('SEGMENT:ERROR', 'Some image data is missing');
    end
    end
       
    %-----------------------
    function dimsizes = getdimensionsizes(self, cropbox)
    %-----------------------
    % Gets the size of each dimension, i.e number of frames,
    % depth, x-size, y-size.
    
    % Change coordinate system of spacetimepos to linephasetime subspace
    spacetimepos = cat(2, self.images(:).spacetimepos);
    linephasetimepos = [self.normal; 0 0 0 0 1; 0 0 0 1 0]*spacetimepos;
    angs = self.angles;
    
    % get nSlices and nPhases
    nSlices = sum(diff(sort(angs)) > 1e-4)+1;
    nPhases = max(linephasetimepos(2, :));
    
    % get nFrame using nSlices and nPhases
    if mod(numel(self.images), nSlices*nPhases) ~= 0
      error('SEGMENT:ERROR', ...
        dprintf('nSlices*nPhases (%d*%d) don''t divide numel(images)=%d',nSlices,nPhases,numel(self.images))); %#ok<SPERR>
    end
    nFrames = numel(self.images)/nSlices/nPhases;
    
    % make dimsizes
    dimsizes = [ ...
      (diff(cropbox(1, :)) + 1) ...
      (diff(cropbox(2, :)) + 1) ...
      nFrames ...
      nSlices ...
      nPhases ];
    end
    
    %-----------------------
    function r = hascollisions(self, imbase, offset)
    %-----------------------
    % Returns true if two images in the object gets
    % the same coordinates using imbase and offset.
    global DATA
    
    silent = DATA.Silent;
    
    spacetimepos = cat(2, self.images(:).spacetimepos);
    
    % calc IMPos and make sure we get integers
    offsetall = offset*ones(1, numel(self.images));
    impos = pinv(imbase)*(spacetimepos - offsetall);
    if norm(round(impos)-impos, inf) > 0.1
      if silent
        %Automated loader routines
        error('SEGMENT:ERROR', hascollisionshelper(self,impos))
      else
        %User loader a file
        stri = hascollisionshelper(self,impos);
        if ~yesno(['Error: ' stri ...
            dprintf('.\n\nReported slice thickness is %0.5g. ',self.getslicethickness) ...
            dprintf('Overriding will impact subsequent quantification by setting thickness to %0.5g. ',self.getslicethickness) ...
            'Do you want to override? (not recommended unless you know exactly what you do)']);
          error('SEGMENT:ERROR', hascollisionshelper(self,impos));
        end;
      end;
    end
    impos = round(impos);
    
    % Check for duplicates
    r = not( isequal(sortrows(impos'), unique(impos', 'rows')) );
    end
    
  end
end