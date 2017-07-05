classdef segdicomtags < handle
  
  properties
    tags
    newtags
  end
  
  methods
    %---------
    function self = segdicomtags(tags)
    %---------
      % Constructor. Sets the tags property.
      
      self.tags = tags;
      self.newtags = [];
    end
    
    %---------
    function createnewtags(self)
    %---------
    %Create a new struct of tags
      fnames = fieldnames(self.tags);
      for i=1:numel(fnames)
        if(not(isequal(fnames{i}, 'PixelData')))
          try
            self.newtags.(fnames{i}) = self.tags.(fnames{i});
          catch e
          end
        end
      end
    end
    
    %---------
    function switchtags(self)
    %---------
      % Switch to new tags
      clear self.tags;
      self.tags = self.newtags;
      clear self.newtags;
    end
    
    %---------
    function r = isduplicate(self, other)
    %---------
      % Checks if other is the same as self
      r = true;
      s = self.getimages();
      o = other.getimages();
      
      if not(isequal(s.spacetimepos, o.spacetimepos))
        r = false;
        return
      end
      
      if not(isequal(s.triggertime, o.triggertime))
        r = false;
        return
      end
      
      if not(isequal(s.inversiontime, o.inversiontime)) %SBT151206
        r = false;%SBT151206
        return%SBT151206
      end%SBT151206
      
      if single(max(abs(s.im(:) - o.im(:)))) / single(max(s.im(:))) > 1/255 %Was 5/255. does not make sense changed...
        r = false;
        return
      end
    end
    
    %---------
    function r = spacetimepos(self)
    %---------
      % returns the spacetimepos of the images
      s = self.getimages();
      r = s(1).spacetimepos;
      r(4) = s(1).triggertime;
    end
    
    %---------
    function r = ignoreme(self)
    %---------
      % Returns true if this image should be ignored from the loading
      % process
      r = false;
      if not(self.haspixeldata())
        r = true;
        return;
      end
      if isfield(self.tags, 'ImageType')
        if(strfind(char(self.tags.ImageType), '\M\PCA')) %was 'ORIGINAL\PRIMARY\M_PCA\M\PCA'
          r = true;
          return;
        end
        if(strcmp(char(self.tags.ImageType), 'DERIVED\PRIMARY\PROJECTION IMAGE\COLLAPSE '))
          r = true;
          return;
        end
      end
    end
    
    %---------
    function r = haspixeldata(self)
    %---------
      % Checks to see that there is pixeldata
      
      r = isfield(self.tags, 'PixelData');
    end

    %---------
    function normal = getnormal(self)
    %---------
      % Returns the normal to the picture, that i
      % the cross product of the orientation vectors
      
      ori = self.getorientation();
      normal = cross(ori(1, :), ori(2, :));
    end
    
    %---------
    function orientation = getorientation(self)
    %---------
      % Get the orinentation vectors as a 2x3 matrix
      
      % Set default value and return if the not(isfield)
      orientation = [1 0 0; 0 1 0];
      if not(isfield(self.tags, 'ImageOrientationPatient'))
        return
      end
      
      % Try to parse the orientation
      tt = sscanf(char(self.tags.ImageOrientationPatient), ...
        '%f\\%f\\%f\\%f\\%f\\%f');
      if ischar(tt) || length(tt) ~= 6
        % sscanf failed
        return
      end
      orientation = reshape(tt, 3, 2)';
      
      % Normalize orientations and make sure they are othogonal
      orientation(1, :) = orientation(1, :)./norm(orientation(1, :));
      orientation(2, :) = orientation(2, :)./norm(orientation(2, :));
      if abs(orientation(1, :)*orientation(2, :)') > 0.1
        error('SEGMENT:ERROR', ...
          'The ImageOrientation vectors are not orthogonal');
      end
    end
    
    %---------
    function pos = getposition(self)
    %---------  
      
      % Returns the image position in 1x3 vector
      
      % Set default value and check if field exist
      pos = [0 0 0];
      if not(isfield(self.tags, 'ImagePositionPatient'))
        return
      end
      
      % Parse the position
      pos = sscanf(char(self.tags.ImagePositionPatient), '%f\\%f\\%f');
      pos = reshape(pos, 1, 3);     
      
      %Check if contains NaN's or Inf. We have seen this
      %from some VERY odd Bruker DICOM files.
      
      if any(isnan(pos) | isinf(pos))
        pos = [0 0 0];
      end;
      
    end
    
    %---------
    function z = getsliceposition(self)
    %---------
      %Returns slice position
      
      z = sum(self.getnormal.*self.getposition); %Use this instead to look for the tag. It is buggy under certain conditions, this is more stable.
    end;
    
    %---------
    function type = gettype(self)
    %---------
      % Find the image type of a dicom image
      % Old values for type is
      % 0: 'mag', 1: 'phase', 2: Unknown, 3: Unknown     
      
      % Get imagetype and seqname
      if isfield(self.tags, 'ImageType')
        imagetype = char(self.tags.ImageType);
      else
        imagetype = '';
      end
      if isfield(self.tags, 'SequenceName')
        seqname = char(lower(self.tags.SequenceName));
      else
        seqname = '';
      end
      
      % Find the type
      type = 'mag';
      switch self.getscanner()
        case 'Philips'
          if strfind(imagetype, 'P\PCA');
            type = 'phase';
          end;
          if strfind(imagetype, 'M\IR');
            type = 'psir' ;
          end
          if strfind(imagetype, 'CR\IR');
            type = 'phase' ;
          end
          if strfind(imagetype, 'M\FFE');
            type = 'mag';
          end
          
        case 'Siemens'
          % Stacks with seqname that starts with 'p' are phase images
          % except for Siemens 4D
          % flow data, where the magnitude stack starts with 'p' anyway. Therefore,
          % we have to check for '\M\' in imagetype - this should
          % indicate that it's a magnitude stack.
          if strfind(imagetype, '\M\')
            type = 'mag';
          elseif strfind(imagetype, '\P\')
            type = 'phase';
          elseif numel(seqname) > 1 && seqname(1) == 'p'
            type = 'phase';
          elseif strfind(imagetype, '\MAG\')    
            type = 'angio';
          end;
        case 'GE'
          % It's not pretty but it will do for now
          if isfield(self.tags, 'VelocityEncoding') && ...
              any(self.getpixeldata() < 0)
            type = 'phase'; %For now. EH:
            if isfield(self.tags, 'VelocityEncoding') %Added EH:
              if isequal(self.tags.VelocityEncoding,[0 0]) %EH:
                type = 'mag'; %EH:
              end; %EH:
            end; %EH:
          end
        case 'Bruker'
          %Added EH: No documentation exists if this is the appropriate
          %test or not. Just my guess.
          if  isfield(self.tags,'ProtocolName')
            if any(strfind(self.tags.ProtocolName,'Velocity'))
              if any(self.getpixeldata()<0)
                type = 'mag';
              end;
            end;
          end;
      end
    end
    
    %---------
    function scanner = getscanner(self)
    %---------
      % Return the scanner based on the Manufacturer tag
      
      % Set the default scanner
      scanner = 'GE';
      
      % Get the manufacturer string
      if not(isfield(self.tags, 'Manufacturer'))
        return
      end
      stri = lower(char(self.tags.Manufacturer));
      if isempty(stri)
        stri = lower(char(self.tags.ManufacturersModelName));
      end
      
      silent = true;
      scanner = extractscanner(stri,silent);
          
    end
    
    %---------
    function line = getline(self)
    %---------
      % Gets the line of the image. A line consist
      % of the orientation vectors and the projection
      % of the position onto the orientation vectors.            
      
      ori = self.getorientation();
      pos = self.getposition();
      
      linepos = pos*ori';
      
      line = [];
      line.orientation = ori;
      line.linepos = linepos;
    end
    
    %---------
    function images = getimages(self)
    %---------
      % Return the images contained in the dicom file.
      % The return data will be a struct with the fields
      % 'im' - The pixeldata as NxM single matrix
      % 'spacetimepos' - the position of the image in space and time
      % 'triggertime' - The trigger time of the image
      % 'instancenumber' - The instance number of the image
      % 'multiframenumber' - The position in PixelData buffer of 
      % the image.
      
      % Get nFrames
      if isfield(self.tags, 'NumberOfFrames')
        nFrames = segdicomtags.parsefloatstr(self.tags.NumberOfFrames);
      else
        nFrames = 1;
      end
      
      %self.tags.ImagePositionPatient
      
      % Check if it was a single frame
      if nFrames == 1
        images = self.getimagessingleframe();
        %disp(sprintf('TT:%0.5g',images.triggertime));
        return
      end

      % Check that FrameIncrementPointer exist
      if not(isfield(self.tags, 'FrameIncrementPointer'))
        disp('Warning: MultiFrame image is missing FrameIncrementPointer');
        frameincrementpointer = [];
      else
        frameincrementpointer = self.tags.FrameIncrementPointer;
      end
      
      % SliceVector
      if isequal(frameincrementpointer, [84 0 128 0])
        images = self.getimagesmultislice(nFrames);
        return
      end
      
      % RRVector, TimeSlotVector, SliceVector
      if isequal(frameincrementpointer, ...
          [84 0 96 0 84 0 112 0 84 0 128 0])
        images = self.getimagesmultirrtimeslice(nFrames);
        return
      end
      
      % Use FrameTime as default
      images = self.getimagesframetime(nFrames);
    end
    
    %---------
    function images = getimagessingleframe(self)
    %---------
      % Used by get images when there is only 
      % one frame.
      
      % Make sure Rows and Columns are set
      if not(isfield(self.tags, 'Rows')) || ...
          not(isfield(self.tags, 'Columns'))
        error('SEGMENT:ERROR', 'Missing Rows or Columns');
      end
      
      % Make sure we have the correct number of pixels
      pixelData = self.getpixeldata();
      cols = segdicomtags.parseuint16(self.tags.Columns);
      rows = segdicomtags.parseuint16(self.tags.Rows);
      if numel(pixelData) ~= cols*rows
        self.tags
        error('SEGMENT:ERROR', 'Wrong number of pixels');
      end

      % Get the pixeldata, reshape it and then permute x and y
      pixelData = reshape(pixelData, [cols rows])';
      
      % Make the images
      images = [];
      images.im = pixelData(:, :);
      images.spacetimepos = [ ...
        self.getposition()'; ...
        0; ...
        self.getvencpos()];
      images.triggertime = self.gettriggertime();
      images.acquisitiontime = self.getacquisitiontime();
      images.echotime = self.getechotime();
      images.inversiontime = self.getinversiontime(); %SBT151116
      images.sliceposition = self.getsliceposition(); %SBT151116
      images.instancenumber = self.getinstancenumber();
      images.multiframenumber = 1;
    end
    
    %---------
    function images = getimagesmultislice(self, nFrames)
    %---------
      % Used by getimages as a default multiframe parser.
      
      % Make sure Rows and Columns are set
      if not(isfield(self.tags, 'Rows')) || ...
          not(isfield(self.tags, 'Columns'))
        error('SEGMENT:ERROR', 'Missing Rows or Columns');
      end

      % Get the pixeldata, reshape it and then permute x and y
      pixelData = permute(reshape(self.getpixeldata(), ...
        [segdicomtags.parseuint16(self.tags.Columns) ...
        segdicomtags.parseuint16(self.tags.Rows) ...
        nFrames]), [2 1 3]);
      
      % Get the SliceVector and SpacingBetweenSlices
      if not(isfield(self.tags, 'SliceVector')) || ...
          not(isfield(self.tags, 'SpacingBetweenSlices'))
        error('SEGMENT:ERROR', ...
          'Missing SliceVector or SpacingBetweenSlices');
      end
      slicevector = typecast(self.tags.SliceVector, 'uint16');
      spacingbetweenslices = segdicomtags.parsefloatstr( ...
        self.tags.SpacingBetweenSlices);
      slicedist = self.getnormal()*spacingbetweenslices;
      
      % Make the images
      images = repmat(struct('im', {}, 'spacetimepos', {}, ...
        'triggertime', {}, 'instancenumber', {}), 1, nFrames);
      for i=1:nFrames
        images(i).im = pixelData(:, :, i);
        images(i).spacetimepos = [ ...
          self.getposition()' + slicedist'*double((slicevector(i)-1)); ...
          0; ...
          self.getvencpos()];
        images(i).triggertime = self.gettriggertime();
        images(i).acquisitiontime = self.getacquisitiontime();
        images(i).echotime = self.getechotime();
        images(i).inversiontime = self.getinversiontime();
        images(i).sliceposition = self.getsliceposition(); %SBT151116
        images(i).instancenumber = self.getinstancenumber();
        images(i).multiframenumber = 1;
      end
    end

    %---------
    function images = getimagesmultirrtimeslice(self, nFrames)
    %---------
      % Used by get images when numberofframes are greater
      % then one and the frame incremental pointer
      % is set to rr timeslice.
      
      % Make sure Rows and Columns are set
      if not(isfield(self.tags, 'Rows')) || ...
          not(isfield(self.tags, 'Columns'))
        error('SEGMENT:ERROR', 'Missing Rows or Columns');
      end

      % Get the pixeldata, reshape it and then permute x and y
      pixelData = permute(reshape(self.getpixeldata(), ...
        [segdicomtags.parseuint16(self.tags.Columns) ...
        segdicomtags.parseuint16(self.tags.Rows) ...
        nFrames]), [2 1 3]);
      
      % Get the RRIntervalVector, TimeSlotVector, SliceVector and
      % SpacingBetweenSlices
      if not(isfield(self.tags, 'RRIntervalVector')) || ...
          not(isfield(self.tags, 'TimeSlotVector')) || ...
          not(isfield(self.tags, 'SliceVector')) || ...
          not(isfield(self.tags, 'SpacingBetweenSlices'))
        self.tags
        error('SEGMENT:ERROR', ...
          'Missing some vector or SpacingBetweenSlices');
      end
      rrintervalvector = ...
        typecast(self.tags.RRIntervalVector, 'uint16');
      timeslotvector = typecast(self.tags.TimeSlotVector, 'uint16');
      slicevector = typecast(self.tags.SliceVector, 'uint16');
      spacingbetweenslices = segdicomtags.parsefloatstr( ...
        self.tags.SpacingBetweenSlices);
      slicedist = self.getnormal()*spacingbetweenslices;
      
      % Check that rrintervalvector isn't used
      if not(all(rrintervalvector == 1))
        error('SEGMENT:ERROR', 'Unsupported use of RRIntervalVector');
      end
      
      % Make the images
      images = repmat(struct('im', {}, 'spacetimepos', {}, ...
        'triggertime', {}, 'instancenumber', {}), 1, nFrames);
      for i=1:nFrames
        images(i).im = pixelData(:, :, i);
        images(i).spacetimepos = [ ...
          self.getposition()' + slicedist'*double((slicevector(i)-1)); ...
          0; ...
          self.getvencpos()];
        images(i).triggertime = self.gettriggertime();
        images(i).acquisitiontime = self.getacquisitiontime();        
        images(i).echotime = self.getechotime();
        images(i).inversiontime = self.getinversiontime();
        images(i).sliceposition = self.getsliceposition(); %SBT151116
        images(i).instancenumber = self.getinstancenumber();
        images(i).multiframenumber = double(timeslotvector(i));
      end
    end
    
    %---------
    function images = getimagesframetime(self, nFrames)
    %---------
      % Used by get images when numberofframes are greater
      % then one and the frame incremental pointer
      % is set to frametime.
      
      % Make sure Rows and Columns are set
      if not(isfield(self.tags, 'Rows')) || ...
          not(isfield(self.tags, 'Columns'))
        error('SEGMENT:ERROR', 'Missing Rows or Columns');
      end

      % Get the pixeldata, reshape it and then permute x and y
      pixelData = permute(reshape(self.getpixeldata(), ...
        [segdicomtags.parseuint16(self.tags.Columns) ...
        segdicomtags.parseuint16(self.tags.Rows) ...
        nFrames]), [2 1 3]);
      
      % Make the images
      images = repmat(struct('im', {}, 'spacetimepos', {}, ...
        'triggertime', {}, 'instancenumber', {}), 1, nFrames);
      for i=1:nFrames
        images(i).im = pixelData(:, :, i);
        images(i).spacetimepos = [ ...
          self.getposition()'; ...
          0; ...
          self.getvencpos()];
        images(i).triggertime = self.gettriggertime();
        images(i).acquisitiontime = self.getacquisitiontime();
        images(i).echotime = self.getechotime();
        images(i).inversiontime = self.getinversiontime();
        images(i).sliceposition = self.getsliceposition(); %SBT151116
        images(i).instancenumber = self.getinstancenumber();
        images(i).multiframenumber = i;
      end
    end

    %---------
    function r = gettriggertime(self)
    %---------
      % Returns the trigger time of the image
      
      r = 0; % Default value
   
      %---------------------------------------------------------
      function r = readttdata(tt)
      %---------------------------------------------------------
        % Auxiliary function to read TriggerTime data. Internal to
        % gettriggertime.
        tryasstring = sscanf(char(tt), '%f'); % try to read as string
        
        if length(tryasstring) ~= 1 % we got nothing, or several numbers
          if isa(tt, 'uint8') && length(tt) == 8 % might be a double
            r = typecast(tt, 'double');
          else
            r = 0;
          end
        else
          r = tryasstring;
        end
      end
      %---------------------------------------------------------
      
      if isfield(self.tags, 'TriggerTime') && ...
          ~isempty(self.tags.TriggerTime)

        r = readttdata(self.tags.TriggerTime);
        return
      end
      
      if isfield(self.tags,'CardiacTriggerDelayTime') && ...
          ~isempty(self.tags.CardiacTriggerDelayTime)

        r = readttdata(self.tags.CardiacTriggerDelayTime);
        return
      end
      
    end
    
    %---------
    function r = getinstancenumber(self)
    %---------
      % Returns the instance number of the dicom
      
      if isfield(self.tags, 'InstanceNumber')
        r = segdicomtags.parsefloatstr(self.tags.InstanceNumber);
      else
        r = 0;
      end
    end
    
    %---------
    function r = getvencpos(self)
    %---------
      % Returns the vencpos.
      %
      % Magnitude images get number 1, through-plane flow gets number 2.
      % The other directions are more tricky - Philips and Siemens do it
      % differently.
      %
      % Velocity directions:
      %
      % PhilipsVENC       Siemens         Flow    Segment
      % example           SequenceName    dir.    vencpos
      % ------------------------------------------------------
      % [100 0   0  ]     ends in 'rl'    RL/L    2
      % [0   100 0  ]     ends in 'ap'    AP/P    3
      % [0   0   100]     ends in 'fh     FH/S    4
      %    -------        ends in 'in'    (*)     (*)
      %
      % *: Siemens datasets so far contain 'in', 'ap' and 'fh' stacks.
      %    The 'in' stack is just 'through-plane', so we have to figure
      %    out which direction it is by looking at ImageOrientation.
      switch self.gettype()
        case 'mag'
          r = 1;
          return

        case 'phase'
          switch self.getscanner()
            case 'Philips'
              if isfield(self.tags, 'PhilipsVENC')
                tt = typecast(self.tags.PhilipsVENC, 'single');
                if numel(tt) ~= 3
                  error('SEGMENT:ERROR', 'Malformed PhilipsVENC tag - (DICOM 2001 101a)');
                end
                r = find(tt, 1) + 1;
                if isempty(r)
                  r = 2;
                end
              else
                r = 2;
              end
            case 'Siemens'
              seqname = strtrim(char(lower(self.tags.SequenceName)));

              if length(seqname) >= 2
                switch seqname(end-1:end)
                  case 'rl'
                    r = 2;
                  case 'ap'
                    r = 3;
                  case 'fh'
                    r = 4;
                  case 'in'
                    % Get through-plane direction in LPS
                    ori = self.getorientation();
                    tp = cross(ori(1,:), ori(2,:));
                    % The largest element in abs(tp) determines
                    % the flow direction.
                    [~, maxind] = max(abs(tp));
                    r = maxind + 1;
                  otherwise
                    r = 2;
                end
              else
                r = 2;
              end

            % We don't know how to handle multidirectional flow
            % in the other scanners.
            case 'GE'
              r = 2;
            case 'TOSHIBA_MEC'
              r = 2;
            otherwise
              r = 2;
          end
          return

        case 'psir'
          r = 3;
          return
        case 'angio'
          r = 5;
          return
          
        otherwise
          error('SEGMENT:ERROR', dprintf('Unknown image type ''%s''',self.gettype())); 
      end
    end

    %---------
    function r = getphotometricinterpretation(self)
    %---------
    % Returns the photometric interpretation
      if isfield(self.tags, 'PhotometricInterpretation')
        r = char(self.tags.PhotometricInterpretation);
      else
        r = 'MONOCHROME2 ';
      end
      if isequal(r, 'MONOCHROME1 ')
        r = 'MONOCHROME2 ';
      end
    end
    
    %---------
    function pixelData = getpixeldata(self)
    %---------
      % Parses and returns the PixelData field.
      
      % Make sure fields are present
      if not(isfield(self.tags, 'BitsAllocated')) || ...
          not(isfield(self.tags, 'BitsStored')) || ...
          not(isfield(self.tags, 'HighBit')) || ...
          not(isfield(self.tags, 'PixelData')) || ...
          not(isfield(self.tags, 'PixelRepresentation'))
        self.tags
        error('SEGMENT:ERROR', 'Missing fields for reading pixelData');
      end
      
      % Read pixelData
      switch segdicomtags.parseuint16(self.tags.BitsAllocated)
        case 8
          pixelData = typecast(self.tags.PixelData, 'uint8');
        case 16
          pixelData = typecast(self.tags.PixelData, 'uint16');
        otherwise
          self.tags
          error('SEGMENT:ERROR', 'BitsAllocated is not 8 or 16');
      end
      
      % Fix if there are unused bits
      bitsStored = segdicomtags.parseuint16(self.tags.BitsStored);
      highBit = segdicomtags.parseuint16(self.tags.HighBit);
      pixelData = bitshift(pixelData, -(highBit+1 - bitsStored));
      pixelData = bitand(pixelData, 2^bitsStored-1);
      
      % Fix if pixelRepresentation = 1 (signed)
      pixelData = single(pixelData);
      pixelRepresentation = segdicomtags.parseuint16( ...
        self.tags.PixelRepresentation);
      if pixelRepresentation == 1
        signedPixels = (pixelData >= 2^(bitsStored-1));
        pixelData(signedPixels) = pixelData(signedPixels)-2^bitsStored;
      end
      
      % Rescale according to RescaleSlope and RescaleIntercept
      if isfield(self.tags, 'RescaleSlope') && ~isempty(self.tags.RescaleSlope)
        pixelData = pixelData * ...
          segdicomtags.parsefloatstr(self.tags.RescaleSlope);
      end
      if isfield(self.tags, 'RescaleIntercept') && ~isempty(self.tags.RescaleIntercept)
        pixelData = pixelData + ...
          segdicomtags.parsefloatstr(self.tags.RescaleIntercept);
      end
      
      % Set to int16 for CT images.
      if strcmp(self.getmodality,'CT')
        pixelData = int16(pixelData);
      end
      
      % Fix special values for GE
      if strcmp(self.getscanner(), 'GE')
        pixelData(pixelData == -32768) = 0;
      end
      
      % Fix RGB -> grayscale
      switch self.getphotometricinterpretation()
        case 'RGB '
          if mod(numel(pixelData), 3) ~= 0
            self.tags
            mywarning('PhotometricInterpretation is RGB but number of pixels is not divisible by 3');
          end
          pixelData = conv(pixelData, [1 1 1]/3, 'valid');
          pixelData = pixelData(1:3:end);
        case {'MONOCHROME2 ','MONOCHROME 2','MONOCHROME2','PALETTE COLOR','PALETTE COLOR '}
          % nothing need to be done in this case.
        otherwise
          error('SEGMENT:ERROR', ...
            'Unsupported PhotometricInterpretation: %s', ...
            self.getphotometricinterpretation());
      end
    end
    
    %---------
    function r = getslicethickness(self)
    %---------
      % Returns slicethickness
      
      if not(isfield(self.tags, 'SliceThickness'))
        r = 0;
      else
        r = segdicomtags.parsefloatstr(self.tags.SliceThickness);
      end
    end
    
    %---------
    function r = getresolutionx(self)
    %---------
      % Returns the resolution in the x direction
      
      % Try the PixelSpacing element
      if isfield(self.tags, 'PixelSpacing')
        r = sscanf(char(self.tags.PixelSpacing), '%f\\%f');
        if not(numel(r) == 2)
          self.tags
          error('SEGMENT:ERROR', 'Couldn''t parse PixelSpacing');
        end
        r = r(1);
        return
      end
    
      r = 0; %It is very dangerous to return 1
    end
    
    %---------
    function r = getresolutiony(self)
    %---------
      % Returns the resolution in the y direction
      
      % Try the PixelSpacing element
      if isfield(self.tags, 'PixelSpacing')
        r = sscanf(char(self.tags.PixelSpacing), '%f\\%f');
        if not(numel(r) == 2)
          error('SEGMENT:ERROR', 'Couldn''t parse PixelSpacing');
        end
        r = r(2);
        return
      end
    
      r = 0; %It is very dangerous to return 1
    end
    
    %---------
    function r = getechotime(self)
    %---------
      % Returns the Echo Time
      
      if isfield(self.tags, 'EchoTime')
        r = segdicomtags.parsefloatstr(self.tags.EchoTime);
      else
        r = 0;
      end
    end
    
    %Inversiontime here /SBT151115
    %---------
    function r = getinversiontime(self)
    %---------
      % Returns the Inversion Time
      
      if isfield(self.tags, 'InversionTime')
        r = segdicomtags.parsefloatstr(self.tags.InversionTime);
      else
        r = 0;
      end
    end
    
    
    %---------
    function r = getrepetitiontime(self)
    %---------
      % Returns the RepetitionTime
      
      if isfield(self.tags, 'RepetitionTime')
        r = segdicomtags.parsefloatstr(self.tags.RepetitionTime);
      else
        r = 0;
      end
    end

    %---------
    function r = getaccessionnumber(self)
    %---------
      % Returns the accessionnumber field
      
      if isfield(self.tags,'AccessionNumber')
        r = char(self.tags.AccessionNumber);
      else
        r  = '';
      end;
    end;
 
    %---------
    function r = getstudyid(self)
    %---------
      % Returns study id field
      
      if isfield(self.tags,'StudyID')
        r = char(self.tags.StudyID);
      else
        r = '';
      end;
    end;
    
    %---------
    function r = getflipangle(self)
    %---------
      % Returns the FlipAngle field
      
      if isfield(self.tags, 'FlipAngle')
        r = segdicomtags.parsefloatstr(self.tags.FlipAngle);
      else
        r = 0;
      end
    end
    
    %---------
    function r = getnumberofaverages(self)
    %---------
      % Returns the NumberOfAverages field
      
      if isfield(self.tags, 'NumberOfAverages')
        r = segdicomtags.parsefloatstr(self.tags.NumberOfAverages);
      else
        r = 1;
      end
    end
    
    %---------
    function r = getseriesdescription(self)
    %---------
      % Returns the SeriesDescription
      
      if isfield(self.tags, 'SeriesDescription')
        r = char(self.tags.SeriesDescription);
      else
        r = '';
      end
    end
    
    %---------
    function r = getacquisitiontime(self)
    %---------
      % Returns the Acquisition Time
      
      if isfield(self.tags, 'AcquisitionTime')
        r = segdicomtags.parsetime(self.tags.AcquisitionTime);
      else
        r = 0;
      end
      
      %Try with ContentTime instead
      if (r==0) && isfield(self.tags, 'ContentTime')
        r = segdicomtags.parsetime(self.tags.ContentTime);
      end;
      
      %r = 0; %OBS only for debugging! Delete!
      
    end
    
    %---------
    function r = getseriesnumber(self)
    %---------
      % Returns the Series Number
      
      if isfield(self.tags, 'SeriesNumber')
        r = segdicomtags.parsefloatstr(self.tags.SeriesNumber);
      else
        r = 0;
      end
    end
    
    %---------
    function r = getimagetype(self)
    %---------
      % Returns the image type
      
      if isfield(self.tags, 'ImageType')
        r = char(self.tags.ImageType);
      else
        r = '';
      end
    end
    
    %---------
    function r = getstudyuid(self)
    %---------
      % Returns the study instance uid
      
      if isfield(self.tags, 'StudyInstanceUID')
        r = char(self.tags.StudyInstanceUID);
      else
        r = '';
      end
    end
    
    %---------
    function r = getsequencename(self)
    %---------
      % Returns the Sequence Name
      
      if isfield(self.tags, 'SequenceName')
        r = char(self.tags.SequenceName);
      else
        r = '';
      end
    end
    
    %---------
    function r = getbitsstored(self)
    %---------
      % Returns the BitsStored field
      
      if isfield(self.tags, 'BitsStored')
        r = segdicomtags.parseuint16(self.tags.BitsStored);
      else
        error('SEGMENT:ERROR', 'BitsStored is missing');
      end
    end
    
    %---------
    function r = getmodality(self)
    %---------
      % Returns the Modality
      
      if isfield(self.tags, 'Modality')
        r = char(self.tags.Modality);
      else
        r = 'MR';
      end
    end
    
    %---------
    function r = getpatientinfo(self)
    %---------
      % Returns a Patientinfo struct with fields
      % 'PatientName'
      % 'PatientID'
      % 'PatientBirthDate'
      % 'PatientSex'
      % 'PatientAge'
      % 'AcquisitionDate'
      % 'PatientWeight'
      
      r = [];
      if isfield(self.tags, 'PatientName') && ...
          numel(self.tags.PatientName) > 0
        r.Name = char(self.tags.PatientName);
      else
        r.Name = '';
      end
      if isfield(self.tags, 'PatientID') && ...
          numel(self.tags.PatientID) > 0
        r.ID = char(self.tags.PatientID);
      else
        r.ID = [];
      end
      if isfield(self.tags, 'PatientBirthDate') && ...
          numel(self.tags.PatientBirthDate) > 0
        r.BirthDate = char(self.tags.PatientBirthDate);
      else
        r.BirthDate = [];
      end
      if isfield(self.tags, 'PatientSex') && ...
          numel(self.tags.PatientSex) > 0
        r.Sex = char(self.tags.PatientSex);
      else
        r.Sex = [];
      end
      if isfield(self.tags, 'PatientAge') && ...
          numel(self.tags.PatientAge) > 0
        r.Age = double(self.tags.PatientAge);
        if numel(self.tags.PatientAge) > 1
          age = char(self.tags.PatientAge);
          r.Age = str2double(age(isstrprop(age,'digit'))); %Remove nondigits
        end
      else
        r.Age = [];
      end
            
      if isfield(self.tags, 'AcquisitionDate') && ...
          numel(self.tags.AcquisitionDate) > 0
        r.AcquisitionDate = char(self.tags.AcquisitionDate);
      else
        r.AcquisitionDate = [];
      end      
      %Try with Content date instead
      if isempty(r.AcquisitionDate) && isfield(self.tags, 'ContentDate') && ...
          numel(self.tags.ContentDate) > 0
        r.AcquisitionDate = char(self.tags.ContentDate);
      end;
      
      if isfield(self.tags, 'PatientSize') && ...
          numel(self.tags.PatientSize) > 0
        try
          r.Length = 100*segdicomtags.parsefloatstr(self.tags.PatientSize); %convert from m to cm
        catch
          r.Length = [];
        end
      else
        r.Length = [];
      end
      if isfield(self.tags, 'PatientWeight') && ...
          numel(self.tags.PatientWeight) > 0
        try
          r.Weight = segdicomtags.parsefloatstr(self.tags.PatientWeight);
        catch
          r.Weight = [];
        end
      else
        r.Weight = [];
      end
      if ~isempty(r.Length) && ~isempty(r.Weight)
        r.BSA = calcfunctions('calcbsa',r.Weight,r.Length);
      else
        r.BSA = [];
      end
    end
    
    %---------
    function r = getheartrate(self)
    %---------
    %Return heart rate
    if isfield(self.tags, 'HeartRate')
      r = segdicomtags.parsefloatstr(self.tags.HeartRate);
    else
      r = 0;
    end
    end
    
    %---------
    function r = getvenc(self)
    %---------
      % Returns the venc
      
      switch(self.getscanner())
        case 'GE'
          if isfield(self.tags, 'VelocityEncoding')
            r = segdicomtags.parseuint16(self.tags.VelocityEncoding)/10;
          else
            r = 0;
          end
        case 'Philips'
          if isfield(self.tags, 'PhilipsVENC')
            philipsVENC = typecast(self.tags.PhilipsVENC, 'single');
            r = max(philipsVENC);
          else
            r = 0;
          end
        case 'Siemens'
          % For Siemens, the VENC is stored in text in SequenceName.
          % Get sequenceName and check that it's longer than 2.
          if not(isfield(self.tags, 'SequenceName'))
            r = 0;
            return
          end
          sequenceName = lower(char(self.tags.SequenceName));
          if length(sequenceName) <= 2
            r = 0;
            return
          end

          % Find VENC in sequenceName. 'pos' is the position two characters
          % *before* the VENC string.
          pos = [];
          if not(isempty(findstr(sequenceName, 'fl')))
            pos = findstr(sequenceName,'_v');
          end
          if isempty(pos) && sequenceName(1)=='v'
            pos=0;
          end
          if isempty(pos)
            if (sequenceName(1)=='p') && not(isempty(strfind(sequenceName,'cms')))
              sequenceName = ['..' segdicomtags.removechars(sequenceName)];
              pos = 1;
            end
          end

          if isempty(pos) % Siemens 4D flow data
            if sequenceName(1) == 'p' && not(isempty(strfind(sequenceName, '_v')))
              pos = strfind(sequenceName, '_v');
            end
          end

          % Return if we didn't find pos
          if isempty(pos)
            r = 0;
            return
          end

          % Get the VENC
          sequenceName = ...
            segdicomtags.removechars(sequenceName((pos+2):end));
          [r, ok] = str2num(sequenceName); %#ok<ST2NM>
          if not(ok)
            r = 0;
          end
        case 'TOSHIBA_MEC'
          if isfield(self.tags, 'ToshibaVENCMaxVal')
            toshibaVENCmax = typecast(self.tags.ToshibaVENCMaxVal, 'single');
            toshibaVENCmin = typecast(self.tags.ToshibaVENCMaxVal, 'single');
            r = max(toshibaVENCmax) - max(toshibaVENCmin);
          else
            r = 0;
          end
        otherwise
          r = 0;
      end
    end
    
    %---------
    function r = hassegmentdata(self)
    %---------
      % Returns true if the tag SegmentData is present
      
      r = isfield(self.tags, 'SegmentData');
    end
    
    %---------
    function r = getsegmentdata(self)
    %---------
      % Parses and returns the SegmentData
      
      r = segdicomfile.unserialize(self.tags.SegmentData);
    end
    
    %---------
    function r = hastriggertime(self)
    %---------
      % Returns true if has trigger time information
      
      r = isfield(self.tags, 'TriggerTime') || isfield(self.tags, 'CardiacTriggerDelayTime');
    end
    
    %---------
    function r = getspectspecialtag(self)
    %---------
    % Returns a special tag used in SPECT images
      r = [];
      if isfield(self.tags, 'UnknownTag269615121')
        r = char(self.tags.UnknownTag269615121);
      end
    end
    
    %---------
    function r = hasrepetitiontime(self)
    %---------
      % Returns the Repetition Time
      
      r = isfield(self.tags, 'RepetitionTime');
    end
    %---------
    function r = hasvelocityencodescale(self)
    %---------
      % Returns true if VelocityEncodeScale field is set
      
      r = isfield(self.tags, 'VelocityEncodeScale');
    end
    %---------
    function r = getvelocityencodescale(self)
    %---------
      % Returns the VelocityEncodeScale field
      
      if not(self.hasvelocityencodescale())
        error('SEGMENT:PANIC', 'No VelocityEncodeScale field');
      end
      r = segdicomtags.parsefloatstr(self.tags.VelocityEncodeScale);
    end
    
    %---------
    function r = unpack(self)
    %---------
    %Unpack data from dicom file
      r = self;
      
      % Check nFrames
      if not(isfield(self.tags, 'NumberOfFrames'))
        return
      end
      nFrames = segdicomtags.parsefloatstr(self.tags.NumberOfFrames);      
      if nFrames == 1
        return
      end
      
      % Check the frameinc isn't set
      if isfield(self.tags, 'FrameIncrementPointer')
        return
      end
      
      % check that FunctionalGroupsSequence's exist
      if not(isfield(self.tags, 'SharedFunctionalGroupsSequence'))
        return
      end
      if not(isfield(self.tags, 'PerFrameFunctionalGroupsSequence'))
        return
      end
      
      %Define groupnames we need to parse deeper into
      groupnames = {...
        'CardiacTriggerSequence', ... 
        'FrameContentSequence', ...        
        'MRImageFrameTypeSequence', ...              
        'MRTimingAndRelatedParametersSeq', ...                
        'PlaneOrientationSequence', ...
        'PixelMeasuresSequence', ...        
        'PlanePositionSequence', ...        
        'PixelValueTransformationSequence', ...
        'RealWorldValueMappingSequence' ,...
        'UnknownTag336535557'}; %This is a Philips special!! Here all the details were included.

      %If think we can ignore these as we are not using information inside
      %them. /EH      
      %  'FrameAnatomySequence', ...
      %  'FrameVOILUTSequence', ...
      %  'MRFOV-GeometrySequence', ...        
      %  'MRReceiveCoilSequence', ...
      %  'MRTransmitCoilSequence', ...
      %  'MRMetaboliteMapSequence', ...
      %  'MRImagingModifierSequence', ...
      %  'MREchoSequence', ...
      %  'MRAveragesSequence', ...
      %  'MRModifierSequence', ...        
      
      % Check the groups if special interpretation
      needspecial = false;
      needspecial = needspecial | self.checkFunctionalGroup(...
        self.tags.SharedFunctionalGroupsSequence, groupnames);
      needspecial = needspecial | self.checkFunctionalGroup(...
        self.tags.PerFrameFunctionalGroupsSequence, groupnames);           
      
      % Unpack the dicom
      if ~needspecial
        %Jonatan Wulcan standard code
        r = [];
        pixelsPerFrame = numel(self.tags.PixelData)/nFrames;
        for i=1:nFrames
          frameTags = self.tags;
          frameTags = rmfield(frameTags, 'SharedFunctionalGroupsSequence');
          frameTags = rmfield(frameTags, 'PerFrameFunctionalGroupsSequence');
          frameTags.PixelData = frameTags.PixelData(pixelsPerFrame*(i-1) + 1:pixelsPerFrame*i);
          frameTags.NumberOfFrames = '1';
          for j=1:numel(groupnames)
            if isfield(self.tags.SharedFunctionalGroupsSequence{1}, groupnames{j})
              frameTags = self.addstructs(...
                frameTags, ...
                self.tags.SharedFunctionalGroupsSequence{1}.( groupnames{j} ){1});
            end
            
            if isfield(self.tags.PerFrameFunctionalGroupsSequence{i}, groupnames{j})
              frameTags = self.addstructs(...
                frameTags, ...
                self.tags.PerFrameFunctionalGroupsSequence{i}.( groupnames{j} ){1});
            end
            
          end
               
          %Convert to segdicomtags class and store to temporary variable
          temptags = segdicomtags(frameTags);
          
          %Some tags are called differently if they are frame base in
          %PerFrameFunctionalGroupsSequence compare to single DICOMs.
          %Simply rename them.
          %if isfield(temptags.tags,'TriggerDelayTime')
          %  temptags.tags.TriggerTime = temptags.tags.TriggerDelayTime; %Copy field info
          %end;
          %if isfield(temptags.tags,'FrameType')
          %  temptags.tags.ImageType= temptags.tags.FrameType; %Copy field info
          %end;
          
          r = [r temptags];
          
          %This code should work, but does not, should be faster than the
          %above, but I am leaving that for now /EH:
          %if i==1
          %  r = segdicomtags(frameTags);
          %  r(nFrames) = r(1);
          %else
          %  r(i) = segdicomtags(frameTags);
          %end;          
          
        end
        
      else
        %%% Special unpack Einar Heiberg, PhilipsEnhancedMR
        
        tags2find = [];
        tags2find(1).code = '0020,0032';  tags2find(1).name =  'ImagePositionPatient'; tags2find(1).vr = 'DS'; %0020,0032
        tags2find(2).code = '0020,0037';  tags2find(2).name =  'ImageOrientationPatient'; tags2find(2).vr = 'DS'; %0020,0037       
        tags2find(3).code = '0020,9153';  tags2find(3).name =  'CardiacTriggerDelayTime'; tags2find(3).vr = 'UN'; %0020,9153
        tags2find(4).code = '0008,0008';  tags2find(4).name =  'ImageType'; tags2find(4).vr = 'CS'; %0008,0008
        tags2find(5).code = '2001,0008';  tags2find(5).name =  'PhilipsVENC'; tags2find(5).vr = 'FL'; %2001,0008        
        tags2find(6).code = '0018,1060';  tags2find(6).name =  'TriggerTime'; tags2find(6).vr = 'UN';
        tags2find(7).code = '0018,0081';  tags2find(7).name =  'EchoTime'; tags2find(7).vr ='UN';
        tags2find(8).code = '0018,1088';  tags2find(8).name =  'HeartRate'; tags2find(8).vr ='UN'; 
        tags2find(9).code = '0020,0013';  tags2find(9).name =  'InstanceNumber'; tags2find(9).vr ='UN';
        tags2find(10).code = '0018,1314'; tags2find(10).name = 'FlipAngle'; tags2find(10).vr = 'UN';
        tags2find(11).code = '0018,0082'; tags2find(11).name = 'InversionTime'; tags2find(11).vr = 'UN';
        %tags2find(12).code = '0028,1052'; tags2find(12).name = 'RescaleIntercept'; tags2find(12).vr = 'UN'; %was 0028,1052
        %tags2find(13).code = '0028,1053'; tags2find(13).name = 'RescaleSlope'; tags2find(13).vr = 'UN';
        tags2find(12).code = '0018,9197'; tags2find(12).name = 'ToshibaVENCSeq'; tags2find(12).vr = 'SQ';
        tags2find(13).code = '0028,0030'; tags2find(13).name = 'PixelSpacing'; tags2find(13).vr = 'DS';       
       
        %Nested sequence groups to find
        nested2find = [];
        %nested2find(1).code = '0018,9090'; nested2find(1).name = 'ToshibaVENCDir'; nested2find(1).vr = 'UN';
        %nested2find(2).code = '0018,9091'; nested2find(2).name = 'ToshibaVENCMinVal'; nested2find(2).vr = 'UN';
        %nested2find(3).code = '0018,9217'; nested2find(3).name = 'ToshibaVENCMaxVal'; nested2find(3).vr = 'UN';

        tags2find(1).value = []; %Reservs memory for all positions...       
        %nested2find(1).value = []; %Reservs memory for all positions...

        %Convert code numbers that corresponds how tags are stored [32 0 55 0] = '0020,0037'
        for loop = 1:length(tags2find)
          stri = tags2find(loop).code;
          num1 = hex2dec(stri(1:4));
          num2 = hex2dec(stri(6:end));
          tags2find(loop).code = [typecast(uint16(num1),'uint8') typecast(uint16(num2),'uint8')]; %#ok<AGROW>
        end;
        
        for loop = 1:length(nested2find)
          stri = nested2find(loop).code;
          num1 = hex2dec(stri(1:4));
          num2 = hex2dec(stri(6:end));
          nested2find(loop).code = [typecast(uint16(num1),'uint8') typecast(uint16(num2),'uint8')]; %#ok<AGROW>
        end;        
        
        sequence = self.tags.PerFrameFunctionalGroupsSequence;
        %seguence = self.tags.SharedFunctionalGroupsSequence;
        
        %Check if implicit or explicit vr
        explicit = true;
        pos = strfind(sequence,uint8([24 0 24 145])); %(0018,9118) cardiacsyncronizationsequence
        if ~isempty(pos)
          pos = pos(1); %Take first
          if isequal(char(sequence(pos+4:pos+5)),'SQ')
            explicit = true;
          else
            explicit = false;
          end;
        end;
        
        %Call helper function to unpack
        tags2find = self.unpackhelper(tags2find,nested2find,sequence,explicit);
        
        %Create new dicom structs
        r = [];
        pixelsPerFrame = numel(self.tags.PixelData)/nFrames;
        for i=1:nFrames
          frameTags = self.tags;
          frameTags = rmfield(frameTags, 'SharedFunctionalGroupsSequence');
          frameTags = rmfield(frameTags, 'PerFrameFunctionalGroupsSequence');
          frameTags.PixelData = frameTags.PixelData(pixelsPerFrame*(i-1) + 1:pixelsPerFrame*i);
          frameTags.NumberOfFrames = '1';

          %Here I need to add extracted tags from above
          for j=1:numel(tags2find)
            
            %Calculate potential duplicate factor
            duplicatefactor = 1;
            numvalues = length(tags2find(j).value);
            if (~isequal(numvalues,nFrames)) && (numvalues>0)
              if numvalues<nFrames
                error('Too few PerFrame items.');
              end;
              if isequal(numvalues,nFrames*2)
                duplicatefactor = 2;
              end;
            end;
            
            if ~isempty(tags2find(j).value)
              if ~strcmp(tags2find(j).vr,'SQ')
                frameTags.( tags2find(j).name ) = tags2find(j).value{1+(i-1)*duplicatefactor};
              else
                tagnest = tags2find(j).value{i};
                for k = 1:numel(tagnest)
                  if ~isempty(tagnest(k).value)                  
                    frameTags.( tagnest(k).name ) = tagnest(k).value{1};
                  end
                end
              end
            end;
          end
          
          r = [r segdicomtags(frameTags)];
        end
        
      end
    end 
  end
  
  methods(Static)
    
    %---------
    function tags2find = unpackhelper(tags2find,nested2find,sequence,explicit)
    %---------
      %Helper function to unpack data. EH:
      
      %Loop over the tags to find
      for loop=1:length(tags2find)

        pos = findstr(sequence,uint8(tags2find(loop).code)); %#ok<FSTR>
          
        if ~isempty(pos)
        
          %Reserve memory
          tags2find(loop).value = cell(1,length(pos));
            
          %Loop over found tags
          for ploop = 1:length(pos)
              
            ppos = pos(ploop)+4; %(group,element)
              
            if explicit
              vr = char(sequence((ppos):(ppos+1))); %Need to use later if extract more.
              ppos = ppos+2;
            else
              vr = tags2find(loop).vr;
            end;
            
            switch vr
              case {'DS','IS'}
                %DS String with double(s)
                %IS String with double(s)
                taglength = double(typecast(sequence(ppos:ppos+1),'uint16'));
                if explicit
                  stri = char(sequence((ppos+2):(ppos+2+taglength-1)));
                else
                  stri = char(sequence((ppos+4):(ppos+4+taglength-1)));
                end;
                tags2find(loop).value{ploop} = stri; %Store
                %disp(sprintf('%s VR:%s Length:%d %s',tags2find(loop).name,vr,taglength,stri));
              case 'CS'
                taglength = double(typecast(sequence(ppos:ppos+1),'uint16'));
                if explicit
                  stri = char(sequence((ppos+2):(ppos+2+taglength-1)));
                else
                  stri = char(sequence((ppos+4):(ppos+4+taglength-1)));
                end;
                tags2find(loop).value{ploop} = stri; %Store
              case 'FD'
                %Added EH:
                taglength = double(typecast(sequence(ppos:ppos+1),'uint16'));
                if explicit
                  val = typecast(sequence((ppos+2):(ppos+2+taglength-1)),'DOUBLE');
                else
                  val = typecast(sequence((ppos+4):(ppos+4+taglength-1)),'DOUBLE');
                end;
                tags2find(loop).value{ploop} = val; %Store
              case 'UN'
                %Unknown
                if explicit
                  taglength = double(typecast(sequence(ppos+2:ppos+5),'uint32'));
                  val = sequence(ppos+6:ppos+6+taglength-1);
                else
                  taglength = double(typecast(sequence(ppos:ppos+3),'uint32'));
                  val = sequence(ppos+4:ppos+4+taglength-1);
                end;
                %disp(sprintf('%s VR:%s Length:%d Val:%s',tags2find(loop).name,vr,taglength,sprintf('%d ',val)));
                
                %Store it
                if isequal(length(val),8)
                  %Store, guess on double if 8 bits
                  tags2find(loop).value{ploop} = typecast(val,'double');
                  
                  %Check if reasonable number, otherwise go back
                  if (tags2find(loop).value{ploop}<0) || (tags2find(loop).value{ploop}>1e5)
                    tags2find(loop).value{ploop} = val;
                  end;
                else
                  tags2find(loop).value{ploop} = val; %Store
                end
              case 'SQ'
                %Aarhgg I hate nested sequence groups EH:
                
                if explicit
                  taglength = double(typecast(sequence(ppos+2:ppos+5),'uint32'));
                  val = sequence(ppos+6:ppos+6+taglength-1);
                else
                  taglength = double(typecast(sequence(ppos:ppos+3),'uint32'));
                  val = sequence(ppos+4:ppos+4+taglength-1);
                end;
                
                nestedtags = segdicomtags.unpackhelper(nested2find,[],val,explicit); %[] no nested nested tags please!
                tags2find(loop).value{ploop} = nestedtags;
              otherwise
                error('SEGMENT:ERROR',dprintf('Unknown VR %s',vr)); 
            end;
            
          end;

        end;
        
      end;
            
    end;
        
    %---------
    function a = addstructs(a, b)
    %---------
    % Add fields from struct b to struct a
    fnames = fieldnames(b);
      for i=1:numel(fnames)
        a.( fnames{i} ) = b.( fnames{i} );
      end
    end
    
    %---------
    function implicit = checkFunctionalGroup(group, groupnames)
    %---------
    % Check functional groups to see if they are implicit
      implicit = false;
      
      if not(isa(group, 'cell'))
        implicit = true;
        return;
      end
      
      for i=1:numel(group)
        for j=1:numel(groupnames)
          if not(isfield(group{i}, groupnames{j}))
            continue
          end
          if not(isa(group{i}.( groupnames{j} ), 'cell'))
            implicit = true;
            return;
          end
        end
      end
    end
    
    %---------
    function r = parseuint16(data)
    %--------- 
      % Parses a uint8 matrix with two elements as
      % a single uint16 number.
      
      r = double(typecast( ...
        data, 'uint16'));
      if numel(r) ~= 1
        data
        error('SEGMENT:ERROR', 'Couldn'' parse uint16');
      end
    end
    
    %---------
    function r = parsesingle(data)
    %---------  
      % Parses a uint8 matrix with four elements as
      % a single number of type single.
      
      r = double(typecast( ...
        data, 'single'));
      if numel(r) ~= 1
        data
        error('SEGMENT:ERROR', 'Couldn'' parse single');
      end
    end
    
    %---------
    function r = parsefloatstr(data)
    %---------
      % Parses a string with a number as a
      % number of type double.
      
      data = strtrim(char(data));
      
      if isempty(data) || isequal(data,0) %EH:
        r = 0; %EH:
      else
        r = sscanf(data, '%f');
        if numel(r) ~= 1
          data
          error('SEGMENT:ERROR', 'Couldn''t parse floatstr');
        end
      end
    end
    
    %---------
    function timenum = parsetime(timestr)
    %---------
      % Parse a dicom time string as number of seconds.
      
      % get timestr and reset
      timestr = char(timestr);
      timestr = timestr(regexp(timestr,'\d')); %include only numerics

      hour = 0;
      minute = 0;
      seconds = 0;
      partseconds = 0;

      % Exctract hour, minute, seconds and partseconds
      % Divide partseconds to make it fractions of a second.
      if length(timestr)>=2
        hour = str2double(timestr(1:2));
      end
      if length(timestr)>=4
        minute = str2double(timestr(3:4));
      end
      if length(timestr)>=6
        seconds = str2double(timestr(5:6));
      end
      if length(timestr)>6
        partseconds = str2double(timestr(7:end));
        partseconds = partseconds*10^(-(length(timestr)-6));
      end
      
      % Calc timenum
      timenum = hour*3600+minute*60+seconds+partseconds;
    end
    
    %---------
    function [s,ind] = removechars(stri)
    %---------
      % Removes everything except numbers and '.'
      % from a string.
      
      ind = (stri=='0');
      ind = ind|(stri=='1');
      ind = ind|(stri=='2');
      ind = ind|(stri=='3');
      ind = ind|(stri=='4');
      ind = ind|(stri=='5');
      ind = ind|(stri=='6');
      ind = ind|(stri=='7');
      ind = ind|(stri=='8');
      ind = ind|(stri=='9');
      ind = ind|(stri=='.');
      s = stri(ind);
    end
       
  end
  
  
  methods(Static) % Public static methods
    %---------
    function dicoms = readfiles(files)
    %---------
      % Returns a matrix of segdicomtags objects
      % contaion DICOM info from 'files'.
      
      segloaderprogressbar('update', ...
        struct('name', 'readfiles', 'numfiles', numel(files)));
      try
        dicomsraw = segdicomread_mex(files);
      catch e
        clear segdicomread_mex
        rethrow(e);
      end
      clear segdicomread_mex
      
      % Wrap the dicom struct in segdicomtags class
      dicoms = segdicomtags(dicomsraw{1}).unpack();
      for i=2:numel(dicomsraw)
        dicoms = [dicoms segdicomtags(dicomsraw{i}).unpack()];
      end
      clear dicomsraw
    end
  end
end

