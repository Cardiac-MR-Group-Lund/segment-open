classdef segdicomfile
  
  %Need to write:
  %0008,0020 = StudyDate 
  %0008,0023 = ContentDate (set to same as studydate)
  %0008,0030 = StudyTime (from AcquisitionTime)
  %0008,0033 = ContentTime (set to same as studytime)
  %0008,0050 = AccessionNumber (undersökningsnummer) OK.
  %0020,0010 = StudyID (Remissnummer)
  
  methods(Static) % Public
    function create(filename, data, study_uid, pat_name, pat_id, ...
                    pat_birth, pat_sex, switchtags)
      % Creats a dicom file. Serialize data arg and
      % store it in special segment tag.

      % Open file
      mem = filebasket(filename, 'w');
      onCleanup(@()mem.close());

      % Generate UID:s and load image
      instance_uid = segdicomfile.generate_uid();
      if isempty(study_uid)
        study_uid = segdicomfile.generate_uid();
      else
        study_uid = uint8(study_uid);
      end
      series_uid = segdicomfile.generate_uid();
      load('image.mat', 'image');

      % Write 128 zeros and DICM
      mem.add([uint8(zeros(1, 128)) uint8('DICM')]);

      % Write meta header
      mh = segdicomfile.create_metaheader(instance_uid);
      segdicomfile.write_tag(mem, 'FileMetaInfoGroupLength', 'UL', ...
        typecast(uint32(length(mh)), 'uint8') );
      mem.add(mh);

      % Write other tags
      segdicomfile.write_tag(mem, ...
        'SOPClassUID', 'UI', uint8('1.2.840.10008.5.1.4.1.1.7'));
      segdicomfile.write_tag(mem, 'SOPInstanceUID', 'UI', instance_uid);
      segdicomfile.write_tag(mem, 'StudyDate', 'DA', uint8(data(1).PatientInfo.AcquisitionDate));
      segdicomfile.write_tag(mem, 'ContentDate', 'DA', uint8(data(1).PatientInfo.AcquisitionDate));
      segdicomfile.write_tag(mem, 'StudyTime', 'TM', uint8(segdicomfile.secondtostring(data(1).AcquisitionTime)));
      segdicomfile.write_tag(mem, 'ContentTime', 'TM', uint8(segdicomfile.secondtostring(data(1).AcquisitionTime)));
      if switchtags
        segdicomfile.write_tag(mem, 'AccessionNumber', 'SH', uint8(data(1).StudyID)); %StudyID and AccessionNumber are switched in Sectra PACS
      else
       segdicomfile.write_tag(mem, 'AccessionNumber', 'SH', uint8(data(1).AccessionNumber)); %StudyID and AccessionNumber are switched in Sectra PACS
      end
      segdicomfile.write_tag(mem, 'Modality', 'CS', uint8('OT'));
      segdicomfile.write_tag(mem, 'ConversionType', 'CS', uint8('WSD'));
      segdicomfile.write_tag(mem, ...
        'ReferringPhysicianName', 'PN', uint8([]));
      segdicomfile.write_tag(mem, 'PatientName', 'PN', uint8(pat_name));
      segdicomfile.write_tag(mem, 'PatientID', 'LO', uint8(pat_id));
      segdicomfile.write_tag(mem, ...
        'PatientBirthDate', 'DA', uint8(pat_birth));
      segdicomfile.write_tag(mem, 'PatientSex', 'CS', uint8(pat_sex));
      segdicomfile.write_tag(mem, ...
        'SecondaryCaptureDeviceManufacturer', 'LO', uint8('Medviso '));
      segdicomfile.write_tag(mem, ...
        'SecondaryCaptureDeviceModelName', 'LO', uint8('Segment '));
      segdicomfile.write_tag(mem, 'StudyInstanceUID', 'UI', study_uid);
      segdicomfile.write_tag(mem, 'SeriesInstanceUID', 'UI', series_uid);
      if switchtags
        segdicomfile.write_tag(mem, 'StudyID', 'SH', uint8(data(1).AccessionNumber));
      else
        segdicomfile.write_tag(mem, 'StudyID', 'SH', uint8(data(1).StudyID));
      end
      segdicomfile.write_tag(mem, 'SeriesNumber', 'IS', uint8([]));
      segdicomfile.write_tag(mem, 'InstanceNumber', 'IS', uint8([]));
      segdicomfile.write_tag(mem, 'PatientOrientation', 'CS', uint8([]));
      segdicomfile.write_tag(mem, ...
        'SamplesPerPixel', 'US', typecast(uint16(1), 'uint8'));
      segdicomfile.write_tag(mem, ...
        'PhotometricInterpretation', 'CS', uint8('MONOCHROME2 '));
      segdicomfile.write_tag(mem, ...
        'Rows', 'US', typecast(uint16(128), 'uint8'));
      segdicomfile.write_tag(mem, ...
        'Columns', 'US', typecast(uint16(128), 'uint8'));
      segdicomfile.write_tag(mem, ...
        'BitsAllocated', 'US', typecast(uint16(16), 'uint8'));
      segdicomfile.write_tag(mem, ...
        'BitsStored', 'US', typecast(uint16(16), 'uint8'));
      segdicomfile.write_tag(mem, ...
        'HighBit', 'US', typecast(uint16(15), 'uint8'));
      segdicomfile.write_tag(mem, ...
        'PixelRepresentation', 'US', typecast(uint16(0), 'uint8'));
      segdicomfile.write_tag(mem, ...
        'SmallestImagePixelValue', 'US', typecast(uint16(0), 'uint8'));
      segdicomfile.write_tag(mem, ...
        'LargestImagePixelValue', 'US', typecast(uint16(65536), 'uint8'));
      segdicomfile.write_tag(mem, 'SegmentData', 'OB', ...
        segdicomfile.serialize(data));
      segdicomfile.write_tag(mem, ...
        'PixelData', 'OW', image);
    end
    
    function r = serialize( data )
      % Serializes a matlab variable. Return a [1xn] uint8 array.

      % Error check
      if issparse(data)
        error('SEGMENT:PANIC', ...
          'Serialize doesn''t support sparse matrices');
      end

      % Create data_size
      data_size = typecast(uint32(size(data)), 'uint8');

      % Create data_content
      switch class(data)
        case {'double', 'single', 'int8', 'uint8', ...
            'int16', 'uint16', 'int32', 'uint32', ...
            'int64', 'uint64'}
          if ~isreal(data)
            error('SEGMENT:PANIC', 'Serialize doesn''t support complex data');
          end
          data_content = typecast(data(:), 'uint8');
        case {'logical', 'char', 'cell', 'struct'}
          data_content = eval(...
            ['segdicomfile.serialize_' class(data) '(data)']);
        otherwise
          error('SEGMENT:PANIC', 'Unknown class');
      end

      % Pack it all togheter
      r = segdicomfile.create_chunk_va(...
        data_size, ...
        segdicomfile.serialize_char(class(data)), ...
        data_content);
    end

    function data = unserialize( r )
      % Unserialize [1xn] uint8 array r to a matlab variable

      % Parse r as chunk
      [data_size, data_class, data_content] = ...
        segdicomfile.parse_chunk_va(r);

      % Parse data_dims and data_class
      data_dims = typecast(data_size, 'uint32');
      data_class = segdicomfile.unserialize_char(data_class);

      % Parse elements
      switch data_class
        case {'double', 'single', 'int8', 'uint8', ...
            'int16', 'uint16', 'int32', 'uint32', ...
            'int64', 'uint64'}
          data = typecast(data_content, data_class);
        case {'logical', 'char', 'cell', 'struct'}
          data = eval(...
            ['segdicomfile.unserialize_' data_class '(data_content)']);
        otherwise
          error('SEGMENT:ERROR', ...
            'Invalid serialization stream: Unknown class');
      end

      % Reshape to correct size
      data = reshape(data, data_dims);
    end
  end
  
  methods(Static = true, Access = private)
    function r = create_metaheader(instance_uid)
      % Returns a memorybasket containing a dicom
      % metaheader.

      r = memorybasket();
      segdicomfile.write_tag(r, 'FileMetaInfoVersion', 'OB', uint8([0 1]));
      segdicomfile.write_tag(r, 'MediaStorageSOPClassUID', ...
        'UI', uint8('1.2.840.10008.5.1.4.1.1.7'));
      segdicomfile.write_tag(r, ...
        'MediaStorageSOPInstanceUID', 'UI', instance_uid);
      segdicomfile.write_tag(r, ...
        'TransferSyntaxUID', 'UI', uint8('1.2.840.10008.1.2.1'));
      segdicomfile.write_tag(r, 'ImplementationClassUID', ...
        'UI', uint8('1.3.6.1.4.1.9590.100.1.0.100.4.0'));
      segdicomfile.write_tag(r, 'ImplementationVersionName', ...
        'SH', uint8('SEGMENT COUNTOURDICOM 0.1'));
    end
    
    function write_tag(mem, tag_name, vr, data)
      % Writes a dicom tag (little endian explicit VR
      % style) to a basket. Note that this function
      % only support a subset of VR:s.
      segdicomfile.name_to_tag(tag_name);

      switch vr
        case {'OB', 'OW'}
          mem.add([segdicomfile.name_to_tag(tag_name) uint8(vr) 0 0 ...
            typecast(uint32(length(data)+mod(length(data), 2)), 'uint8')]);
          mem.add(data);
          if mod(length(data), 2) ~= 0
            mem.add(uint8(0));
          end
        case {'UI'}
          mem.add([segdicomfile.name_to_tag(tag_name) uint8(vr) ...
            typecast(uint16(length(data)+mod(length(data), 2)), 'uint8')]);
          mem.add(data);
          if mod(length(data), 2) ~= 0
            mem.add(uint8(0));
          end
        case {'SH', 'TM', 'CS', 'PN', 'LO', 'IS'}
          mem.add([segdicomfile.name_to_tag(tag_name) uint8(vr) ...
            typecast(uint16(length(data)+mod(length(data), 2)), 'uint8')]);
          mem.add(data);
          if mod(length(data), 2) ~= 0
            mem.add(uint8(20));
          end
        case {'DA', 'UL', 'US'}
          mem.add([segdicomfile.name_to_tag(tag_name) uint8(vr) ...
            typecast(uint16(length(data)), 'uint8')]);
          mem.add(data);
        otherwise
          error('SEGMENT:PANIC', 'Unknown VR, got %s', vr);
      end
    end
    
    function tags = get_tags()
      % Returns a tags struct with tag names as fieldnames and tag as value

      % Create the tags struct
      tags = [];
      tags.FileMetaInfoGroupLength = uint8([2 0 0 0]);
      tags.FileMetaInfoVersion = uint8([2 0 1 0]);
      tags.MediaStorageSOPClassUID = uint8([2 0 2 0]);
      tags.MediaStorageSOPInstanceUID = uint8([2 0 3 0]);
      tags.TransferSyntaxUID = uint8([2 0 16 0]);
      tags.ImplementationClassUID = uint8([2 0 18 0]);
      tags.ImplementationVersionName = uint8([2 0 19 0 ]);
      tags.SOPClassUID = uint8([8 0 22 0]);
      tags.SOPInstanceUID = uint8([8 0 24 0]);
      tags.StudyDate = uint8([8 0 32 0]);
      tags.ContentDate = uint8([8 0 35 0]);
      tags.StudyTime = uint8([8 0 48 0]);
      tags.ContentTime = uint8([8 0 51 0]);
      tags.AccessionNumber = uint8([8 0 80 0]);
      tags.Modality = uint8([8 0 96 0]);
      tags.ConversionType = uint8([8 0 100 0]);
      tags.ReferringPhysicianName = uint8([8 0 144 0]);
      tags.PatientName = uint8([16 0 16 0]);
      tags.PatientID = uint8([16 0 32 0]);
      tags.PatientBirthDate = uint8([16 0 48 0]);
      tags.PatientSex = uint8([16 0 64 0]);
      tags.SecondaryCaptureDeviceManufacturer = uint8([24 0 22 16]);
      tags.SecondaryCaptureDeviceModelName = uint8([24 0 24 16]);
      tags.StudyInstanceUID = uint8([32 0 13 0]);
      tags.SeriesInstanceUID = uint8([32 0 14 0]);
      tags.StudyID = uint8([32 0 16 0]);
      tags.SeriesNumber = uint8([32 0 17 0]);
      tags.InstanceNumber = uint8([32 0 19 0]);
      tags.PatientOrientation = uint8([32 0 32 0]);
      tags.SamplesPerPixel = uint8([40 0 2 0]);
      tags.PhotometricInterpretation = uint8([40 0 4 0]);
      tags.Rows = uint8([40 0 16 0]);
      tags.Columns = uint8([40 0 17 0]);
      tags.BitsAllocated = uint8([40 0 0 1]);
      tags.BitsStored = uint8([40 0 1 1]);
      tags.HighBit = uint8([40 0 2 1]);
      tags.PixelRepresentation = uint8([40 0 3 1]);
      tags.SmallestImagePixelValue = uint8([40 0 6 1]);
      tags.LargestImagePixelValue = uint8([40 0 7 1]);
      tags.PixelData = uint8([224 127 16 0]);
      tags.SegmentData = uint8([251 55 71 25]);
      tags.StartOfItem = uint8([254 255 0 224]);
      tags.EndOfItem = uint8([254 255 13 224]);
      tags.EndOfSequence = uint8([254 255 221 224]);
      tags.Manufacturer = uint8([8 0 112 0]);
      tags.ImageOrientation = uint8([32 0 55 0]);
      tags.ImagePosition = uint8([32 0 50 0]);
    end
    

    function tag = name_to_tag(name)
      % Convert a tag name to a tag.

      tags = segdicomfile.get_tags();
      tag = tags.(name);
    end
    
    function r = serialize_logical(data)
      % Serialize a logical array, doesn't mind shape
      
      r = reshape(uint8(data), 1, numel(data));
    end

    function data = unserialize_logical(r)
      % Unserialize a logical array, doesn't mind shape
      
      data = logical(r);
    end

    function r = serialize_char(data)
      % Serialize a char array, doesn't mind shape

      r = typecast(uint16(data), 'uint8');
    end

    function data = unserialize_char(r)
      % Unserialize a char array, doesn't mind shape

      data = char(typecast(r, 'uint16'));
    end

    function r = serialize_cell(data)
      % Serialize a cell array, doesn't mind shape

      tt = {};
      for n=1:numel(data)
        tt{n} = segdicomfile.serialize(data{n});
      end
      r = segdicomfile.create_chunk(tt);
    end

    function data = unserialize_cell(r)
      % Unserialize a cell array, doesn't mind shape

      tt = segdicomfile.parse_chunk(r);

      data = {};
      for n=1:numel(tt)
        data{n} = segdicomfile.unserialize(tt{n});
      end
    end

    function r = serialize_struct(data)
      % Serialize a struct array, doesn't mind shape

      % Get fieldnames
      fnames = fieldnames(data);

      % Serialize content
      content = {};
      for i=1:numel(data)
        element = cell([0]);
        for j=1:length(fnames)
          element{end+1} = segdicomfile.serialize( data(i).(fnames{j}) );
        end
        content{i} = segdicomfile.create_chunk(element);
      end

      % Serialize fieldnames and pack it all in r
      fnames_ser = {};
      for n=1:length(fnames)
        fnames_ser{n} = segdicomfile.serialize_char(fnames{n});
      end
      r = segdicomfile.create_chunk_va(...
        segdicomfile.create_chunk(fnames_ser), ...
        segdicomfile.create_chunk(content));
    end

    function data = unserialize_struct(r)
      % Unserialize struct array, doesn't mind shape

      % Read fnames and content
      [fnames_raw, content] = segdicomfile.parse_chunk_va(r);
      fnames_raw = segdicomfile.parse_chunk(fnames_raw);
      content = segdicomfile.parse_chunk(content);

      % Parse fnames_raw
      fnames = {};
      for n=1:length(fnames_raw)
        fnames{n} = segdicomfile.unserialize_char(fnames_raw{n});
      end

      % generate empty struct
      data = struct();
      for i=1:numel(fnames)
        data.(fnames{i}) = [];
      end
      data(1) = [];

      % parse content
      for i=1:numel(content)
        field_content = segdicomfile.parse_chunk(content{i});
        for j=1:length(field_content)
          data(i).(fnames{j}) = segdicomfile.unserialize(field_content{j});
        end
      end
    end

    function mem = create_chunk(indata)
      % Serializes a 1xn cell array of uint8 1xn arrays.
      % This is a basic serialization used in alot of places.

      mem = memorybasket();
      mem.add(typecast(uint32(length(indata)), 'uint8'));
      for n=1:length(indata)
        mem.add(typecast(uint32(length(indata{n})), 'uint8'));
        mem.add(indata{n});
      end
    end

    function str = create_chunk_va(varargin)
      % varargin shortcut to create_chunk

      str = segdicomfile.create_chunk(varargin);
    end

    function outdata = parse_chunk(data)
      % Inverse of create_chunk

      pos = 1;
      outdata = {};

      % get size
      s = typecast(data(pos:pos+3), 'uint32');
      pos = pos+4;

      % Read the data
      for n=1:s
        % read len
        len = typecast(data(pos:pos+3), 'uint32');
        pos = pos+4;

        % read data
        outdata{n} = data(pos:pos+len-1);
        pos = pos + len;
      end
    end

    function varargout = parse_chunk_va(data)
      % varargout shortcut to parse_chunk

      varargout = segdicomfile.parse_chunk(data);
      if length(varargout) ~= nargout
        error('SEGMENT:ERROR', ...
          'Invalid serialization stream: Wrong number of parts in chunk');
      end
    end
    
  end
  
  methods(Static = true, Access = public)
   function r = generate_uid()
      % Returns a new random UID (uses the matlab root UID).
      
      % This is the UID root matlab uses
      r = uint8('1.3.6.1.4.1.9590.100.1.2.');
      nums = uint8('1234567890');
      s = RandStream.create('mt19937ar','seed','shuffle');
      r = [r nums( min([floor(s.rand(1)*10) 9]) + 1 )];
      while length(r) < 64
        r = [r nums( min([floor(s.rand(1)*10) 9]) + 1 )];
      end
    end
    
    function stri = secondtostring(t)
      %Converts from seconds to a timestring with hhmmss.sss
      hour = floor(t/3600);
      min = floor((t-hour*3600)/60);
      sec = t-hour*3600-min*60;
      stri = sprintf('%02d%02d%09.6f',hour,min,sec);
    end;    
        
  end
end
