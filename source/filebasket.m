classdef filebasket < handle
  properties
    fid
    len
  end
  
  methods
    function self = filebasket(filename, mode)
      % append  '\\?\' to avoid MAX_PATH_LENGTH problem
      %first check that it is Windows, long filepath and that the
      %files is local (does not start with \\)
      if iscell(filename)
        filename = char(filename);
      end
       if ispc && length(filename) > 250 && ~strcmp(filename(1:2), '\\')
         filename = append('\\?\',filename);
       end
      self.fid = fopen(filename, mode);
      if self.fid == -1
        errorstr = dprintf('Could not open file.');
        error('SEGMENT:ERROR', errorstr);
      end
      fseek(self.fid, 0, 'eof');
      self.len = ftell(self.fid);
      fseek(self.fid, 0, 'bof');
    end
    
    function add(self, data)
      if strcmp(class(data), 'memorybasket')
        data = data.render();
      end
      fwrite(self.fid, data);
      self.len = self.len + length(data);
    end
    
    function r = pop(self, len)
      len = double(len);
      r = fread(self.fid, len, '*uint8')';
      if len > length(r)
        error('SEGMENT:ERROR', ...
          'Invalid stream: Not Enough bytes left');
      end
    end
    
    function r = bytes_left(self)
      r = self.len - ftell(self.fid);
    end
    
    function close(self)
      fclose(self.fid);
    end
  end
end