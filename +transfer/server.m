classdef server < handle
  properties(SetAccess = private, GetAccess = private)
    user
    password
    url
  end
  
  methods(Access = public)
    function self = server(user, password, url)
      self.user = user;
      self.password = password;
      self.url = url;
    end
    
    function id = new_study(self, uid, name)
      r = self.send('new_study', {'uid', uid, 'name', name});
      if(~isequal(r{1}, 'ok'))
        error('SEGMENT:ERROR', ['Couldn''t create study: ' r{1}]);
      end
      id = r{2};
    end
    
    function id = new_serie(self, desc, study)
      r = self.send('new_serie', {'desc', desc, 'study', study});
      if(~isequal(r{1}, 'ok'))
        error('SEGMENT:ERROR', ['Couldn''t create serie: ' r{1}]);
      end
      id = r{2};
    end
    
    function fileinfos = new_dicomfiles(self, fileinfos, serie)
      r = self.send('new_dicomfile', {...
          'hash', self.join( {fileinfos.hash} ), ...
          'serie', serie});
      if(~isequal(r{1}, 'ok'))
        error('SEGMENT:ERROR', ['Couldn''t create file: ' r{1}]);
      end
      ids = self.split(r{2});
      for i=1:numel(ids)
        fileinfos(i).id = ids{i};
      end
    end
    
    function update_dicomfile(self, data)
      r = self.send('update_dicomfile', {...
        'file', self.join( {data.fileid} ), ...
        'pos', self.join( arrayfun(@(x) {sprintf('%d', x{1})}, {data.pos}) ), ... 
        'data', self.join( {data.data} ) ...
        } ) ;
      if(~isequal(r{1}, 'ok'))
        error('SEGMENT:ERROR', ['Couldn''t upload part of file: ' r{1}]);
      end
    end
    
    function finalize_dicomfile(self, ids)
      r = self.send('finalize_dicomfile', {'file', self.join(ids)});
      if(~isequal(r{1}, 'ok'))
        error('SEGMENT:ERROR', ['Couldn''t finalize dicomfile: ' r{1}]);
      end
    end

    function id = new_thumbnail(self, thumb_big_hash, thumb_small_hash, thumb_anim_hash, serie, prio)
      r = self.send('new_thumbnail', { ...
        'serie', serie, ...
        'hash_big', self.join(thumb_big_hash), ...
        'hash_small', self.join(thumb_small_hash), ...
        'hash_anim', self.join(thumb_anim_hash), ...
        'prio', self.join(prio)});
      if(~isequal(r{1}, 'ok'))
        error('SEGMENT:ERROR', ['Couldn''t create thumbnail: ' r{1}]);
      end
      id = self.split(r{2});
    end
    
    function update_thumbnail(self, data, type)
      r = self.send('update_thumbnail', {...
        'thumb', self.join( {data.fileid} ), ...
        'type', type, ...
        'pos', self.join( arrayfun(@(x) {sprintf('%d', x{1})}, {data.pos}) )...
        'data', self.join( {data.data} )});
      if(~isequal(r{1}, 'ok'))
        error('SEGMENT:ERROR', ['Couldn''t upload part of thumbnail: ' r{1}]);
      end
    end
    
    function finalize_thumbnail(self, id)
      r = self.send('finalize_thumbnail', {'thumb', self.join(id)});
      if(~isequal(r{1}, 'ok'))
        error('SEGMENT:ERROR', ['Couldn''t finalize thumbnail: ' r{1}]);
      end
    end
    
    function set_thumbnail(self, serie, thumb)
      r = self.send('set_thumbnail', ...
        {'serie', serie, 'thumb', thumb});
      if(~isequal(r{1}, 'ok'))
        error('SEGMENT:ERROR', ['Couldn''t set thumbnail: ' r{1}]);
      end
    end
    
    function finalize_study(self, study)
      r = self.send('finalize_study', ...
        {'study', study});
      if(~isequal(r{1}, 'ok'))
        error('SEGMENT:ERROR', ['Couldn''t finalize study: ' r{1}]);
      end
    end
    
    function study = get_waiting_for_mat(self)
      r = self.send('get_waiting_for_mat', {});
      if isequal(r{1}, 'ok') && length(r)>1
        study = r{2};
      elseif isequal(r{1}, 'nothing')
        study = [];
      elseif isequal(r{1},'noAccess')
        error('SEGMENT:ERROR','NoAccess returned. Probably incorrect credentials.');
      else
        error('SEGMENT:ERROR','An unknown error occurd while trying to fetch series waiting for mat.');
      end
    end
    
    function [pattern] = get_pattern_for_user(self)
      r = self.send('get_pattern_for_user', {});
      if isequal(r{1}, 'ok') && length(r)>1
        pattern = r(2); %return all studies
        return;
      end;
      
      if isequal(r{1},'patternNotSet')
        error('SEGMENT:ERROR','PatternNotSet returned. Pattern is not set for this user');
      end;
      
      if isequal(r{1},'noAccess')
        error('SEGMENT:ERROR','NoAccess returned. Probably incorrect credentials.');
      end;
      
      %Should not get here, should return above in ok clause. Throw error.
      error('SEGMENT:ERROR','An unknown error occurd while trying to get pattern.');

    end;
    
    function [study,names] = get_studies_by_status(self,status)
      r = self.send('get_studies_by_status', {'status', status});
      if isequal(r{1}, 'ok') && length(r)>1
        study = r(4:2:end); %return all studies
        names = r(3:2:end); %return all names
      elseif isequal(r{1}, 'nothing')
        study = [];
        names = {};
      elseif isequal(r{1},'noAccess')
        error('SEGMENT:ERROR','NoAccess returned. Probably incorrect credentials.');
      else
        error('SEGMENT:ERROR',sprintf('An unknown error occurd while trying to fetch series %s.',status)); %#ok<SPERR>
      end
    end;
    
    function [study,names] = get_studies_by_status_and_scantime(self, status, scantime)
      r = self.send('get_studies_by_status_and_scantime', {'status', status, 'scantime', scantime});
      if isequal(r{1}, 'ok') && length(r)>1
        study = r(4:2:end); %return all studies
        names = r(3:2:end); %return all names
      elseif isequal(r{1}, 'nothing')
        study = [];
        names = {};
      elseif isequal(r{1},'noAccess')
        error('SEGMENT:ERROR','NoAccess returned. Probably incorrect credentials.');
      else
        error('SEGMENT:ERROR',sprintf('An unknown error occurd while trying to fetch series %s.',status)); %#ok<SPERR>
      end
    end;
        
    function [study,names] = get_waiting_for_consensus(self)
      r = self.send('get_waiting_for_consensus', {});
      if isequal(r{1}, 'ok') && length(r)>1
        study = r(4:2:end); %return all studies
        names = r(3:2:end); %return all names
      elseif isequal(r{1}, 'nothing')
        study = [];
        names = {};
      elseif isequal(r{1},'noAccess')
        error('SEGMENT:ERROR','NoAccess returned. Probably incorrect credentials.');
      else
        error('SEGMENT:ERROR','An unknown error occurd while trying to fetch series waiting for consensus.');
      end
    end    

    function r = get_intresting_series(self, study)
      r = self.send('get_intresting_series', {'study', study});
      if not(isequal(r{1}, 'ok'))
        error('SEGMENT:ERROR', ['Couldn''t get intresting series: ' r{1}]);
      end
      r(1:2) = [];
    end
    
    function r = get_dicom_files(self, series)
      %Probably obsoleted??? EH:
      r = self.send('get_dicom_files', {'series', series});
      if not(isequal(r{1}, 'ok'))
        error('SEGMENT:ERROR', ['Couldn''t get dicom files: ' r{1}]);
      end
      r(1:2) = [];
    end
    
    function r = get_mat_files(self, study, owner)
      r = self.send('get_mat_files', {'study', study, 'owner', owner});
      if not(isequal(r{1}, 'ok'))
        error('SEGMENT:ERROR', ['Couldn''t get mat files: ' r{1}]);
      end
      r(1:2) = [];
    end            
    
    function r = get_dicom_data(self, dicom)
      r = self.rawsend('get_dicom_data', {'dicom', dicom});
      if isempty(r)
        error('SEGMENT:ERROR', 'Couldn''t get dicom data');
      end
    end
    
    function id = new_matfile(self, study, hash)
      r = self.send('new_matfile', {'study', study, 'hash', hash});
      if not(isequal(r{1}, 'ok'))
        error('SEGMENT:ERROR', ['Couldn''t create new matfile: ' r{1}]);
      end
      id = r{2};
    end
    
    function update_matfile(self, id, pos, data)
      r = self.send('update_matfile', {'matfile', id, 'pos', sprintf('%d', pos), 'data', data});
      if(~isequal(r{1}, 'ok'))
        error('SEGMENT:ERROR', ['Couldn''t upload part of matfile: ' r{1}]);
      end
    end
    
    function update_status(self,id,newstatus)
      r = self.send('update_status', ...
        {'study', id, 'newstatus', newstatus});
      if not(isequal(r{1}, 'ok'))
        error('SEGMENT:ERROR', ['Couldn''t update status: ' r{1}]);
      end
    end;
    
    function finalize_matfile(self, matfile)
      r = self.send('finalize_matfile', ...
        {'matfile', matfile});
      if(~isequal(r{1}, 'ok'))
        error('SEGMENT:ERROR', ['Couldn''t finalize matfile: ' r{1}]);
      end
    end

  end
  
  methods(Access = private)
    function r = send(self, action, args)
      r = self.unpack_result(self.rawsend(action, args));
      %r = self.unpack_result(urlread(self.url, 'post', ...
      %  [{'user', self.user, 'password', self.password, 'action', action} args]));
    end
  end

  methods(Access = public)
    function r = rawsend(self, action, args)
      failed_attempts = 0;
      while true
        try
          r = urlread(self.url, 'post', ...
               [{'user', self.user, 'password', self.password, 'action', action} args]);
        catch e
          if strcmp(e.message, 'Error downloading URL.')
            if failed_attempts > 3
              error('SEGMENT:ERROR', 'Couldn''t connect to core lab server');
            else
              failed_attempts = failed_attempts + 1;
              continue;
            end
          else
            rethrow(e);
          end
        end
        break
      end
      
    end
  end

  methods(Static, Access = public)
    function r = unpack_result(res)
      r = {};
      while(numel(res) ~= 0)
        [A, ~, ~, nextindex] = sscanf(res, '%s', 1);
        r{end+1} = A;
        res = res(nextindex:end);
      end
    end
    
    function r = join(str_list)
      r = '';
      for i=1:numel(str_list)
        if i == 1
          r = str_list{1};
        else
          r = [r '|' str_list{i}];
        end
      end
    end
    
    function r = split(data)
      r = regexp(data, '\|', 'split');
    end
  end
end