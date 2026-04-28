classdef saveoptionfunctions
  % class containing all functions related to save options
  % initial implementation with focus on image/movie

  methods (Static)
    %----------------------
    function outoptions = getdefaultsaveoptions(functionname)
      %----------------------
      % get save options for provided function name

      opt = saveoptionfunctions.getallsaveoptions;
      switch functionname
        case {'flow4d'}
          savekey = {'pacs','avi','image','gif'};
        otherwise
          outoptions = {''};
          return
      end
      outoptions = values(opt,savekey);
    end

    %----------------------
    function opt = getallsaveoptions
      %----------------------
      % key value pairs for all save options
      saveopt = {
        'pacs', dprintf('Save to PACS');... %key - value pairs
        'avi', 'avi';...
        'image', dprintf('Images');...
        'gif',  dprintf('animated gif');...
        };
      opt = containers.Map(saveopt(:,1), saveopt(:,2));
    end

    %------------------------------------------------------
    function optkey = getchosensaveoption(functionname,ind)
      %----------------------------------------------------
      % get default save options as container
      defopt = saveoptionfunctions.getallsaveoptions;
      % Get the current save options for input function name (provided in
      % the order as used in the corresponding function)
      opt = saveoptionfunctions.getdefaultsaveoptions(functionname);
      % Get the default index of the value in 'opt'
      defind = matches(defopt.values,opt{ind});
      % get key name for the found index
      defkeys = defopt.keys;
      optkey = defkeys{defind};
    end

    %-------------------------------------------------
    function settingslist = getsettingslist(mediatype)
      %-------------------------------------------------

      switch mediatype
        case 'pacs'
          pacslist = pacsaccess('getpacslist','send');
          settingslist = {pacslist.DescriptiveName};
        case 'avi'
          settingslist = {
            'None' ...
            'Motion JPEG' ...
            'Motion JPEG 2000' ...
            'Motion JPEG 2000 (lossless)' ...
            };
        case 'image'
          settingslist = {
            '.png',...
            '.jpg',...
            '.bmp',...
            '.tif', ...
            };
        case 'gif'
          settingslist = {
            dprintf('Dithering'),...
            dprintf('No dithering')
            };
        otherwise
          settingslist = {'             '};
      end
    end

  end
end
