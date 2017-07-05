function resulturl = sendtogensvar(url, username, password, data)
    % This function creates a study in gensvar and returns the url of the
    % study. For a list of valid fieldnames in the data structure see
    % the function getStudyFieldGroup() in common.php.
    
    % Prepare data 
    fnames = fieldnames(data);
    sdata = {};
    for i=1:numel(fnames)
        value = data.(fnames{i});
        if isa(value, 'char')
            sdata = [sdata {fnames{i}} {value}];
            continue
        end
        if isa(value, 'double') || isa(value, 'single')
            sdata = [sdata {fnames{i}} {sprintf('%f', value)}];
            continue
        end
        error('SEGMENT:ERROR', 'Bad value type in data parameter');
    end
    
    % Send
    rr = urlread([url '/api.php'], 'POST', [{'action', 'create_study', 'username', username, 'password', password} sdata]);
    r = regexp(rr, '\n', 'split');
    if(not(isequal(r{1}, 'gensvar')))
        error('SEGMENT:ERROR', ['Couldn''t connect to gensvar: ' rr]);
    end
    if(isequal(r{2}, 'ok'))
        resulturl = [url '/editstudy.php?id=' r{3}];
        return
    end
    if(isequal(r{2}, 'noaccess'))
        error('SEGMENT:ERROR', 'Bad combination of username and password');
    end
    if(isequal(r{2}, 'baddata'))
        error('SEGMENT:ERROR', ['Bad data sent. Gensvar says: ''' r{3} '''.']);
    end
    error('SEGMENT:ERROR', 'Unknown result sent from gensvar');
end