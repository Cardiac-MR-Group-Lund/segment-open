function names = getperfusionnbrs(fieldkids,field)
      global DATA
      gui = DATA.GUI.Perfusion;
      
      upslopes = cell(16,1);
      smoothupslopes = cell(16,1);
      myored = cell(16,1);
      myoblack = cell(16,1);
      tangents = cell(16,1);
      startline = get(gui.handles.([field 'startline']),'ydata');
      endline = get(gui.handles.([field 'endline']),'ydata');
      
      bloodpool = get(gui.handles.([field 'bloodpool']), 'ydata');
      smoothbloodpool = get(gui.handles.([field 'smoothbloodpool']), 'ydata');
      
      bloodred = get(gui.handles.([field 'bloodpoolred']), 'ydata');
      bloodblack = get(gui.handles.([field 'bloodpoolblack']), 'ydata');
      bptangent = get(gui.handles.([field 'bptangents']), 'ydata');
    
      
      for val = 6
        upslopes{val} = 0;
        smoothupslopes{val} = get(gui.handles.([field 'smoothupslopes'])(val),'ydata');
        myored{val} = get(gui.handles.([field 'myored'])(val),'ydata');
        myoblack{val} = get(gui.handles.([field 'myoblack'])(val),'ydata');
        tangents{val} = get(gui.handles.([field 'tangents'])(val),'ydata');
      end
      
      names = cell(length(fieldkids),1);
      for n = 1:length(fieldkids)
        %Loop through all kids:
        ykidstemp = get(fieldkids(n),'YData');
        foundupslopes = nnz(cellfun(@(x) isequal(x,ykidstemp),upslopes));
        
        if foundupslopes==0
          foundsmoothupslopes = nnz(cellfun(@(x) isequal(x,ykidstemp),smoothupslopes));
          if foundsmoothupslopes==0
            foundmyored = nnz(cellfun(@(x) isequal(x,ykidstemp),myored));
            if foundmyored==0
              foundmyoblack = nnz(cellfun(@(x) isequal(x,ykidstemp),myoblack));
              if foundmyoblack==0
                foundtangents = nnz(cellfun(@(x) isequal(x,ykidstemp),tangents));
                if foundtangents==0
                  foundstart = nnz(isequal(startline,ykidstemp));
                  if foundstart==0
                    foundend = nnz(isequal(endline,ykidstemp));
                    if foundend==0
                       foundblood = nnz(isequal(bloodpool,ykidstemp));
                       if foundblood==0
                         foundsmoothblood = nnz(isequal(smoothbloodpool,ykidstemp));
                         if foundsmoothblood==0
                           foundbloodred = nnz(isequal(bloodred,ykidstemp));
                           if foundbloodred==0
                             foundbloodblack = nnz(isequal(bloodblack,ykidstemp));
                             if foundbloodblack==0
                               foundbptangent = nnz(isequal(bptangent,ykidstemp));
                               if foundbptangent==0
                                 names{n} = [];
                               else
                                 names{n} = 'bptangent';
                               end
                             else
                               names{n} = 'bloodblack';
                             end
                             
                           else
                             names{n} = 'bloodred';
                           end
                           
                         else
                           names{n} = 'smoothbloodpool';
                         end
                       else
                         names{n} = 'bloodpool';
                       end
                    else
                      names{n} = 'endline';
                    end
                  else
                    names{n} = 'startline';
                  end
                else
                  names{n} = 'tangents';
                end
              else
                names{n} = 'myoblack';
              end
            else
              names{n} = 'myored';
            end
          else
            names{n} = 'smoothupslope';
          end
        else
          names{n} = 'upslope';
        end
      end
    end