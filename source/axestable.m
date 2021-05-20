classdef axestable < handle % inherrit from handle to have persistant objects
  % axestable is an object that can attatch to a axes object in matlab
  % and draw table in it. It was writting for segment CVQ component and is
  % therefor limited in apearence.
  %
  % There are 3 public properties that control the look of the table:
  %    -backgroundcolor
  %    -fontcolor
  %    -defaultclearstring
  %
  % Init a axestable by attaching it to an axes object.
  %     table=axestable(axes_handle);
  %
  % One axestable can draw several tables in a axes object.
  % Each table can then hold several keys. They are fields with
  % accompanating value field(s). Each key holds:
  %   Key, Name, unit, values.
  %
  % Add a table by calling addTable method
  % table.addTable('Title',total_columns,nbr_values [,colsize]);
  %
  %     Total_columns: How many total colums the table conist of
  %     nbr_values: How many value fields each key in table hold.
  %                 mod(total_columns,nbr_value)=0 needs to be furfilled
  %                 for a table. Else it wont be able to fit keys evenly on
  %                 each row in table.
  %
  % Add a key to the latest created table:
  % table.addKey(key,'Name','unit',{'values'})
  %
  % Once the table is created it must be drawn:
  % table.draw();
  %
  % When the table is created there are several update methods to update
  % values in keys:
  %
  %     method                             function
  %     t.clearAllValues(clearStr)    Set all values to a clearStr.
  %
  %     t.clearByUnit(unit,clearStr)  Set all values to clearStr where
  %                                   Unit==key.Unit.
  %
  %     t.updateKey(key,value)        Set key.values=values.
  %
  %     t.setVisible(key,visible)     Toggle the visibility of a key.
  %
  % Most methods take a optional argument draw=true for instant update of
  % graphics. draw is by default false. This is to change several methods
  % and then do a graphics update by caling table.draw();
  % The visibilty of the complete axes including tables can be toggled by
  % the table.hide() table.show() methods.
  %
  % This code is copyrighted and is the property of Medviso all rights are
  % resereved.
  % Created by Erik Södervall
  %
  % Einar made some changes to better reflect coding standard
  % and include fontsize.
  %
  % Helen Fransson made some changes 2017 to extend the usage of the table
  % class to be used for Segment CMR, Segment CT and Segment

  properties
    backgroundcolor='black';
    fontcolor='white';
    fontsize=8;
    defaultclearstring='---';
    ystep=20;
  end
  properties (SetAccess = 'private',Hidden)
    axes=[];
    xstart=10;
    xpos=10;
    ypos=0;
    FigWith=0;
    FigHeight=0;
    NbrTables=0;
    Key={};
    Table={};
    createHandles=true;
    hidden=true;
    linex=[];
    liney=[];
  end
  methods
    
    %---------------------------------
    function obj=axestable(axesHandle)
    %--------------------------------
    
      obj.axes=axesHandle;
      obj.xstart=10;
      obj.xpos=obj.xstart;
      obj.ypos=0;
      units=get(axesHandle,'Units');
      set(axesHandle,'Units','Pixels');
      pos=get(axesHandle,'Position');
      vis=get(axesHandle,'Visible');
      if strcmp(vis,'off')
        obj.hidden=true;
      else
        obj.hidden=false;
      end
      set(axesHandle,'Units',units);

      obj.FigWith=300;%pos(3);
      obj.FigHeight=pos(4);
      obj.NbrTables=0;
      obj.Key={};
      obj.Table={};
    end

    %-----------------------------------------
    function clearAllValues(obj,clearStr,draw)
    %-----------------------------------------
      if nargin<2
        clearStr='---';
      end
      if nargin<3
        draw=false;
      end
      %set all values to clearStr in table.
      for loop=1:length(obj.Key) %Update all values.
        for vloop=1:length(obj.Key{loop}.values)
          obj.Key{loop}.values{vloop}=clearStr;
        end
        if draw&&not(obj.createHandles)
          for vloop=1:length(obj.Key{loop}.values)
            set(obj.Key{loop}.valuesHandle{vloop},'String',clearStr);
          end
        end
      end
    end

    %-------------------------------------------
    function clearByUnit(obj,unit,clearStr,draw)
    %-------------------------------------------
      if nargin<3
        clearStr='---';
      end

      if nargin<4&&islogical(clearStr)
        draw=true;
        clearStr='---';
      elseif nargin<4
        draw=false;
      end

      for loop=1:length(obj.Key) %Update all values.
        if strcmp(obj.Key{loop}.Unit,unit)
          for vloop=1:length(obj.Key{loop}.values)
            obj.Key{loop}.values{vloop}=clearStr;
          end
          if draw&&not(obj.createHandles)
            for vloop=1:length(obj.Key{loop}.values)
              set(obj.Key{loop}.valuesHandle{vloop},'String',clearStr);
            end
          end
        end
      end
    end

    %-----------------
    function show(obj)
    %-----------------
      if not(obj.hidden)
        return
      end
      obj.hidden=false;
      obj.setVisible('on');
      obj.draw();
      axis(obj.axes,'off')
    end

    %-----------------
    function hide(obj)
    %-----------------
      if obj.hidden
        return
      end

      obj.hidden=true;
      obj.setVisible('off');
    end

    %-----------------
    function draw(obj)
    %-----------------
      if obj.NbrTables==0
        disp('No tables to draw');
        return
      end

      if obj.hidden && not(obj.createHandles)
        return %dont do anything when object is invisible.
      end

      if obj.createHandles %first time. Create all handles
        obj.createHandles=false;
        firstdraw(obj);
        if(obj.hidden)
          obj.setVisible('off');
        else
          obj.setVisible('on');
        end
        return;
      end
      
      set(obj.axes,'Color',obj.backgroundcolor);
      set(obj.linex,'color',obj.backgroundcolor)
      set(obj.liney,'color',obj.backgroundcolor)
      
      for loop=1:length(obj.Table) %Update all titles
        stri = translation.dictionary(obj.Table{loop}.Title);
        set(obj.Table{loop}.Handle,'String',stri,'color',obj.fontcolor)
      end
      for loop=1:length(obj.Key) %Update all values.
        for vloop=1:length(obj.Key{loop}.values)
          stri=obj.parseValue(obj.Key{loop}.values{vloop});
          set(obj.Key{loop}.valuesHandle{vloop},...
            'String',stri,...
            'Position',[obj.Key{loop}.valuesXpos(vloop) obj.FigHeight-obj.Key{loop}.valuesYpos(vloop)],'color',obj.fontcolor);
%           stri=obj.parseValue(obj.Key{loop}.Name);
          if ~isempty(obj.Key{loop}.Unit)
            stri=sprintf('%s (%s)',translation.dictionary(obj.Key{loop}.Name),obj.Key{loop}.Unit);
          else
            stri=translation.dictionary(obj.Key{loop}.Name);
          end
          set(obj.Key{loop}.keyHandle,'String',stri,'color',obj.fontcolor);
        end
      end
    end

    %-------------------------------------------------
    function addTable(obj,title,col,nbrValues,colsize)
    %-------------------------------------------------
      if isempty(obj.axes)
        error('No handle, call setupfigtable');
      end
      if ~obj.createHandles
        error('Forbidden to add Tables once table is drawn')
      end

      if nargin<3
        col=4;
        nbrValues=1;
      end

      if mod(col,nbrValues+1)~=0
        error('Need to be able to fit values evenly per row');
      end

      obj.NbrTables=obj.NbrTables+1;
      obj.Table{obj.NbrTables}=[];
      obj.Table{obj.NbrTables}.Cols=col;
      obj.Table{obj.NbrTables}.CurCols=1;
      obj.Table{obj.NbrTables}.nbrValues=nbrValues;
      obj.Table{obj.NbrTables}.ystep=obj.ystep;

      if nargin<5
        obj.Table{obj.NbrTables}.xstep=ones(1,col)*obj.FigWith/col;
      elseif sum(colsize)<1+1e-5&&sum(colsize)>1-1e-5 &&...
          length(colsize)==col
        obj.Table{obj.NbrTables}.xstep=colsize*(obj.FigWith-obj.xstart);
      else
        obj.Table{obj.NbrTables}.xstep=ones(1,col)*obj.FigWith/col;
        disp('Ignoring colsize, does not sum to 1');
      end

      if obj.NbrTables~=1
        %add some space before title if more than one table in axes
        obj.ypos=obj.ypos+obj.Table{obj.NbrTables}.ystep;
%         obj.Table{obj.NbrTables-1}.ypos = obj.Table{obj.NbrTables-1}.ypos+obj.Table{obj.NbrTables}.ystep;
      end
      if strcmp(title,'LV') && obj.ypos == 0 %this is for fixing the incorrect placement of the LV table in Segment CMR
        obj.ypos = 20;
      end

      obj.Table{obj.NbrTables}.xpos=obj.xstart;
      obj.Table{obj.NbrTables}.ypos=obj.ypos;
      obj.Table{obj.NbrTables}.Title=title;
      obj.ypos=obj.ypos+obj.Table{obj.NbrTables}.ystep+5;%extra space for line
      obj.xpos=obj.xstart;
    end

    %---------------------------------------
    function addTableHead(obj,column_names)
    %--------------------------------------
    %add the head of the table
      if obj.NbrTables==0
        disp('addTableHeader: Table must exist');
        return
      end

      if obj.Table{obj.NbrTables}.nbrValues+1~=length(column_names)
        disp('addTableHeader: Wrong number or fields.');
        return;
      end

      if obj.Table{obj.NbrTables}.CurCols~=1
        disp('addTableHeader: Call this method directly after addTable');
        return;
      end

      obj.Table{obj.NbrTables}.Head.Title=column_names;
      rep=obj.Table{obj.NbrTables}.Cols/length(column_names);

      obj.Table{obj.NbrTables}.Head.Repeat=rep;
      obj.Table{obj.NbrTables}.Head.ypos=obj.ypos;
      obj.Table{obj.NbrTables}.Head.lineypos=obj.ypos+5;
      obj.ypos=obj.ypos+obj.Table{obj.NbrTables}.ystep+5;
      obj.Table{obj.NbrTables}.Head.TitleHandle={};
    end

    %----------------------------------------
    function addKey(obj,key,name,unit,values)
    %----------------------------------------
    %add a key (line) to the table

      if isempty(obj.axes)
        error('No handle, call setupfigtable');
      end

      if ~obj.createHandles %TODO make it handle proper redraw.
        error('Forbidden to add Keys once table is drawn')
      end

      if obj.NbrTables==0
        error('No table, call newtable');
      end

      if ~isa(values,'cell')
        values={values};
      end

      if length(values)~=obj.Table{obj.NbrTables}.nbrValues
        error('Wrong number of values for current table');
      end

      %TODO: ensure no dublicate keys
      nbr_keys=length(obj.Key);
      nbr_keys=nbr_keys+1;%new key:

      obj.Key{nbr_keys}.Key=key;
      obj.Key{nbr_keys}.Name=name;
      obj.Key{nbr_keys}.Unit=unit;
      obj.Key{nbr_keys}.showKey=true;
      obj.Key{nbr_keys}.xpos=obj.xpos;
      obj.Key{nbr_keys}.ypos=obj.ypos;

      obj.xpos=obj.xpos+obj.Table{obj.NbrTables}.xstep(obj.Table{obj.NbrTables}.CurCols);

      obj.Table{obj.NbrTables}.CurCols=obj.Table{obj.NbrTables}.CurCols+1;
      obj.Key{nbr_keys}.valuesHandle=[];
      obj.Key{nbr_keys}.values=values;
      obj.Key{nbr_keys}.valuesXpos=[];
      obj.Key{nbr_keys}.valuesYpos=[];

      for loop=1:length(values) %figure out postition
        obj.Key{nbr_keys}.valuesXpos(loop)=obj.xpos;
        obj.Key{nbr_keys}.valuesYpos(loop)=obj.ypos;
        obj.xpos=obj.xpos+obj.Table{obj.NbrTables}.xstep(obj.Table{obj.NbrTables}.CurCols);
        obj.Table{obj.NbrTables}.CurCols=obj.Table{obj.NbrTables}.CurCols+1;
      end

      if obj.Table{obj.NbrTables}.CurCols>...
          obj.Table{obj.NbrTables}.Cols
        obj.ypos=obj.ypos+obj.Table{obj.NbrTables}.ystep;
        obj.Table{obj.NbrTables}.CurCols=1;
        obj.xpos=obj.xstart;
      end
    end

    %---------------------
    function addSpace(obj)
    %---------------------
      if obj.NbrTables==0
        obj.ypos=obj.ypos+obj.ystep;
      else
        obj.xpos=obj.xstart;
        obj.ypos=obj.ypos+obj.Table{obj.NbrTables}.ystep;
        obj.Table{obj.NbrTables}.CurCols=1;
      end
    end

    %--------------------------------------
    function updateKey(obj,key,values,draw)
    %--------------------------------------
      %get key update value then return..
      %also redraw when draw==true
      %key can be cell of keys for multiple updates
      if nargin<4
        draw=false;
      end
      if ~isa(values,'cell')
        values={values};
      end

      for loop=1:length(obj.Key)
        if strcmp(obj.Key{loop}.Key,key)
          if length(values)~=length(obj.Key{loop}.values)
            error('Nbr of values does not match');
          end
          obj.Key{loop}.values=values;

          if draw&&~obj.createHandles
            for vloop=1:length(values)
              stri=obj.parseValue(obj.Key{loop}.values{vloop});
              set(obj.Key{loop}.valuesHandle{vloop},'String',stri);
            end
          end

          if ~isa(key,'cell')
            return
          end
        end
      end

      if ~isa(key,'cell')
        error('Could not find key in table');
        return
      end
    end
    
    %--------------------------------------
    function updateName(obj,key,names,draw)
    %--------------------------------------
      %update names then return
      %also redraw when draw==true
      %key can be cell of keys for multiple updates
      if nargin<4
        draw=false;
      end
      if isa(names,'cell')
        names=names{1};
      end

      for loop=1:length(obj.Key)
        if strcmp(obj.Key{loop}.Key,key)
%           if length(names)~=length(obj.Key{loop}.Name),
%             error('Nbr of names does not match');
%           end
          obj.Key{loop}.Name=names;
        end
      end
    end
    
    %--------------------------------------
    function updateUnit(obj,key,unit,draw)
    %--------------------------------------
      %update names then return
      %also redraw when draw==true
      %key can be cell of keys for multiple updates
      if nargin<4
        draw=false;
      end
      if isa(unit,'cell')
        unit=unit{1};
      end

      for loop=1:length(obj.Key)
        if strcmp(obj.Key{loop}.Key,key)
%           if length(names)~=length(obj.Key{loop}.Name),
%             error('Nbr of names does not match');
%           end
          obj.Key{loop}.Unit=unit;
        end
      end
    end
    
    %--------------------------------------
    function updateTitle(obj,title,draw)
    %--------------------------------------
      %update names then return
      %also redraw when draw==true
      %key can be cell of keys for multiple updates
      if nargin<4
        draw=false;
      end

      if length(obj.Table) ~= length(title)
        error('Nbr of titles does not match');
      end
      for loop=1:length(obj.Table)
        obj.Table{loop}.Title = title{loop};
      end
    end
    
    %----------------------------------------
    function keyVisible(obj,key,visible,draw)
    %----------------------------------------
      if nargin<4
        draw=false;
      end
      for loop=1:length(obj.Key)
        if strcmp(obj.Key{loop}.Key,key)
          obj.Key{loop}.showKey=visible;

          if visible
            vis='on';
          else
            vis='off';
          end
          if draw&&not(obj.createHandles)
            set(obj.Key{loop}.keyHandle,'Visible',vis);
            for vloop=1:length(obj.Key{loop}.values)
              set(obj.Key{loop}.valuesHandle{vloop},'Visible',vis);
            end
          end
          return;
        end
      end
    end

  end  %end of public methods
    
  methods (Access = 'private')
    
    %---------------------------------
    function h=putstring(obj,stri,x,y)
    %---------------------------------
      h=text('Position',[x obj.FigHeight-y],...
        'String',stri,...
        'VerticalAlignment','Bottom',...
        'FontSize',obj.fontsize,...
        'Color',obj.fontcolor,...
        'parent',obj.axes);
    end

    %------------------------
    function h=putline(obj,y) 
    %------------------------  
      %always will want to put the same line at y
      h=line([obj.FigWith 0],[obj.FigHeight-y obj.FigHeight-y],'Color',obj.fontcolor,'parent',obj.axes);
    end

    %----------------------
    function firstdraw(obj)
    %----------------------
      % draw lines in background to create arena
      set(obj.axes,'XTick', [], 'YTick', []);
      set(obj.axes,'Color',obj.backgroundcolor);
      axis(obj.axes,'off')

      obj.linex=line([0 obj.FigWith],[0 0],'Color',obj.backgroundcolor,'parent',obj.axes);
      obj.liney=line([0 0],[obj.FigHeight 0],'Color',obj.backgroundcolor,'parent',obj.axes);
%       obj.linex=line([0 0],[0 0],'Color',obj.backgroundcolor,'parent',obj.axes);
%       obj.liney=line([0 0],[0 0],'Color',obj.backgroundcolor,'parent',obj.axes);


      for loop=1:obj.NbrTables
        obj.Table{loop}.Handle = putstring(obj,...
          obj.Table{loop}.Title,obj.Table{loop}.xpos,obj.Table{loop}.ypos);
        obj.Table{loop}.LineHandle = putline(obj,obj.Table{loop}.ypos);

        if isfield(obj.Table{loop},'Head')
          obj.Table{loop}.Head.LineHandle=...
            putline(obj,obj.Table{loop}.Head.lineypos);
          col_lenth=length(obj.Table{loop}.Head.Title);
          x_pos=[obj.xstart obj.xstart+cumsum(obj.Table{obj.NbrTables}.xstep(1:end-1))];
          for vloop=1:col_lenth*obj.Table{loop}.Head.Repeat
            ind=mod(vloop-1,col_lenth)+1;
            obj.Table{loop}.Head.TitleHandle{vloop}=putstring(obj,...
              obj.Table{loop}.Head.Title{ind}, x_pos(vloop) ,...
              obj.Table{loop}.Head.ypos);
          end
        end
      end

      for loop=1:length(obj.Key)
        if ~isempty(obj.Key{loop}.Unit)
          stri=sprintf('%s (%s)',obj.Key{loop}.Name,obj.Key{loop}.Unit);
        else
          stri=obj.Key{loop}.Name;
        end

        obj.Key{loop}.keyHandle= putstring(obj,stri, ... %create key
          obj.Key{loop}.xpos,obj.Key{loop}.ypos);
        for vloop=1:length(obj.Key{loop}.values) %create values
          obj.Key{loop}.valuesHandle{vloop}= putstring(obj,obj.Key{loop}.values{vloop},...
            obj.Key{loop}.valuesXpos(vloop),...
            obj.Key{loop}.valuesYpos(vloop));
        end
      end

    end

    %---------------------------------
    function str=parseValue(obj,value)
    %---------------------------------
      switch class(value)
        case {'double','single'}
          if isnan(value)
            str=obj.defaultclearstring;
          elseif value==round(value)
            str=sprintf('%3d',value);
          elseif value < 10 && value > -10
            str=sprintf('%3.2f',value);
          else
            str=sprintf('%3.1f',value);
          end
        case 'string'
          str=value;
        case {'int8','uint8','int16','uint16',...
            'int32','uint32','int64','uint64', 'logical'}
          str=sprintf('%d',value);
        case 'char'
          str=sprintf('%c',value);
        otherwise
          str='???';
      end
    end

    %-----------------------------
    function setVisible(obj,value)
    %-----------------------------      
      %TODO: Fewer set calls with more handles in each.
%       set(obj.axes,'Visible',value);
      if obj.createHandles
        return; %other handles does not exist yet
      end

%       set(obj.linex,'Visible',value);
%       set(obj.liney,'Visible',value);

      for loop=1:obj.NbrTables
        set(obj.Table{loop}.Handle,'Visible',value);
        set(obj.Table{loop}.LineHandle,'Visible',value);
        if isfield(obj.Table{loop},'Head')
          set(obj.Table{loop}.Head.LineHandle,'Visible',value);
          for vloop=1:length(obj.Table{loop}.Head.TitleHandle)
            set(obj.Table{loop}.Head.TitleHandle{vloop},'Visible',value)
          end
        end
      end

      for loop=1:length(obj.Key) %Update all values.
        if obj.Key{loop}.showKey
          set(obj.Key{loop}.keyHandle,'Visible',value);
          for vloop=1:length(obj.Key{loop}.values)
            set(obj.Key{loop}.valuesHandle{vloop},'Visible',value);
          end
        else
          set(obj.Key{loop}.keyHandle,'Visible','off');
          for vloop=1:length(obj.Key{loop}.values)
            set(obj.Key{loop}.valuesHandle{vloop},'Visible','off');
          end
        end
      end
    end

  end
end
