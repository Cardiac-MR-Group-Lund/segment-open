classdef presetview
  
  properties
    
    name
    hotkey
    
    nbrofims = [];
    matrix = [];
    panelstype = [];
    
    ImageOrientation %(this will only be used if it contains ones and zeros)
    ImageType
    ImageViewPlane
    SequenceName
    SeriesDescription
    DICOMImageType
    IsEmpty %to handle empty frames
    IsTimeResolved %(this is not stored, but calculated as TSize>1)
    IsSingleSlice %(this not stored but calculated as ZSize==1)
    IsCine %(true if returned as cine in findno call)
    IsScar %(true if returned as scar in findno call)
    IsFlow %(true if returned as flow in findno call)
    IsMar %(true if returned as MaR in findno call)
    IsPerfusion %(true if returned as Perfusion in findno call) Reserved for future use
    IsT2Star %(true if...) Reserved for future use
    %Maybe some more deduced image type 
    
  end
  
  methods
    
    %-----------------------------------
    function v = presetview(nmin,hkeyin)
    %-----------------------------------
    % Constructor
    global DATA SET
    
    if nargin < 2
      hkeyin = '';
    end
    
    v.name = nmin;
    v.hotkey = hkeyin;
    v.matrix = DATA.ViewMatrix;
    v.panelstype = DATA.ViewPanelsType;
    nbrofims = prod(v.matrix); %#ok<*PROP>
    v.nbrofims = nbrofims;
    
    prp = setdiff(properties(v),{'name','hotkey','nbrofims','matrix','panelstype'});
    isprp = prp(strncmp('Is',prp,2));
    fdprp = setdiff(prp,isprp); %properties directly taken from fields of SET struct
    [cineno,scarno,flowno,~,marno] = findfunctions('findno');
    
    for i = 1:numel(isprp)
      v.(isprp{i}) = zeros(1,nbrofims);      
    end
    for i = 1:numel(fdprp)
      v.(fdprp{i}) = cell(1,nbrofims);
    end
    
    for i = 1:nbrofims
      no = DATA.ViewPanels(i);
      if no == 0
        v.IsEmpty(i) = 1;        
      else
        setstr = SET(no);
        for fi = 1:numel(fdprp)
          field = fdprp{fi};
          v.(field){i} = setstr.(field);
        end
        v.IsTimeResolved(i) = setstr.TSize > 1;
        v.IsSingleSlice(i) = setstr.ZSize == 1;
        v.IsCine(i) = ismember(no,cineno);
        v.IsScar(i) = ismember(no,scarno);
        v.IsFlow(i) = ismember(no,flowno);
        v.IsMar(i) = ismember(no,marno);
        v.IsPerfusion(i) = 0; %future use
        v.IsT2Star(i) = 0; %future use
      end
    end
    
    end
    
    %----------------------
    function nos = match(v)
    %----------------------
    % Find nos that best match stored view
    global SET
    nos = zeros(1,v.nbrofims);
    prp = setdiff(properties(v),{'name';'hotkey';'nbrofims';'matrix';'panelstype'});
    isprp = prp(strncmp('Is',prp,2));
    fdprp = setdiff(prp,isprp); %properties directly taken from fields of SET struct
    w3prp = {'ImageType',...
      'ImageViewPlane',...
      'IsTimeResolved',...
      'IsSingleSlice'};
    w1prp = setdiff(prp,w3prp); %Properties that are weighted at one point
    [cineno,scarno,flowno,~,marno] = findfunctions('findno');
    
    for i = 1:v.nbrofims
      if v.IsEmpty(i)
        nos(i) = 0;
      else
        points = zeros(1,numel(SET));
        for no = 1:numel(SET)
          setstr = SET(no);
          
          %Go through field properties to look for matchings
          for j = 1:numel(fdprp)
            field = fdprp{j};
            if isequal(setstr.(field),v.(field){i})
              if ismember(field,w3prp)
                points(no) = points(no)+3;
              else
                points(no) = points(no)+1;
              end
            end
          end
          
          %Now for the 'Is' properties
          isvec = zeros(numel(isprp),1);
          isvec(1) = v.IsTimeResolved(i) == (setstr.TSize > 1);
          isvec(2) = v.IsSingleSlice(i) == (setstr.ZSize == 1);
          isvec(3) = v.IsCine(i) == ismember(no,cineno);
          isvec(4) = v.IsScar(i) == ismember(no,scarno);
          isvec(5) = v.IsFlow(i) == ismember(no,flowno);
          isvec(6) = v.IsMar(i) == ismember(no,marno);
          isvec(isvec == 0) = -1;
          isvec(7) = 0; %future use
          isvec(8) = 0; %future use
          
          weights = ones(numel(isprp),1)+2*ismember(isprp,w3prp);
          points(no) = points(no) + weights'*isvec; %scalar product
        end
        
        [~,nos(i)] = max(points);
      end
    end
    
    end
    
  end
  
end