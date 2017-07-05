function mesh2stl(x,y,z,resolution,closeapex,fignr,pathname,filename)
%MESH2STL(X,Y,Z,Resolution,CloseApex,Pathname,Filename)
%
%or
%
%MESH2STL(X,Y,Z,Resolution,CloseApex,fid)
%
%Takes a mesh (x,y,z) and exports it to an STL file.
%Assumes that the mesh is cyclic (see comment below)

%Einar Heiberg 2014-08-01.

%---For debugging purposes create the mesh
%o = linspace(0,2*pi,20);
%o = repmat(o,[10 1]);
%x = sin(o);
%y = cos(o);
%z = 1:10;
%z = repmat(z(:),[1 20]);
%x = x.*z;
%y = y.*z;
%pathname = pwd;
%filename = 'test.stl';
%resolution = 3;
%figure(45);
%mesh(x,y,z);
%closeapex = false;

if (~ischar(pathname)) && (nargin==7)
  %Not a pathname, assume it is a fid.
  fidsupplied = true;
  fid = pathname;
else
  if nargin<8
    myfailed('Too few input arguments.');
    return;
  end;
  fidsupplied = false;
end;

if isempty(fignr)
  doplot = false;
else
  doplot = true;
  c = rand(1,1);  
end;

if size(x,1)<2
  disp('Only one slice, ignored'); 
end;

%If should plot
if doplot
  figure(fignr);
  hold on;
end;

%Open the file
if ~fidsupplied
  fid = fopen([pathname filesep filename],'w');
  if isequal(fid,-1)
    myfailed(dprintf('Could not open the file %s for writing.',[pathname filename]));
    return;
  end;
end;

%Start the file
if ~fidsupplied
  fprintf(fid,'solid segment stl \n');
end;

n = 1:size(x,2); %Vector 1-n elements along the contour

for row = 1:(size(x,1)-1)
  
  %Extract row1
  xrow1 = x(row,:);
  yrow1 = y(row,:);
  zrow1 = z(row,:);
  
  %Extract row2
  xrow2 = x(row+1,:);
  yrow2 = y(row+1,:);
  zrow2 = z(row+1,:);  
    
  %Calc length of row1 and row2
  drow1 = sum(sqrt((xrow1(2:end)-xrow1(1:(end-1))).^2 + (yrow1(2:end)-yrow1(1:(end-1))).^2 + (zrow1(2:end)-zrow1(1:(end-1))).^2));
  drow2 = sum(sqrt((xrow2(2:end)-xrow2(1:(end-1))).^2 + (yrow2(2:end)-yrow2(1:(end-1))).^2 + (zrow2(2:end)-zrow2(1:(end-1))).^2));
  
  %Make sure row2 is the longest
  if drow1>drow2
    %Switch them
    oxrow2 = xrow2;
    oyrow2 = yrow2;
    ozrow2 = zrow2;
    odrow2 = drow2;
    
    xrow2 = xrow1;
    yrow2 = yrow1;
    zrow2 = zrow1;
    drow2 = drow1;
    
    xrow1 = oxrow2;
    yrow1 = oyrow2;
    zrow1 = ozrow2;
    drow1 = odrow2;
    
    switched = true;
  else
    switched = false;
  end;
  
  %Resample row1 according to length 
  newn = round(drow1/resolution);
  if newn<1
    newn = 1;
  end;
  if newn>size(x,2)
    newn = size(x,2);
  end;
  ni = linspace(1,size(x,2)-(size(x,2)-1)/newn,newn); %As the mesh is cyclic it should not interpolate to the end of the line.
  xrow1 = interp1(n,xrow1,ni);
  yrow1 = interp1(n,yrow1,ni);
  zrow1 = interp1(n,zrow1,ni);
  
  %Resample row2 according to length
  newn = round(drow2/resolution);
  if newn<1
    newn = 1;
  end;
  if newn>size(x,2)
    newn = size(x,2);
  end;
  ni = linspace(1,size(x,2)-(size(x,2)-1)/newn,newn);
  xrow2 = interp1(n,xrow2,ni);
  yrow2 = interp1(n,yrow2,ni);
  zrow2 = interp1(n,zrow2,ni);  
  
  closestpoint = round(linspace(1,length(xrow1),length(xrow2)));
  
  %Define vector of next indices
  nextind = [1:length(xrow2) 1]; %This is to close the circles
  for loop = 1:length(xrow2)
    
    %Extract point 1
    x1 = xrow2(nextind(loop));
    y1 = yrow2(nextind(loop));
    z1 = zrow2(nextind(loop));

    %Extract point 2 (one row up)
    x2 = xrow1(closestpoint(nextind(loop)));
    y2 = yrow1(closestpoint(nextind(loop)));
    z2 = zrow1(closestpoint(nextind(loop)));    
    
    %Extract point 3 (to the right)
    x3 = xrow2(nextind(loop+1));
    y3 = yrow2(nextind(loop+1));
    z3 = zrow2(nextind(loop+1));    
    
    %Write triangle 1->2->3
    normal = -cross([x2-x1,y2-y1,z2-z1],[x3-x1,y3-y1,z3-z1]);
    if switched
      normal = -normal;
    end;
    fprintf(fid, 'facet normal %f %f %f \n', normal );
    fprintf(fid, 'outer loop \n');
    fprintf(fid, 'vertex %f %f %f \n', x1,y1,z1);
    fprintf(fid, 'vertex %f %f %f \n', x2,y2,z2);
    fprintf(fid, 'vertex %f %f %f \n', x3,y3,z3);
    fprintf(fid, 'endloop \n');
    fprintf(fid, 'endfacet \n');    
    
    if doplot
      %Graphical update
      patch([x1;x2;x3],[y1;y2;y3],[z1;z2;z3],c*ones(3,1));
    end;
      
    %Check if not redundant triangle
    if ~isequal(closestpoint(nextind(loop)),closestpoint(nextind(loop+1)))
      
      %Extract point 4 (up and to the right)
      x4 = xrow1(closestpoint(nextind(loop+1)));
      y4 = yrow1(closestpoint(nextind(loop+1)));
      z4 = zrow1(closestpoint(nextind(loop+1)));
      
      %Write triangle 2->4->3
      normal = -cross([x2-x4,y2-y4,z2-z4],[x4-x4,y4-y3,z4-z3]);
      if switched
        normal = -normal;
      end;
      fprintf(fid, 'facet normal %f %f %f \n', normal );
      fprintf(fid, 'outer loop \n');
      fprintf(fid, 'vertex %f %f %f \n', x2,y2,z2);
      fprintf(fid, 'vertex %f %f %f \n', x4,y4,z4);
      fprintf(fid, 'vertex %f %f %f \n', x3,y3,z3);
      fprintf(fid, 'endloop \n');
      fprintf(fid, 'endfacet \n');
    
      if doplot
        %Graphical update
        patch([x2;x4;x3],[y2;y4;y3],[z2;z4;z3],c*ones(3,1));
      end;
      
    end; %If not this triangle is the same as the previous one
    
  end; %loop over number of points in the largest layer

  %Close apex
  if isequal(row,size(x,1)-1) && closeapex
    %Last row and we should close the apex
    %Row2 is the last row.

    %Recalculate if it was previously switched
    xrow2 = x(row+1,:);
    yrow2 = y(row+1,:);
    zrow2 = z(row+1,:);
    drow2 = sum(sqrt((xrow2(2:end)-xrow2(1:(end-1))).^2 + (yrow2(2:end)-yrow2(1:(end-1))).^2 + (zrow2(2:end)-zrow2(1:(end-1))).^2));
    
    %Check if already closed.
    if round(drow2/resolution)>1
      %Not already closed
      
      %--- Close it
      
      %Resample according to length
      newn = round(drow2/resolution);
      if newn<1
        newn = 1;
      end;
      if newn>size(x,2)
        newn = size(x,2);
      end;
      ni = linspace(1,size(x,2)-(size(x,2)-1)/newn,newn);
      xrow2 = interp1(n,xrow2,ni);
      yrow2 = interp1(n,yrow2,ni);
      zrow2 = interp1(n,zrow2,ni);
      
      %Define centerpoint
      x2 = mean(xrow2);
      y2 = mean(yrow2);
      z2 = mean(zrow2);
      
      %Loop over points
      nextind = [1:newn 1]; %This is to close the circles
      for loop = 1:newn
        
        %Extract point 1
        x1 = xrow2(nextind(loop));
        y1 = yrow2(nextind(loop));
        z1 = zrow2(nextind(loop));
                
        %Extract point 3 (to the right)
        x3 = xrow2(nextind(loop+1));
        y3 = yrow2(nextind(loop+1));
        z3 = zrow2(nextind(loop+1));
        
        %Write triangle 1->2->3
        normal = cross([x2-x1,y2-y1,z2-z1],[x3-x1,y3-y1,z3-z1]);
        fprintf(fid, 'facet normal %f %f %f \n', normal );
        fprintf(fid, 'outer loop \n');
        fprintf(fid, 'vertex %f %f %f \n', x1,y1,z1);
        fprintf(fid, 'vertex %f %f %f \n', x2,y2,z2);
        fprintf(fid, 'vertex %f %f %f \n', x3,y3,z3);
        fprintf(fid, 'endloop \n');
        fprintf(fid, 'endfacet \n');
        
        if doplot
          %Graphical update
          patch([x1;x2;x3],[y1;y2;y3],[z1;z2;z3],c*ones(3,1));
        end;
      
      end; %Loop over points when closing
    end; %If indeed to be closed
  end; %Close apex and last row 
end; %Loop over rows

%End the file
if ~fidsupplied  
  fprintf(fid,'endsolid stl\n');
  fclose(fid);
end;

if doplot
  hold off;
end;
