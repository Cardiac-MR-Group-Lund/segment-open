classdef zoomfunctions
  %#ok<*GVMIS>
  methods (Static)

    %-----------------------------------
    function imdim = getxysize(no,panel)
      %---------------------------------
      % function to get image dimensions that are visible in the panel
      global SET DATA
      originalheight = SET(no).XSize;
      originalwidth = SET(no).YSize;
      switch DATA.ViewPanelsType{panel}
        case {'montage','montagerow'}

          zoomstate = SET(no).NormalZoomState;
          if isempty(zoomstate)
            % init it to the whole image
            imdim = zoomfunctions.resetimagedimensions(originalheight,originalwidth);
            return
          end
%           axpos = getpixelposition(DATA.Handles.imageaxes(panel));
%           axwidth = axpos(3);
%           axheight = axpos(4);
%           viewmatrix = DATA.ViewPanelsMatrix{panel};
%           imheight = axheight/viewmatrix(2);
%           imwidth = axwidth/viewmatrix(1);

          %calculate zoom state center point
          yc = mean(zoomstate(1:2));
          xc = mean(zoomstate(3:4));

          %find new zoom region adapted to montage view, using center point
          width = zoomstate(2) - zoomstate(1);
          ystart = max(floor(yc - width/2),1);
          yend = min(ceil(yc + width/2),originalwidth);
          newwidth = yend-ystart+1;

          height = zoomstate(4) - zoomstate(3);
          xstart = max(floor(xc - height/2),1);
          xend = min(ceil(xc + height/2),originalheight);
          newheight = xend-xstart+1;

          if newheight  >= originalheight && newwidth  >= originalwidth
            % reset zoom state and image dimension to the whole image
            SET(no).NormalZoomState = [];
            imdim = zoomfunctions.resetimagedimensions(originalheight,originalwidth);
            return
          end

          if newheight < originalheight
            % height was changed -> set width to the same proportion
            newproportion = newheight/originalheight;
            newwidth = newproportion*originalwidth;
            %adjust zoom window
            ystart = max(floor(yc - newwidth/2),1);
            yend = min(ceil(yc + newwidth/2),originalwidth);
            newwidth = yend-ystart+1;
          elseif newwidth < originalwidth
            % check if width was changed -> set height to the same proportion
            newproportion = newwidth/originalwidth;
            newheight = newproportion*originalheight;
            %adjust zoom window
            xstart = max(floor(xc - newheight/2),1);
            xend = min(ceil(xc + newheight/2),originalheight);
            newheight = xend-xstart+1;
          end

          imdim.XStart = xstart;
          imdim.XEnd = xend;
          imdim.XSize = newheight;

          imdim.YStart = ystart;
          imdim.YEnd = yend;
          imdim.YSize = newwidth;

        otherwise
          imdim = zoomfunctions.resetimagedimensions(originalheight,originalwidth);
      end
    end

    %------------------------------------------------------------------
    function imdim = resetimagedimensions(originalheight,originalwidth)
      %----------------------------------------------------------------
      % reset dimensions to origianl width and height
      imdim.XStart = 1;
      imdim.XEnd = originalheight;
      imdim.XSize = originalheight;

      imdim.YStart = 1;
      imdim.YEnd = originalwidth;
      imdim.YSize = originalwidth;
    end

    %------------------------------------------------------------------
    function imdim = getdefaultimagedimensions(no)
      %----------------------------------------------------------------
      global SET
      originalheight = SET(no).XSize;
      originalwidth = SET(no).YSize;
      imdim = zoomfunctions.resetimagedimensions(originalheight,originalwidth);
    end

    %-----------------------------------------------------------
    function [xt,yt] = gettranslationforslice(no,panel,sliceind)
      %---------------------------------------------------------
      % function to get translation in x and y direction, based on the
      % slice index in montage view (other view types tranlation is 0)
      global DATA
      % row and column number are switched, since DATA.ViewPanelsMatrix has
      % first number of columns and then number of rows stored
      [colnum,rownum] = ind2sub(DATA.ViewPanelsMatrix{panel},sliceind);
      if isempty(colnum) || isempty(rownum)
        rownum = 1;
        colnum = 1;
      end

      imdim = zoomfunctions.getxysize(no,panel);
      xt = (rownum-1)*imdim.XSize - imdim.XStart +1;
      yt = (colnum-1)*imdim.YSize - imdim.YStart +1;
    end

    %------------------------------------------------------------------------------
    function [contourx,contoury] = maskoutzoomedcontour(no,panel,contourx,contoury)
      %----------------------------------------------------------------------------

      imdim = zoomfunctions.getxysize(no,panel);
      margin = 0.5;
      contourx(contourx > imdim.XEnd+margin) = nan;
      contourx(contourx < imdim.XStart-margin) = nan;
      contoury(contoury > imdim.YEnd+margin) = nan;
      contoury(contoury < imdim.YStart-margin) = nan;
    end
  end


end