classdef colorfunctions
  % colllection of functions that handle color

  methods(Static)
    %--------------------------------------------------------------------
    function hexcolor = rgb2hex(rgbcolor)
      %--------------------------------------------------------------------
      % convert RGB values to hexadecimal format
      rgbcolor = round(255*(rgbcolor));
      hexcolor = '#';
      for n = 1:3
        hexvalue = dec2hex(rgbcolor(n), 2);
        hexcolor = [hexcolor, hexvalue]; %#ok<AGROW> 
      end
    end

    %--------------------------------------------------------------------
    function rgbcolor = hex2rgb(hexcolor)
      %--------------------------------------------------------------------
      %Convert an hexadecimal code to a RGB value
      hexcolor(1) = []; %remove hashtag
      rgbcolor = reshape(sscanf(hexcolor,'%2x'),3,[]);
    end

    %--------------------------------------------------------------------
    function brightness = getbrightness(rgbcolor)
      %--------------------------------------------------------------------
      %Return perceived brightness (relative luminance) of an RGB color.
      %Relative luminance of RGB color: L = 0.2126 x R + 0.7152 x G + 0.0722 x B
      brightness = 0.2126 * rgbcolor(1) + 0.7152 * rgbcolor(2) + 0.0722 * rgbcolor(3);
    end

  
  end
end