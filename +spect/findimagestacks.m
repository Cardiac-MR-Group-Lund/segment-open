function [restno,stressno,restgatedno,stressgatedno] = findimagestacks
%Find rest and stress image stacks, gated and ungated.

global SET

%find stress and rest image stack number
restno = [];
stressno = [];
restgatedno = [];
stressgatedno = [];
for noloop = 1:length(SET)
  if strcmp(SET(noloop).ImageType,'Perfusion Rest')
    if SET(noloop).TSize > 1
      restgatedno = [restgatedno noloop]; %#ok<AGROW>
    else
      restno = [restno noloop]; %#ok<AGROW>
    end
  elseif strcmp(SET(noloop).ImageType,'Perfusion Stress')
    if SET(noloop).TSize > 1
      stressgatedno = [stressgatedno noloop]; %#ok<AGROW>
    else
      stressno = [stressno noloop]; %#ok<AGROW>
    end
  end
end