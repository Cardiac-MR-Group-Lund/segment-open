function outdcms = findsliceuniquetimeframes(mydcms)

%unqdcms = [];
%mydir = dir([dicompath filesep '*.dcm']);
%mydcms = {mydir.name};
ttt = zeros(1,numel(mydcms));
numdcms = numel(mydcms);
wb = mywaitbarstart(numdcms,'Searching for files that do not match slice/phase.');
for i = 1:numdcms
  mytags = segdicomread_mex({[mydcms{i}]});
  ttt(i) = str2double(char(mytags{1}.AcquisitionTime));
  wb = mywaitbarupdate(wb);
end
mywaitbarclose(wb);

unqt = unique(ttt);
thist = hist(ttt,unqt);
unqs = find(thist < max(thist));
unqdcms = {};
for i = unqs
  if numel(i) == 1
    unqdcms = [unqdcms mydcms(ttt == unqt(i))];
  end
end

outdcms = setdiff(mydcms,unqdcms);