function outstri = dictionary(stri,language,fromlanguage)
% DICTIONARY
% Translate a string (in English) into a specified language. Returns the
% English string if no translation is available

global DATA
persistent labels %#ok<USENS>

if nargin < 3
  fromlanguage = 'English';
end

if nargin < 2
  if ~isempty(DATA) && isfield(DATA.Pref,'Language')
    language = DATA.Pref.Language;
  else
    language = 'English';
  end
end

outstri = stri;
if strcmp(language,fromlanguage) || size(stri,1) ~= 1
  return
end

if isempty(labels)
  load(fullfile('+translation','dictionary.mat'),'labels')
end

dictind = find(strcmp(stri,{labels.(fromlanguage)}),1);
if ~isempty(dictind) && isfield(labels,language) && ...
    ~isempty(labels(dictind).(language))
  outstri = labels(dictind).(language);
end