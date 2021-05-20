function outstri = dictionary(stri,language,fromlanguage)
% DICTIONARY
% Translate a string (in English) into a specified language. Returns the
% English string if no translation is available

global DATA
persistent labels 

if nargin < 3
  fromlanguage = 'English';
end

if nargin < 2
  if ~isempty(DATA) && isprop(DATA,'Pref') && isfield(DATA.Pref,'Language')
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
  if isdeployed
    transpath = 'dictionary.mat';
  else
    transpath = fullfile([DATA.SegmentFolder filesep '+translation'],'dictionary.mat');
  end  
  load(transpath,'labels')
end

dictind = find(strcmp(stri,{labels(:).(fromlanguage)}),1);
% if ~isempty(dictind) && isfield(labels,language) && ...
%     ~isempty(labels(dictind).(language))
%   outstri = labels(dictind).(language);
% end

% new implementation that if the value does not exist in one language it is
% searched in the english dictionary since it might be an enlish wording
if isempty(dictind)
  % try to find index in the english value
  dictind = find(strcmp(stri,{labels(:).English}),1);
end
if ~isempty(dictind)
  %translate if the value is found either in fromlanguage or in the english
  if isfield(labels,language) && ~isempty(labels(dictind).(language))
    % set value to the found value in the language
    outstri = labels(dictind).(language);
  else
    % set value to the english one
    outstri = labels(dictind).English;
  end
end

