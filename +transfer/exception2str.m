function str = exception2str(e)
  str = '';
  str = [str sprintf('Error message:\n%s\n\n',e.message) ];
  str = [str sprintf('Stack trace:\n') ];
  for i=1:numel(e.stack)
    str = [str sprintf('Line %d in ''%s''\n', e.stack(i).line, e.stack(i).file) ];
  end
  str = [str sprintf('\nEnd\n') ];
end