function h = GetAlt()
%
%
persistent altData     
persistent idx init


if isempty(init)
  load AltData
  idx = 1;
  
  init = 1;
end

h = altData(idx);
  
idx = idx + 1;