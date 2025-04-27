function h = GetAccZ()
%
%
persistent acc_z     
persistent idx init


if isempty(init)
  load acc_z.mat
  idx = 1;
  
  init = 1;
end

h = acc_z(idx);
  
idx = idx + 1;