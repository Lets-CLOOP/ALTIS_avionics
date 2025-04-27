function h = GetGyro()
%
%
persistent gyroMat     
persistent idx init


if isempty(init)
  load gyroMat
  idx = 1;
  
  init = 1;
end

h = gyroMat(idx);
  
idx = idx + 1;