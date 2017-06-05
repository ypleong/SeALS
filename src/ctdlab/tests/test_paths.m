clear all;

ok = 1;

fprintf('\n***test_paths***\n')
try
  tenrand([3,3,3]);
  fprintf('tensor toolbox : ok\n')
catch
  fprintf('tensor toolbox : path not found.\n')
  fprintf('FAILED\n')
  ok = 0;
  return
end
  
try
  krandn(3,5,1,1);
  fprintf('ctd  :  ok\n')
catch
  fprintf('ctd  :  path not found.\n')
  fprintf('FAILED\n')
  ok = 0;
  return
end

fprintf('PASSED\n')