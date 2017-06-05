clear all;  close all;  clc;

addpath ~/myCode/src/tensor_toolbox_2.5/   % sandia tensor toolbox
addpath ~/myCode/repos/ctd/                % my library

test_paths;   if ~ok; return;  end;
test_norms;   if ~ok; return;  end;
test_getterms;  if ~ok;  return;  end;
test_trncsval; if ~ok;  return;  end;
test_trncspmat;  if ~ok;  return;  end;
test_poswts;  if ~ok;  return;  end;
test_trncel;  if ~ok;  return;  end;
test_gram;    if ~ok;  return;  end;
test_k;    if ~ok;  return;  end;
test_krand;    if ~ok;  return;  end;
test_vvpwm;    if ~ok;  return;  end;
test_mvm;   if ~ok; return;  end;
%test_mmm;   if ~ok;  return; end;
%test_als;   if ~ok;  return; end;
test_tenid;  if ~ok;  return;  end;