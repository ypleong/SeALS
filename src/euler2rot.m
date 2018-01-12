function [ R ] = euler2rot( rx, ry, rz )
% change from Euler321 to a rotation matrix (check)!!!
%
% The rotation is in order rx, ry, rz
%

Rx = [1 0 0;
      0 cos(rx) -sin(rx);
      0 sin(rx) cos(rx)];
  
Ry = [cos(ry) 0 sin(ry);
      0 1 0;
      -sin(ry) 0 cos(ry)];
  
Rz = [cos(rz) -sin(rz) 0;
      sin(rz) cos(rz) 0;
      0 0 1];
  
R = Rx*Ry*Rz;
end

