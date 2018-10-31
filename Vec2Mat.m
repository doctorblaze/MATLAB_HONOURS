%%-------------------------------------------------------------------------
% This script takes a 6 element vector containing rotations (in degrees) along x, y and z 
% and translations along x y and z
% It returns an affine 4x4 transformation matrix M
%
% e.g. [M] = Vec2Mat([140 7 2 -2 10 -14]); 
%   note: if using rots and trans from VINCI for elastix, invert signs because 
%   VINCI and elastix have different x, y, z orientation conventions (see maskmake.m notes)
%
% Author: Giorgos Angelis
%%-------------------------------------------------------------------------

function [M] = Vec2Mat(vec)

    rx = vec(1)*pi/180; % input in degrees
    ry = vec(2)*pi/180;
    rz = vec(3)*pi/180;
    tx = vec(4); %input in mm
    ty = vec(5);
    tz = vec(6);
    
    
    T1 = [1 0 0 tx;
          0 1 0 ty;
          0 0 1 tz;
          0 0 0  1];

    T2 = [       1        0         0     0; %ROTATIONS ALONG X
                 0  cos(rx)  -sin(rx)     0;
                 0  sin(rx)   cos(rx)     0;
                 0        0         0     1];

    T3 = [cos(ry)         0   sin(ry)     0; %ROTATIONS ALONG Y
                0         1         0     0;
         -sin(ry)         0   cos(ry)     0;
                0         0         0     1];

    T4 = [cos(rz) -sin(rz)          0     0; %ROTATIONS ALONG Z
          sin(rz)  cos(rz)          0     0;
                0        0          1     0;
                0        0          0     1];

    [M] = T1*T4*T3*T2*eye(4);                     

    N=M(1:3,1:3)';
    transformixFormat = [N(:)' M(1:3,4)']
    
end

