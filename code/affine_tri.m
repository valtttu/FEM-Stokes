% HJ-FEM - finite element solver : affine_tet.m
% 
% Antti Hannukainen and Mika Juntunen 1.3.2007 
% 
% affine_tet.m constructs vectorized affine mappings.
%
% The vectorized affine mapping consists out of three matrices,
% Ax, and Ay and two vectors bx and by. Row i in Ax corresponds
% to affine mapping of component x from the reference element to
% element i. 
%
% The classical affine mapping from reference coordinates (x,y,z)
% to element i is 
%
%  Ax(i,:) *  x  +  bx(i)
%  Ay(i,:)    y     by(i)
%
% The structure is designed to support operations for
% elements simultaniously. In practice this means, that a mapping
% of integration points from reference element to global element
% can be done simply as
%
% Ax *  x1 x2 x3 ...   + bx
%       y1 y2 y3 ...
%       z1 z2 z3 ... 
%
% 
% The calling syntax of affine_tri is
%
% function [Ax,Ay,bx,by,detA,Px,Py] = affine_tri(mesh)
%
% mesh         = mesh-structure corresponding to the fname. 
% 
%


function [Ax,Ay,bx,by,detA,Px,Py] = affine_tri(mesh)



Ax = [(mesh.p(1,mesh.t(2,:))-mesh.p(1,mesh.t(1,:)))' ...
      (mesh.p(1,mesh.t(3,:))-mesh.p(1,mesh.t(1,:)))'];

Ay = [(mesh.p(2,mesh.t(2,:))-mesh.p(2,mesh.t(1,:)))' ...
      (mesh.p(2,mesh.t(3,:))-mesh.p(2,mesh.t(1,:)))'];


bx = mesh.p(1,mesh.t(1,:))';
by = mesh.p(2,mesh.t(1,:))';


% Determinant
%-------------
detA = -Ax(:,2).*Ay(:,1) + Ax(:,1).*Ay(:,2) ;

% Transpose of inverse of A
%----------------------------
Px = bsxfun(@rdivide,[ Ay(:,2) -Ay(:,1) ],detA);
Py = bsxfun(@rdivide,[ -Ax(:,2)  Ax(:,1)],detA);







