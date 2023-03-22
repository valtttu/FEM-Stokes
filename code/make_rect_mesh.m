% Antti Hannukainen / 17.1.2008 / Otaniemi
%
% make_rect_mesh.m
%
% Generate mesh for a (0,1)x(0,1) rectangle
%
% N is the number of refinements performed on mesh
%  


function mesh = make_rect_mesh(N)


%--------------------------------
% Initialize mesh on a rectangle 
%--------------------------------

% Define cornerpoints
p = [ 0 1 1 0 ; 0 0 1 1];

% Define triangles
t = [ 1 2 3; 1 4 3]';

% Define zero dirichlet boundaries (rest are zero Neumann)
b = [ 1 2; 1 4 ; 4 3 ; 2 3]';

% Initialize mesh
mesh = inittri(p,t);

% refine N-times

for i=1:N
  mesh = refine_tri(mesh);
end
