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


% Add the boundary conditions
out_inds = find(mesh.e2t(2,:) == 0);
out_edges = mesh.edges(:,out_inds);
out_ts = unique(mesh.e2t(1,out_inds(:)));
bdof = [];
bvals = [];
for i = 1:length(out_ts)
    % Go over the edges of triangle i
    for j = 1:3
        if(sum(mesh.t2e(j, i) == out_inds) == 1)
            eind = mesh.t2e(j, i) == out_inds;
            ei = j;
      
            % Append the set of all outer edge indices
            bdof(end+1) = mesh.edof(i,ei,1);
            bdof(end+1) = mesh.edof(i,ei,2);
    
        end
    end

end

mesh.bdof = bdof;
mesh.bvals = zeros(size(bdof));