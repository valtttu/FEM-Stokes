% Run the whole solver for cavity problem

clear all;
close all;

% Build the mesh
mesh = build_mesh('./domains/cavity_domain.txt',1);

% The number of elements
Nt = size(mesh.t, 2);
Ne = size(mesh.edges, 2);
n = Nt + Ne; % The number of dof for u

% Assemble the matrices
[A, B] = assembly_stokes(mesh);


% Construct the saddle point system
bdof = mesh.bdof;
bvals = mesh.bvals;
iidof = setdiff(1:(2*n + Nt), bdof);
M = [A, -B; B',eye(Nt)];
b = zeros([2*n+Nt,1]);
u(iidof) = M(iidof, iidof)\b(iidof);
u(bdof) = bvals;


% Plot the solution
[X, Y, U, V, p] = get_solution(u, mesh);
figure()
quiver(X, Y, U, V);