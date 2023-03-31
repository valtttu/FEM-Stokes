% Run the whole solver for cavity problem

clear all;
close all;

% Build the mesh
mesh = build_mesh('./domains/cavity_domain.txt',3);

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
M = [A, -B; B', 0.0001.*eye(Nt)];
b = [ones([2*n,1]); zeros([Nt,1])]; % Set divergence to zero with latter zeros
u(iidof) = M(iidof, iidof)\b(iidof);
u(bdof) = bvals;


% Plot the solution and the mesh
[X, Y, U, V, p] = get_solution(u, mesh);
figure()
hold on;
quiver(X, Y, U, V);
plot_2Dtri_mesh(mesh);

