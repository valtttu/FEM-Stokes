% Run the solver for some F=0 problem from ./domains

clear all;
close all;

% Build the mesh
mesh = build_mesh(['./domains/cavity_domain.txt'],5);

% The number of elements
Nt = size(mesh.t, 2);
Ne = size(mesh.edges, 2);
n = Nt + Ne; % The number of dof for u

% Define the load function (example 1 from the paper)
f = {@(x) 0;...
     @(x) 0};

% Assemble the matrices
[A, B, F, absT] = assembly_stokes(mesh, f);

% Construct the saddle point system - solve with int_omega approach
bdof = mesh.bdof;
bvals = mesh.bvals;
iidof = setdiff(1:(2*n + Nt), bdof);
u = zeros(2*n + Nt,1);
M = [A, -B; B', zeros(Nt);zeros([1,2*n]),absT];
u(iidof) = M([iidof,2*n+Nt+1], iidof)\([F(iidof);0] - M([iidof,2*n+Nt+1],bdof)*bvals');
u(bdof) = bvals;


% Plot the solution and the mesh
[X, Y, U, V, p] = get_solution(u, mesh);
figure()
hold on;
plot_2Dtri_mesh(mesh);
quiver(X, Y, U, V,1);
xlabel('x');
ylabel('y');
xlim([min(X), max(X)]);
ylim([min(Y), max(Y)]);

% Plot the pressure
figure();
hold on;
for i = 1:Nt
    patch(mesh.p(1,mesh.t(:,i)),mesh.p(2,mesh.t(:,i)),p(i));
end
xlabel('x');
ylabel('y');
colorbar();

% Plot the velocity magnitude
figure();
hold on;
for i = 1:Nt
    patch(mesh.p(1,mesh.t(:,i)),mesh.p(2,mesh.t(:,i)),sqrt(u(mesh.idof(i,1))^2 + u(mesh.idof(i,2))^2));
end
xlabel('x');
ylabel('y');
colorbar();

figure(1);