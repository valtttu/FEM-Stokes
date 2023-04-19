% Run the whole solver

clear all;
close all;

% Build the mesh
mesh = build_mesh('./domains/square_domain.txt',4);

% The number of elements
Nt = size(mesh.t, 2);
Ne = size(mesh.edges, 2);
n = Nt + Ne; % The number of dof for u

% Define the load function (example 1 from the paper)
f = {@(x) -4*pi^3.*sin(pi*x(2,:)).*cos(pi*x(2,:)).*(cos(pi*x(1,:)).^2 - sin(pi*x(1,:)).^2) - pi*cos(pi*x(2,:)).*sin(pi*x(1,:));...
     @(x) 4*pi^3.*sin(pi*x(1,:)).*cos(pi*x(1,:)).*(cos(pi*x(2,:)).^2 - sin(pi*x(2,:)).^2) - pi*cos(pi*x(1,:)).*sin(pi*x(2,:))};

% The analytic solution for f
g = {@(x) 2*pi.*sin(pi*x(1,:)).*sin(pi*x(1,:)).*sin(pi*x(2,:)).*cos(pi*x(2,:));...
     @(x) -2*pi.*sin(pi*x(1,:)).*sin(pi*x(2,:)).*sin(pi*x(2,:)).*cos(pi*x(1,:))};

% Assemble the matrices
[A, B, F] = assembly_stokes(mesh, f);


% Construct the saddle point system
bdof = mesh.bdof;
bvals = mesh.bvals;
bvals = bvals;
iidof = setdiff(1:(2*n + Nt), bdof);
u = zeros(2*n + Nt,1);
M = [A, -B; B', (1e-6).*eye(Nt)];
b = [zeros([2*n,1]); zeros([Nt,1])]; % Set divergence to zero with latter zeros
u(iidof) = M(iidof, iidof)\(F(iidof) - M(iidof,bdof)*bvals');
u(bdof) = bvals;


% Plot the solution and the mesh
[X, Y, U, V, p] = get_solution(u, mesh);
figure()
hold on;
plot_2Dtri_mesh(mesh);
quiver(X, Y, U, V, 1);
title('Stokes solver')

figure();
hold on;
plot_2Dtri_mesh(mesh);
quiver(X, Y, g{1}([X(:)';Y(:)']), g{2}([X(:)';Y(:)']),1);
title('Analytic solution');

figure();
hold on;
plot_2Dtri_mesh(mesh);
quiver(X, Y, f{1}([X(:)';Y(:)']), f{2}([X(:)';Y(:)']));
title('Load function');

% Plot the pressure - still needs some work to plot triangles
figure();
imagesc(unique(X(1:Nt)), unique(Y(1:Nt)), reshape(p,[2*sqrt(Nt/2),sqrt(Nt/2)]));
set(gca,'YDir','normal')
xlabel('x');
ylabel('y');
colorbar();
axis tight equal

figure(1);