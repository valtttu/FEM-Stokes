%% Code for testing the mesh

clear all;
close all;

%% Init a (0,1)^2 rect mesh
mesh = make_rect_mesh(0);

% The number of elements
Nt = size(mesh.t, 2);
Ne = size(mesh.edges, 2);
n = 2*Nt + 2*Ne; % The number of dof for u

% Extract the node coordinates
px = mesh.p(1,:);
py = mesh.p(2,:);

%% Test plotting the interior edges
figure();
for i = 1:Ne
    hold on
    if(mesh.e2t(2,i) ~= 0)
        node = mesh.edges(:,i);
        plot(px(node), py(node), 'k-')
    end
end


%% Demo assembly and indexing
% Loop over the elements
A = zeros([n, n]);
B = zeros([n, Nt]);
for i = 1:Nt
    
    % Compute interior contribution to A from T_i
    A(mesh.idof(i,1), mesh.idof(i,1)) = 1; % u_x
    A(mesh.idof(i,2), mesh.idof(i,2)) = 1; % u_y

    % Compute edge contribution to A from T_i
    for j = 1:3
        A(mesh.edof(i,j,1), mesh.edof(i,j,1)) = A(mesh.edof(i,j,1), mesh.edof(i,j,1)) + 1; % u_x
        A(mesh.edof(i,j,2), mesh.edof(i,j,2)) = A(mesh.edof(i,j,2), mesh.edof(i,j,2)) + 1; % u_y
    end

    % Compute interior to B
    B(mesh.idof(i,1), mesh.idof(i,1)) = 1; % p (scalar)

end

% Construct the saddle point system
M = [A, -B; B',zeros(Nt)];
b = zeros([n+Nt,1]);
u = M\b;
