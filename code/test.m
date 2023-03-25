%% Code for testing the mesh

clear all;
close all;

%% Initialize a mesh
% Option 1: Rect mesh on (0,1)^2
%mesh = make_rect_mesh(1);
% Option 2: A more complex mesh from input file
mesh = build_mesh('blocks.txt',2);

% The number of elements
Nt = size(mesh.t, 2);
Ne = size(mesh.edges, 2);
n = Nt + Ne; % The number of dof for u

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
% Loop over the elements (interiors in tau_h)
A = zeros([2*n, 2*n]);
B = zeros([2*n, Nt]);
idof = mesh.idof;
edof = mesh.edof;
for i = 1:Nt
    
    for k = 1:2 % Both components of u
        % Compute interior contribution to A from T_i
        A(idof(i,k), idof(i,k)) = A(idof(i,k), idof(i,k)) + 1; % nabla phi_i * nabla phi_i
    
        % Compute edge contribution to A from T_i and edges of T_i
        for j = 1:3
            A(idof(i,k), edof(i,j,k)) = A(idof(i,k), edof(i,j,k)) + 1; % nabla phi_i * nabla psi_j
            A(edof(i,j,k), idof(i,k)) = A(edof(i,j,k), idof(i,k)) + 1; % nabla psi_i * nabla phi_j

            % Compute edge edge contribution to A
            A(edof(i,j,k), edof(i,j,k)) = A(edof(i,j,k), edof(i,j,k)) + 1; % nabla psi_i * nabla psi_i
        end

    end

    % Compute the contributions to B (edge and interior contribution only)
    % only for k=1, since p is scalar
    for j = 1:3
        B(mesh.edof(i,j,1), mesh.idof(i,1)) = B(mesh.edof(i,j,1), mesh.idof(i,1)) + 1; % (nabla * psi_i) phi_j
    end

end

% Construct the saddle point system
bdof = mesh.bdof;
bvals = mesh.bvals;
iidof = setdiff(1:2*n, bdof);
M = [A, -B; B',zeros(Nt)];
b = zeros([2*n+Nt,1]);
u(iidof) = M(iidof, iidof)\b(iidof);
u(bdof) = bvals;

