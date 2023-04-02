function [A,B,F] = assembly_stokes(mesh, f)
    % The number of elements
    Nt = size(mesh.t, 2);
    Ne = size(mesh.edges, 2);
    n = Nt + Ne; % The number of dof for u

    A = sparse(2*n, 2*n);
    B = sparse(2*n, Nt);
    F = zeros([2*n + Nt, 1]);
    idof = mesh.idof;
    edof = mesh.edof;
    for i = 1:Nt
        GPhi = grad_phi(mesh, i);
        inner1 = inner_prod(GPhi, GPhi, mesh, i);
        for k = 1:2 % Both components of u
            % Compute interior contribution to A from T_i
            A(idof(i,k), idof(i,k)) = A(idof(i,k), idof(i,k)) + inner1; % nabla phi_i * nabla phi_i
        
            % Compute edge contribution to A from T_i and edges of T_i
            for j = 1:3
                GPsi = grad_psi(mesh, i, mesh.t2e(j, i));
                inner2 = inner_prod(GPhi, GPsi, mesh, i);
                inner3 = inner_prod(GPsi, GPhi, mesh, i);
                inner4 = inner_prod(GPsi, GPsi, mesh, i);
                A(idof(i,k), edof(i,j,k)) = A(idof(i,k), edof(i,j,k)) + inner2; % nabla phi_i * nabla psi_j
                A(edof(i,j,k), idof(i,k)) = A(edof(i,j,k), idof(i,k)) + inner3; % nabla psi_i * nabla phi_j
    

                % Compute edge edge contribution to A
                A(edof(i,j,k), edof(i,j,k)) = A(edof(i,j,k), edof(i,j,k)) + inner4; % nabla psi_i * nabla psi_i
            end

            % Compute the RHS
            F(idof(i,k)) = inner_prod(f{k}, @(x) ones(size(x,2)), mesh, i);

        end
        
        % Compute the contributions to B (edge and interior contribution only)
        % only for k=1, since p is scalar
        for j = 1:3
            for k = 1:2
                DPsi_j = div_psi(mesh, i, mesh.t2e(j, i), k);
                inner5 = inner_prod(DPsi_j, @(x) ones(size(x,2)), mesh, i);
                B(edof(i,j,k), idof(i,1)) = B(edof(i,j,k), idof(i,1)) + inner5; % (nabla * psi_i) phi_j
            end
        end
    end
end

% TODO: Assembly loop

% Assembly happens about like this:

% GP1 = grad_psi(mesh,T1,E1)
% GP2 = grad_psi(mesh,T2,E2)
% A_{i,j} = inner_prod(GP1,GP2,mesh,T1)

% GP1 and GP2 are function handles, and T1, T2, E1, E2 correspond to basis functions i and j

function res = div_phi()
    res = @(x) 0;
end

function res = div_psi(mesh,T,E,j) % non-zero only for triangles that share edge e_i
    % j \in {1,2} gives the velocity's different components:
    % j=1 gives x, j=2 gives y
    n = edge_normal(mesh,T,E);
    res = @(x) 0*norm(x) + 1 / triangle_area(mesh,T) * n(j) * edge_length(mesh,E);
end

function res = grad_phi(mesh,T) % triangle inside
    x_T = centroid(mesh,T);
    C_T = @(x) 2*triangle_area(mesh,T) / norm(x-x_T)^2;
    res = @(x) -C_T(x) * (x-x_T);
end

function res = grad_psi(mesh,T,E)
    x_T = centroid(mesh,T);
    absT = triangle_area(mesh,T);
    C_T = @(x) 2*absT / norm(x-x_T)^2;
    res = @(x) C_T(x)/3 * (x-x_T) + edge_length(mesh,E)/absT * edge_normal(mesh,T,E);
end

function c = centroid(mesh,T)
    edge_coord = mesh.p(:,mesh.t(:,T));
    c = mean(edge_coord,2); % parameter 2 gives a row-wise mean
end

function absT = triangle_area(mesh,T)
    edge_coord = mesh.p(:,mesh.t(:,T));
    edge_vector1 = edge_coord(:,2) - edge_coord(:,1);
    edge_vector2 = edge_coord(:,3) - edge_coord(:,1);
    absT = abs(det([edge_vector1, edge_vector2]));
end

function outward_normal = edge_normal(mesh,T,E)
    edge_coord = mesh.p(:,mesh.edges(:,E));
    % find the third point, which is not on the edge
    % start by list of all triangle indices, and remove E one by one
    triangle_third_coord = mesh.t(mesh.t(:,T) ~= mesh.edges(1,E),T);
    triangle_third_coord = triangle_third_coord(triangle_third_coord ~= mesh.edges(2,E));
    triangle_third_coord = mesh.p(:,triangle_third_coord);
    % evaluate some projections to find inward normal
    edge_vector = edge_coord(:,2) - edge_coord(:,1);
    edge_to_third_vector = triangle_third_coord - edge_coord(:,1);
    inward_normal = edge_to_third_vector - edge_to_third_vector' * edge_vector / (edge_vector' * edge_vector) * edge_vector;
    outward_normal = -inward_normal / norm(inward_normal);
end

function len = edge_length(mesh,E)
    edge_coord = mesh.p(:,mesh.edges(:,E));
    len = norm(edge_coord(:,2) - edge_coord(:,1));
end

function res = inner_prod(a,b,mesh,T) % int_T a' * b dx 
    % simple center point integral
    % could made better with quadrature rules?
    corner_coord = mesh.p(:,mesh.t(:,T));
    qpts = corner_coord*[2/3  1/6  1/6;
                         1/6  2/3  1/6;
                         1/6  1/6  2/3];
    A = a(qpts(:,1))' * b(qpts(:,1)) +...
        a(qpts(:,2))' * b(qpts(:,2)) +...
        a(qpts(:,3))' * b(qpts(:,3));
    res = A * triangle_area(mesh,T) / 3;
end

function res = hannukainen_quad2(a, b, mesh, T)
    % Calculate quad points and weights
    vertices = mesh.p(:, mesh.t(:, T));
    next = [vertices(:, 2:end), vertices(:, 1)];
    t = vertices + 1/2 * (next-vertices);
    w = triangle_area(mesh, T)*[1/3 1/3 1/3];
    % Operate
    res = w*(a(t).*b(t))';
end