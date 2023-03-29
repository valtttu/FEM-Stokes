function [A,B] = assembly_stokes(mesh,T,E)
    %%% TODO
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
    res = @(x) 1 / triangle_area(mesh,T) * n(j) * edge_length(mesh,E);
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
    absT = det([edge_vector1, edge_vector2]);
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
    M = [2/3  1/6  1/6;
         1/6  2/3  1/6;
         1/6  1/6  2/3];
    qpts = [(M*corner_coord(1,:)')' ; (M*corner_coord(2,:)')'];
    A = a(qpts(:,1))' * b(qpts(:,1)) +...
        a(qpts(:,2))' * b(qpts(:,2)) +...
        a(qpts(:,3))' * b(qpts(:,3));
    res = A * triangle_area(mesh,T) / 3;
end