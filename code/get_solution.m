% Function for plotting solutions

% Parameters:
% u         Solution given by the linear system
% mesh      The mesh object defining the domain etc.
% Returns:
% X         Vector of x coordinates
% Y         Vector of y coordinates
% U         Velocity x-components
% V         Velocity y-components
% p         Pressure at the given locations

function [X, Y, U, V, p] = get_solution(u, mesh)

    % The number of elements
    Nt = size(mesh.t, 2);
    Ne = size(mesh.edges, 2);
    n = Nt + Ne; % The number of dof for u
   
    % Loop over the elements and compute coordinates and components
    X = [];
    Y = [];
    U = [];
    V = [];
    p = u(end-Nt:end); % This way order of u correponds to element inds

    for ti = 1:Nt
        % Interior -> centroid
        edge_coord = mesh.p(:,mesh.t(:,ti));
        c = mean(edge_coord,2); % parameter 2 gives a row-wise mean
        X(end+1) = c(1);
        Y(end+1) = c(2);
        U(end+1) = u(mesh.idof(ti,1));
        V(end+1) = u(mesh.idof(ti,2));

    end

    for ti = 1:Nt
        % Edges -> mid-point of the edge
        for ei = 1:3
            eind = mesh.edges(:,mesh.t2e(ei,ti));
            cp = mean([mesh.p(:,eind(1)), mesh.p(:,eind(2))], 2);
            X(end+1) = cp(1);
            Y(end+1) = cp(2);
            U(end+1) = u(mesh.edof(ti,ei,1));
            V(end+1) = u(mesh.edof(ti,ei,2));
        end
    end


end