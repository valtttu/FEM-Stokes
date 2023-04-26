% Function for debugging outer normals

function debug_normals(mesh, T, E, n)
    % Find the center point of E
    cp = mean([mesh.p(:,mesh.edges(:,E)), mesh.p(:,mesh.edges(:,E))], 2);

    % Plot the normal and element T
    quiver(cp(1), cp(2), n(1), n(2))
    hold on
    xp = [mesh.p(1,mesh.t(:,T)),mesh.p(1,mesh.t(1,T))];
    yp = [mesh.p(2,mesh.t(:,T)),mesh.p(2,mesh.t(1,T))];
    plot(xp, yp, 'k-');
    axis equal

end