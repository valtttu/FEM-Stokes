%
% Antti Hannukaine / 23.2.2011 / Otaniemi
%
% return the length of the longest edge in the mesh.
%

function h=give_h(mesh)

if( size( mesh.p,1) == 2)

h = max( sqrt( sum( (mesh.p(1,mesh.edges(1,:)) - mesh.p(1,mesh.edges(2,:))).^2 +  (mesh.p(2,mesh.edges(1,:)) - mesh.p(2,mesh.edges(2,:))).^2,1)));

end

if( size( mesh.p,1) == 3)

h = max( sqrt( sum( (mesh.p(1,mesh.edges(1,:)) - mesh.p(1,mesh.edges(2,:))).^2 ... 
                 +  (mesh.p(2,mesh.edges(1,:)) - mesh.p(2,mesh.edges(2,:))).^2 ...
                 +  (mesh.p(3,mesh.edges(1,:)) - mesh.p(3,mesh.edges(2,:))).^2 ,1)));

end
