% Antti Hannukainen / 11.1.2008 / Otaniemi
%
% plot_2Dtri_mesh.m
%
%
% Plot 2D triangular mesh. The plotting is done by using line or patch -
% functions.
%
%
% The calling syntax of plot_2Dtri_mesh is
%
% H = plot_2Dtri_mesh(mesh,method) 
%
% mesh         = mesh-structure to be plotted (must be 2D triangles)
%
% method       =  method can be used to choose between lines / patch. "patch" is the default. It can produce bad
%                 quality .eps - files, but it is faster.  "lines" option is recommended for
%                 exporting .eps - files
%
% RETURNS
% 
% H            = Handle to drawn graphic object ( can be used to change
%                                                 color etc. see set(H) )
%



function H=plot_2Dtri_mesh(mesh,method) 

if(nargin < 3)
    material = 0;
end

if(nargin ==1)
    method = '';
end

if( strcmp(method,'lines'))

    X(1,:) = mesh.p(1,mesh.edges(1,:));
    X(2,:) = mesh.p(1,mesh.edges(2,:));
    
    Y(1,:) = mesh.p(2,mesh.edges(1,:));
    Y(2,:) = mesh.p(2,mesh.edges(2,:));
   
    H = line(X,Y);
else
    
    px = mesh.p(1,:);
    py = mesh.p(2,:);
    
    X = px(mesh.t);
    Y = py(mesh.t);
    
    H= patch(X,Y,0*X);set(H,'FaceColor','none');

end


