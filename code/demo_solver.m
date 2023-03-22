
% make the mesh
mesh = make_rect_mesh(1);

% Specify the bilinear form
bilin = @(U,V,dU,dV,gX)(dU{1}.*dV{1} + 0.001*dU{2}.*dV{2} ); % K = [1 0 ; 0 0.001]
linf = @(V,dV,gX)(V); % f = 1.

% assembly matrices 
[Ahat,bhat] = simple_assembly(mesh,bilin,linf);

% find boundary and interior nodes 
be = find( mesh.e2t(2,:) == 0);
bind = mesh.edges(:,be);
bind = unique(bind(:));
iind = setdiff(1:size(mesh.p,2),bind);

% solve the problem with zero Dirichlet bc.
u =zeros( size(mesh.p,2),1);
u(iind) = Ahat(iind,iind)\bhat(iind);

% plot the solution 
X = mesh.p(1,:);
Y = mesh.p(2,:);
t = mesh.t;
patch(X(t),Y(t),u(t),u(t));


