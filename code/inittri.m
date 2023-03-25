%----------------------------------------------------------------- 
%
% Antti Hannukainen and Mika Juntunen 3.11.2005
%
% Vectorization by Antti Hannukainen 29.8.2007
% Revision by Antti Hannukainen 27.10.2008
%
%-----------------------------------------------------------------
%
% Construct mesh - structure from a given triangulation (p,t
% explained in detail alongside the mesh structure) 
% 
% TRIMESH STRUCTURE :
% 
%   p       = nodes in a 2xNp-matrix
%
%   t       = triangles in a 3xNt-matrix.
% 
%   edges   = a matrix of all edges in the mesh. Each column is an edge :
%             [n1 ; n2] where n1 < n2;
%
%   t2e     = a matrix connecting  triangle's and edges's. 
%             Each column corresponds to a triangle and has
%             triangle's edges in the order n1->n2, n2->n3,
%             n1->n3. 
%
%   e2t     = inverse of t2e.
%
%   idof    = a Ntx2 matrix that gives the interior dof indices, use this
%             to fill the A and B matrices and access elements from them
%
%   edof    = a Ntx3x2 matrix that gives the edge dof indices, use this to
%             fill the A matrix
%
% CALLING SYNTAX IS 
%
% function mesh = inittri(p,t)
%
%   p       = nodes
%   t       = triangles
%     
%   mesh    = trimesh structure corresponding to (p,t)
%


function mesh = inittri(p,t)

t = sort(t(1:3,:));

mesh.p  = p;
mesh.t  = t;

% Initalize size variables
Nt = size(t,2);

e = [1 2; 2 3; 1 3]';

edges = [];
for i=1:size(e,2)
  edges = [edges [ sort( [t(e(1,i),:); t(e(2,i),:)],1)]];
end

[mesh.edges,~,mesh.t2e] = sc_unique(edges);
mesh.t2e = reshape(mesh.t2e,Nt,3)';

% mesh.e2t
e = [mesh.t2e(1,:) mesh.t2e(2,:)  mesh.t2e(3,:)];
t = repmat([1:Nt],1,3);

[ef,If]= unique(e,'first');
[el,Il]= unique(e,'last');

mesh.e2t(1,ef) = t(If);
mesh.e2t(2,el) = t(Il);

mesh.e2t(2,find( (mesh.e2t(1,:)-mesh.e2t(2,:))==0))=0;


% mesh.idof
mesh.idof = reshape(1:2*Nt, [Nt, 2]);

% mesh.edof
edof = zeros([Nt,3,2]);
edof(:,:,1) = 2*Nt + mesh.t2e';
edof(:,:,2) = 2*Nt + max(mesh.t2e, [], 'all') + mesh.t2e';
mesh.edof = edof;

