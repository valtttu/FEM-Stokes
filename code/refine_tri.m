% HJ-FEM - finite element package : refine_tri.m
%
%----------------------------------------------------------------- 
%
% Antti Hannukainen 30.8.2007 Oulunkylä
%
% Revised by Antti Hannukainen 10.1.2008 Otaniemi
% Brama revision by Antti Hannukainen 27.10.2008
%
%-----------------------------------------------------------------
%
% 
% Uniformly refine triangular mesh.
%
% CHILD STRUCTURE
%
%   child.t_parent = parent element
%   child.e_parent = parent edge
%
% If the edge / face have no parent, we use index 0.
%
%
% CALLING SYNTAX IS 
%
% function [rmesh] = refine_tri(mesh)
%
%           
%

function [rmesh,ch] = refine_tri(mesh)

t = mesh.t;
p = mesh.p;
e = mesh.edges;
t2e = mesh.t2e;

Nt = size(t,2);
Ne = size(e,2);

for n=1:2
  e_nodes(n,:) = sum( [p(n,e(1,:)) ; p(n,e(2,:))])/2;
end

% Create new mesh
rmesh.p = [p  e_nodes];


% Create New elements
eb = size(p,2); 

% Edges as n1->n2, n2->n3, n1->n3.
new_t = [t(1,:) ; t2e(1,:)+eb ; t2e(3,:)+eb ];
rmesh.t = [new_t];

new_t = [t(2,:) ; t2e(1,:)+eb ; t2e(2,:)+eb ];
rmesh.t = [rmesh.t  new_t];

new_t = [t(3,:) ; t2e(3,:)+eb ; t2e(2,:)+eb ];
rmesh.t = [rmesh.t  new_t];

new_t = [t2e(1,:)+eb ; t2e(2,:)+eb ; t2e(3,:)+eb ];
rmesh.t = [rmesh.t  new_t];


% Initialize new mesh
rmesh=inittri(rmesh.p,rmesh.t);

% create child - structure
ch.t_parent = repmat(1:Nt,1,4);

e1 = [ e(1,:) ; [1:Ne]+eb];
e2 = [ e(2,:) ; [1:Ne]+eb];

ch.e_child(1,:) = sc_find(rmesh.edges,e1);
ch.e_child(2,:) = sc_find(rmesh.edges,e2);





  
  
