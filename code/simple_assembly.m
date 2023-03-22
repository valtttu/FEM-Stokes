% Antti Hannukainen 1.10.2010 / Otaniemi
%
%
% Assembly routine for P1-elements
%
% The bilinear form is given as a structure bilin :
%
% bilin = function specifying the bilinear form. The input for
%                  this function is :
%
%                  bilin(U,V,dU,dV,gX)
%
%                  - U,V  basisfunction values
%
%                   - dU,dV = cell array of function derivatives
%                     dU{i} = xi - derivative
%
%                  - global coordinate gX in a cell array
%
% linf = function specifying the load functional. The input for
%        this function is :
%
%                   linf(V,dVx,dVy,gX).
%
% ALL INPUT WILL BE IN MATRIX FORM.
%


function [K,F] = simple_assembly(mesh,bilin,linf)

Ndof = size(mesh.p,2);

% initialize
[Ax,Ay,bx,by,detA,Px,Py] = affine_tri(mesh);

[X,W] = inttri(3);

gX{1} = bsxfun(@plus,Ax*X,bx);
gX{2} = bsxfun(@plus,Ay*X,by);



L{1} = 1-X(1,:)-X(2,:);
L{2} = X(1,:);
L{3} = X(2,:);

Nip = size(X,2);
Nt = size(mesh.t,2);

dL{1} = [ -ones(1,Nip); -ones(1,Nip) ];
dL{2} = [  ones(1,Nip); zeros(1,Nip) ];
dL{3} = [ zeros(1,Nip);  ones(1,Nip) ];


% define arrays for matrix entries.
ffind = [];
ff = zeros(3,Nt);

iind = [];
jind = [];
kk = zeros(9,Nt);

mind = 1;

for i=1:3
    
    Li = repmat(L{i},Nt,1);
    
    dLi{1} = Px*dL{i};
    dLi{2} = Py*dL{i};
    
    ff(i,:) = linf(Li,dLi,gX)*W.*abs(detA);
    

    
    for j=1:3
        
        Lj = repmat(L{j},Nt,1);
        
        
        dLj{1} = Px*dL{j};
        dLj{2} = Py*dL{j};
               
        % Keep track of indeces : can be eliminated ??  !! ??
        iind = [iind i] ; jind = [jind j];
        
        % SIMPLIFIED HERE !!!!
        kk(mind,:) = bilin(Lj,Li,dLj,dLi,gX)*W.*abs(detA);
        mind = mind+1;
    end
end

K = sparse(mesh.t(iind,:),mesh.t(jind,:),kk,Ndof,Ndof);

F = sparse(mesh.t,ones(size(mesh.t)),ff,Ndof,1);

