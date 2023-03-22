% HJ-FEM - finite element package : sc_find.m
%
%----------------------------------------------------------------- 
%
% Brama version Antti Hannukainen 27.10.2008, Otaniemi 
%
%-----------------------------------------------------------------
% 
% This function takes as an argument a matrix of integers A and B and returns
% unique index I, A(:,I) = B.  Columns are treated as sets, i.e. ordering of
% column elements is irrelevant. 
% 
% If element of B is not included in A, index is set to zero
%
%
% CALLING SYNTAX IS
%
% function I = sc_unique(A,B)
%
%   A,B     = input matrices, both are assumed as integers.
%
%   I       = index vector, A(:,I) = B (naturally, "=" is up to ordering)
%
%
% SEE ALSO
%   sc_unique 


function I = sc_find(A,B)

% columnwise sort
B = sort(B,1);

% make sure, that dublicates in B will be correctly handled
[uB,useless,uIB] = unique(B','rows');

[D,iA,iB] = intersect(A',uB,'rows');

uI = repmat(-1,1,size(uB,2));

uI(iB) = iA;
I = uI(uIB);
