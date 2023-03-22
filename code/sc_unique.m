% HJ-FEM - finite element package : sc_unique.m
%
%----------------------------------------------------------------- 
%
% Antti Hannukainen 7.8.2007, Otaniemi
%
% Revised by Antti Hannukainen 9.11.2007 Oulunkylä
% Revised by Antti Hannukainen 8.7.2008  Aachen
% Brama revision by Antti Hannukainen 27.10.2008
%
%-----------------------------------------------------------------
% 
% This function takes as an argument a matrix of integers A and returns
% unique columns B and index I, B(I) = A.  Columns are treated as
% sets, i.e. ordering of column elements is irrelevant
%
%
% CALLING SYNTAX IS
%
% function [B,I,J] = sc_unique(A)
%
%   A       = the original matrix. All elements of this matrix
%             are assumed as integers
%              
%   B       = unique columns of A. columns are treated as sets, so
%             ordering is irrelevant
%   
%   J,I     = index vectors, B = A(I), A(J) = B; (naturally, "=" is up to ordering)
%
%
% SEE ALSO
%   sc_find


function [B,I,J] = sc_unique(A)

% sort columnwise
A = sort(A,1);

% unique
[B,I,J] = unique(A','rows');

% transpose to columns + "unsort"
B = B';

I = I(:)';
J = J(:)';