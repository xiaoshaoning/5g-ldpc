% Name: tanner_LDPC.m
% Description: This script displays different parameters of an LDPC code
% given his A-list format file.
% Copyright (c) 2009. Robert Morelos-Zaragoza. All rights reserved.

% clear all
% 
% name = input('enter file name: ');
name = 'NR_ldpc';
[H] = alist2sparse('NR_ldpc.alist')

figure(1)
spy(H)
title( strcat('H matrix of LDPC code with A-list file name: ', name))

n1 = size(H,1); n2 = size(H,2);

% Adjacency matrix as defined in Matlab
A = [ zeros(n2) H'; H zeros(n1) ];

xy = [];
for i=1:n2
    xy = [ xy; i 2 ];                 % Locations of code nodes
end
for i=1:n1
    xy = [ xy; i 1 ];                 % Locations of check nodes
end

figure(2)
gplot(A,xy)
axis([0 n2+1 0.9 2.1])
title( strcat('Tanner graph of LDPC code with A-list file name: ', name))

figure(3)
for i=1:n1
    cnd(i) = sum(H(i,:));
end
stem(cnd)
title('Check node degree distribution')

figure(4)
for j=1:n2
    vnd(j) = sum(H(:,j));
end
stem(vnd)
title('Variable node degree distribution')

fprintf('Density = %e\n',nnz(H)/(n1*n2));
