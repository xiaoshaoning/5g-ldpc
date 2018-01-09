function [H] = alist2sparse(fname) 
% reads binary parity check matrix in "alist" format from file FNAME and 
% converts it to sparse matrix used in MATLAB routines. 
% This is an interface to matrices at http://wol.ra.phy.cam.ac.uk/mackay/codes/ 
% 
% Example 
%        [H] = alist2sparse('A');   % A is the ascii file in alist format 
 
 
%   Copyright (c) 1999 by Igor Kozintsev igor@ifp.uiuc.edu 
%   $Revision: 1.1 $  $Date: 2000/03/23 $ Bug fixed by Hatim Behairy 
 
fid = fopen(fname); 
n = fscanf(fid,'%d',1); 
m = fscanf(fid,'%d',1); 
maxinrow = fscanf(fid,'%d',1);  
junk = fscanf(fid,'%d',1); % no need 
num = fscanf(fid,'%d',[1 n]); % number of elements in rows 
 
num2(1:n)=maxinrow;     
junk = fscanf(fid,'%d',[1 m]); % no need 
 
position = zeros(n,maxinrow); 
for i=1:n 
   for j=1:num2(i)     
      position(i,j) = fscanf(fid,'%d',1); 
   end 
end 
 
ii = zeros(1,sum(num)); 
jj = ii; 
k = 1; 
for i=1:n 
      for j=1:num(i) 
      jj(k) = i; 
      ii(k) = position(i,j); 
      ss = 1; 
      k = k+1 ;  
   end 
end 
H = sparse(ii,jj,ss,m,n); 
fclose(fid);