function [h, g] = ldpc_h2g(H,q) 
% [h, g] = ldpc_h2g(H,q) 
% converts tentative binary LDPC matrix H into a new matrix h 
% (columns are permuted) and produces the generator matrix g 
% H should be a sparse matrix in MATLAB format. 
% q - Field base (power of 2) now only 2 4 8 1 32 64 128 and 256 
% 
% MEX file 
 
 
%   Copyright (c) 1999 by Igor Kozintsev igor@ifp.uiuc.edu 
%   $Revision: 1.1 $  $Date: 1999/08/23 $ - implementation for GFq 