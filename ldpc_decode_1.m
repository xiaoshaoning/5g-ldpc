function [x_hat, success, k] = ldpc_decode_1(f,H,qq) 
% decoding of LDPC over GFqq, qq = 2,4,8,16,32,64,128 and 256 
% as in Comm. Letters by Davey&MacKay June 1998 with e few modifications. 
% For notations see the same reference. 
% outputs the estimate "x_hat" of the ENCODED sequence for 
% the received vector with channel likelihoods "f". 
% "f" ([2^qq][n]) stores the likelihoods for "n" symbols in natural  
% ordering. E.g., y(3,5) is the probability of 5-th symbol is equal to "2". 
% "H" is the parity check matrix. Success==1 signals 
% successful decoding. Maximum number of iterations is set to 100. 
% k returns number of iterations until convergence. 
% 
% Examples: 
% We assume G is systematic G=[A|I] and G*H'=0 over GFq 
% Binary case 
%         sigma = 1;                          % AWGN noise deviation 
%         x = (sign(randn(1,size(G,1)))+1)/2; % random bits 
%         y = mod(x*G,2);                     % encoding  
%         z = 2*y-1;                          % BPSK modulation 
%         z=z + sigma*randn(1,size(G,2));     % AWGN transmission 
% 
%         f1=1./(1+exp(-2*z/sigma^2));        % likelihoods 
%         f1 = (f1(:))';                      % make it a row vector 
%         f0=1-f1; 
%         [z_hat, success, k] = ldpc_decode([f0;f1],H,2); 
%         x_hat = z_hat(size(G,2)+1-size(G,1):size(G,2)); 
%         x_hat = x_hat';  
% 
% Nonbinary case 
%         sigma = 1;                          % AWGN noise deviation 
%         q = 4;                              % Field parameter 
%         nbits = log2(q);                    % bits per symbol 
%         h = ldpc_generate(400,600,2.5,q,123); % Generate H 
%         [H,G] = ldpc_h2g(h,q);              % find systematic G and modify H 
%         x = floor(rand(1,size(G,1))*q);     % random symbols 
%         y = ldpc_encode(x,G,q);             % encoding  
%         yb = (fliplr(de2bi(y,nbits)))';     % convert total index to binary format 
%         yb = yb(:);                         % make a vector 
%         zb = 2*yb-1;                        % BPSK modulation 
%         zb=zb + sigma*randn(size(zb));      % AWGN transmission 
% 
%         f1=1./(1+exp(-2*zb/sigma^2));        % likelihoods for bits 
%         f1 = f1(:);                         % make it a vector 
%         f1 = reshape(f1,nbits,length(y));   % reshape for finding priors on symbols                     
%         f0=1-f1; 
%         junk = ones(q,length(y));           % this is a placeholder in the next function 
%         [v0, v1, pp] = bits_smbl_msg(f0,f1,junk); 
%         [z_hat, success, k] = ldpc_decode(pp,H,q); 
%         x_hat = z_hat(size(G,2)+1-size(G,1):size(G,2)); 
%         x_hat = x_hat';  
 
 
%   Copyright (c) 1999 by Igor Kozintsev igor@ifp.uiuc.edu 
%   $Revision: 1.2 $  $Date: 1999/11/23 $ 
%   fixed high-SNR decoding 
%   works for GFq, q= 2^m now 
 
if qq==2 % binary case first, just use the old code 
    
   [m,n] = size(H); if m>n, H=H'; [m,n] = size(H); end 
   if ~issparse(H) % make H sparse if it is not sparse yet 
      [ii,jj,sH] = find(H); 
      H = sparse(ii,jj,sH,m,n); 
   end 
    
   f0 = f(1,:); % prob of 0 
   f1 = f(2,:); 
    
   %initialization 
   [ii,jj,sH] = find(H);          % subscript index to nonzero elements of H  
   indx = sub2ind(size(H),ii,jj); % linear index to nonzero elements of H 
   q0 = H * spdiags(f0(:),0,n,n); 
   sq0 = full(q0(indx));  
   sff0 = sq0; 
 
   q1 = H * spdiags(f1(:),0,n,n);  
   sq1 = full(q1(indx)); 
   sff1 = sq1; 
 
   %iterations 
   k=0; 
   success = 0; 
   max_iter = 100; 
   while ((success == 0) & (k < max_iter)), 
      k = k+1; 
    
      %horizontal step 
      sdq = sq0 - sq1; sdq(find(sdq==0)) = 1e-20; % if   f0 = f1 = .5 
      dq = sparse(ii,jj,sdq,m,n); 
      Pdq_v = full(real(exp(sum(spfun('log',dq),2)))); % this is ugly but works :) 
      Pdq = spdiags(Pdq_v(:),0,m,m) * H; 
      sPdq = full(Pdq(indx)); 
      sr0 = (1+sPdq./sdq)./2; sr0(find(abs(sr0) < 1e-20)) = 1e-20; 
      sr1 = (1-sPdq./sdq)./2; sr1(find(abs(sr1) < 1e-20)) = 1e-20; 
      r0 = sparse(ii,jj,sr0,m,n); 
      r1 = sparse(ii,jj,sr1,m,n); 
    
      %vertical step 
      Pr0_v = full(real(exp(sum(spfun('log',r0),1)))); 
      Pr0 = H * spdiags(Pr0_v(:),0,n,n); 
      sPr0 = full(Pr0(indx)); 
      Q0 = full(sum(sparse(ii,jj,sPr0.*sff0,m,n),1))'; 
      sq0 = sPr0.*sff0./sr0; 
    
      Pr1_v = full(real(exp(sum(spfun('log',r1),1)))); 
      Pr1 = H * spdiags(Pr1_v(:),0,n,n); 
      sPr1 = full(Pr1(indx));  
      Q1 = full(sum(sparse(ii,jj,sPr1.*sff1,m,n),1))'; 
      sq1 = sPr1.*sff1./sr1; 
    
      sqq = sq0+sq1; 
      sq0 = sq0./sqq; 
      sq1 = sq1./sqq; 
    
      %tentative decoding 
      QQ = Q0+Q1; 
      Q0 = Q0./QQ; 
      Q1 = Q1./QQ; 
    
      x_hat = (sign(Q1-Q0)+1)/2; 
      if rem(H*x_hat,2) == 0, success = 1; end 
   end 
   % end of binary case 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
else % GFq, nonbinary   
   % our strategy is "divide and concur" - we partition H into several matrices with 
   % the fixed number of variables per function in each of them and the other way around 
    
   [m,n] = size(H); if m>n, H=H'; [m,n] = size(H); end 
   if ~issparse(H) % make H sparse if it is not sparse yet 
      [ii,jj,sH] = find(H); 
      H = sparse(ii,jj,sH,m,n); 
   end 
    
   %initialization 
   [ii,jj,sH] = find(H);          % subscript index to nonzero elements of H  
   W = sparse(ii,jj,ones(size(ii)),m,n); %indicator function 
   nvars = full(sum(W,2));        % number of variables participating each check function 
   minvars = min(nvars);          % min number of variables in a function 
   maxvars = max(nvars);          % max number of variables in a function 
    
   nfuns = full(sum(W,1));        % number of functions per variable 
   minfuns = min(nfuns);          % min number of functions per variable 
   maxfuns = max(nfuns);          % max number of functions per variable 
    
   % the following will be used in solving linear equations over GFq 
   M=log2(qq); % GFq exponent 
   [tuple power] = gftuple([-1:2^M-2]', M, 2);  
   alpha = tuple * 2.^[0 : M - 1]'; 
   beta(alpha + 1) = 0 : 2^M - 1; 
 
 
   % create cell arays which contain sparse matrices with fixed # of variables in rows 
   for nnvars = minvars:maxvars 
      tmp = zeros(size(H)); 
      rows = find(nvars == nnvars); %rows of H having 'nnvars' variables 
      tmp(rows,:) = H(rows,:); 
      [jjj,iii,ssH] = find(tmp');  
      iir{nnvars} = reshape(iii,nnvars,length(iii)/nnvars)'; 
      jjr{nnvars} = reshape(jjj,nnvars,length(jjj)/nnvars)'; 
      Hr{nnvars} = reshape(ssH,nnvars,length(ssH)/nnvars)';% separate parity matrices 
      q{nnvars} = reshape(f(:,jjr{nnvars})',[size(jjr{nnvars}),qq]); %initialize to channel likelihoods 
       
      % Prestore valid configurations in array X 
      if(~isempty(Hr{nnvars})) % make sure the are functions for this case 
         Hleft = Hr{nnvars}(:,1);         % will solve for these varibles 
         Hright = Hr{nnvars}(:,2:nnvars); % while setting these arbitrary 
         for i=0:(qq^(nnvars-1)-1)  % there are qq^(nnvars-1) different combinations 
            xr = (fliplr(de2bi(i,nnvars-1,qq))); % current nonzero combination  
              
            % find the remaining variable to satisfy the parity checks 
            right_part = ones(size(Hleft))*(-Inf); %exponent over GFq 
            for j=1:(nnvars-1) % multiply each column of Hright by the symbol from x and accumulate 
               rr1 = power(beta(xr(j)+1)+1);      % get expon. representation of xr(i) 
               rr2 = power(beta(Hright(:,j)+1)+1);% same for the column of Hright 
               rr3 = gfmul(rr1,rr2,tuple)'; % this is exponential representation of the product 
               right_part = gfadd(right_part,rr3,tuple); 
            end 
            left_part = mod((qq-1)+ right_part - power(beta(Hleft+1)+1),qq-1); 
            xl=zeros(size(left_part)); 
            nzindx = find(isfinite(left_part)); 
            xl(nzindx) = alpha(left_part(nzindx)+2); 
            x = [xl repmat(xr,[length(xl),1])]; %this is a valid configuration 
            X{nnvars}(i+1,:,:) = x; 
         end 
      end 
       
       
   end 
    
   % create cell arays which contain sparse matrices with fixed # of functions in columns 
   for nnfuns = minfuns:maxfuns 
      tmp = zeros(size(H)); 
      cols = find(nfuns == nnfuns); %rows of H having 'nnvars' variables 
      tmp(:,cols) = H(:,cols); 
      [iii,jjj,ssH] = find(tmp);  
      iic{nnfuns} = reshape(iii,nnfuns,length(iii)/nnfuns); 
      jjc{nnfuns} = reshape(jjj,nnfuns,length(jjj)/nnfuns); 
      Hc{nnfuns} = reshape(ssH,nnfuns,length(ssH)/nnfuns);% separate parity matrices 
      ff{nnfuns} = reshape(f(:,jjc{nnfuns})',[size(jjc{nnfuns}),qq]); %  this will not change 
   end 
 
   %iterations 
   k=0; 
   success = 0; 
   max_iter = 100; 
   while ((success == 0) & (k < max_iter)), 
      k = k+1 
       
      buffer = zeros([size(H),qq]); 
       
      % Horizontal step - forming messages to variables from the parity check functions 
      % each Hr is processed separately 
      for nnvars = minvars:maxvars 
         if(~isempty(Hr{nnvars})) % make sure the are functions for this case 
         result = zeros([size(Hr{nnvars}) qq]); % will store the intermediate result 
         for i=0:(qq^(nnvars-1)-1)  % there are qq^(nnvars-1) different combinations 
               x = squeeze(X{nnvars}(i+1,:,:)); %lookup a valid configuration 
                
               %calculate products 
               a = cumsum(ones(size(x)),1); 
               b = cumsum(ones(size(x)),2); 
               idx = sub2ind(size(q{nnvars}),a,b,x+1); %index of current configuration in 3D 
               pp = repmat(prod(q{nnvars}(idx),2),[1,size(x,2)]); %product for this configuration  
         
               denom = q{nnvars}(idx); 
               denom(find(denom==0)) = realmin; 
               result(idx) = result(idx) + pp./denom; 
            end 
             
            %  update global distribution 
            a = repmat(iir{nnvars},[1,1,qq]); 
            b = repmat(jjr{nnvars},[1,1,qq]); 
            c = permute(repmat((1:qq)',[1 size(a,1) size(a,2)]),[2 3 1]); 
            gidx = sub2ind(size(buffer),a,b,c); 
            buffer(gidx) = result; 
         end 
      end 
       
      % initialize r from the global data in buffer 
      for nnfuns = minfuns:maxfuns 
         a = repmat(iic{nnfuns},[1,1,qq]); 
         b = repmat(jjc{nnfuns},[1,1,qq]); 
         c = permute(repmat((1:qq)',[1 size(a,1) size(a,2)]),[2 3 1]); 
         gidx = sub2ind(size(buffer),a,b,c); 
         r{nnfuns} = buffer(gidx); 
      end 
 
       
      %vertical step 
      buffer = zeros([size(H),qq]); 
      QQ = zeros(qq,size(H,2)); 
      for nnfuns = minfuns:maxfuns 
         if(~isempty(Hc{nnfuns})) % make sure the are variables for this case 
            %calculate products 
            pp = repmat( prod ( r{nnfuns},1),[size(r{nnfuns},1),1]).*ff{nnfuns}; %product for this configuration  
            denom = r{nnfuns}; 
            denom(find(denom==0)) = realmin; 
            result = pp./denom; 
            result = result./repmat((sum(result,3)),[1,1,qq]); %normalize to distribution 
            %  update global distribution 
            a = repmat(iic{nnfuns},[1,1,qq]); 
            b = repmat(jjc{nnfuns},[1,1,qq]); 
            c = permute(repmat((1:qq)',[1 size(a,1) size(a,2)]),[2 3 1]); 
            gidx = sub2ind(size(buffer),a,b,c); 
            buffer(gidx) = result; 
             
            Q{nnfuns} = pp.*ff{nnfuns}; 
            b = repmat(jjc{nnfuns}(1,:),[qq,1]); 
            c = repmat((1:qq)',[1, size(b,2)]); 
            qidx =  sub2ind(size(QQ),c,b); 
            QQ(qidx) = squeeze(Q{nnfuns}(1,:,:))'; 
         end 
      end 
       
 
       
      % initialize q from the global data in buffer 
      for nnvars = minvars:maxvars 
         a = repmat(iir{nnvars},[1,1,qq]); 
         b = repmat(jjr{nnvars},[1,1,qq]); 
         c = permute(repmat((1:qq)',[1 size(a,1) size(a,2)]),[2 3 1]); 
         gidx = sub2ind(size(buffer),a,b,c); 
         q{nnvars} = buffer(gidx); 
      end 
    
       
      %tentative decoding 
      QQ = QQ ./ repmat(sum(QQ,1),[qq 1]); %normalize - can be used as soft outputs 
      [xi xj sx] = find(QQ == repmat(max(QQ),[size(QQ,1),1])); 
      x_hat = xi-1; 
      if ldpc_encode(x_hat,H',qq) == 0, success = 1; end 
   end 
end % end of nonbinary case 