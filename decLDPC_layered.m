% =========================================================================
% Title       : LDPC decoder with layered detection schedule
% File        : decLDPC_layered.m
% -------------------------------------------------------------------------
% Description :
%   Layered LDPC detection according to E. Sharon, S. Litsyn, and J.
%   Goldberger 2004. Supports SPA and MPA. Note that this LDPC decoder
%   requires LLRs defined according to L=log(Pr[x=0]/Pr[x=1]).
% -------------------------------------------------------------------------
% Revisions   :
%   Date       Version  Author  Description
%   20-may-11  1.2      studer  cleanup for reproducible research
%   04-jul-07  1.1      studer  some bugfixes
%   02-jul-07  1.0      studer  file created
% -------------------------------------------------------------------------
%   (C) 2006-2011 Communication Theory Group
%   ETH Zurich, 8092 Zurich, Switzerland
%   Author: Dr. Christoph Studer (e-mail: studer@rice.edu)
% =========================================================================

function [bit_output,LLR_D2,NumC,NumV] = decLDPC_layered(TxRx,LDPC,LLR_A2, base_graph_check_node_list, Z_c)
% -- initializations
%   numOfEntries = sum(sum(LDPC.H==1));
%   Rcv = spalloc(LDPC.par_bits,LDPC.tot_bits,numOfEntries); % msg matrix
Rcv = zeros(LDPC.par_bits, LDPC.tot_bits);
LLR_D2 = [zeros(1, 2*Z_c), LLR_A2];
% initialize with input bitslis
%   bit_output = zeros(1, LDPC.inf_bits);
NumC = 0; % number of computed check nodes
NumV = 0; % number of computed variable nodes

% === TEST
%   vMask = ones(LDPC.par_bits,1); % indicates NCU needs work
%   WORKLOAD = 0;

% -- BEGIN loop over LDPC-internal iterations
for iter=1:TxRx.Decoder.LDPC.Iterations
    
    %%% % -- count the number of nonzero entries of H*x
    %     errors = 0;
    
    switch (TxRx.Decoder.LDPC.Type)
        case 'SPA' % -- sum-product algorithm
            % == for all parity check node
            for j=1:LDPC.par_bits
                % idx = find(LDPC.H(j,:)==1); % slow
                idx = base_graph_check_node_list{j};
                S = LLR_D2(idx)-full(Rcv(j,idx));
                Spsi = sum(-log(1e-300+tanh(abs(S)*0.5)));
                Ssgn = prod(sign(S));
                for v=1:length(idx)
                    Qtemp = LLR_D2(idx(v)) - Rcv(j,idx(v));
                    Qtemppsi = -log(1e-300+tanh(abs(Qtemp)*0.5));
                    Qtempsgn = Ssgn*sign(Qtemp);
                    Rcv(j,idx(v)) = Qtempsgn*(-log(1e-300+tanh(abs(Spsi-Qtemppsi)*0.5)));
                    LLR_D2(idx(v)) = Qtemp + Rcv(j,idx(v));
                    % -- count number variable & check nodes computed
                    NumC = NumC + 1;
                    NumV = NumV + 1;
                end
                
            end
        case 'MPA' % -- max-log approximation (max-product algorithm)
            % == for all parity check node
            for j=1:LDPC.par_bits
                %idx = find(LDPC.H(j,:)==1); % slow
                idx = base_graph_check_node_list{j};
                S = LLR_D2(idx)-full(Rcv(j,idx));
                for v=1:length(idx)
                    Stmp = S;
                    Stmp(v) = []; % clear row
                    Rcv(j,idx(v)) = min(abs(Stmp))*prod(sign(Stmp));
                    LLR_D2(idx(v)) = S(v) + Rcv(j,idx(v));
                    % --variable node & check node computed
                    NumC = NumC + 1;
                    NumV = NumV + 1;
                end
            end
        case 'OMS' % -- offset min-sum [Chen05]
            % == for all parity check node
            for j=1:LDPC.par_bits
                idx = find(LDPC.H(j,:)==1); % slow
%                 idx = base_graph_check_node_list{j};
                S = LLR_D2(idx)-Rcv(j,idx);
                Stmp_abs = abs(S);

                [Stmp_abs_min, min_position] = min(Stmp_abs);
                Stmp_abs_without_min = Stmp_abs;
                Stmp_abs_without_min(min_position) = [];
                Stmp_abs_second_min = min(Stmp_abs_without_min);
                               
                Stmp_sign_1 = sign(S);
                Stmp_sign_1(Stmp_sign_1==0) = 1;       
                Stmp_sign_prod = prod(Stmp_sign_1);
                for v=1:length(idx)
                    if v == min_position
                        Magnitude = Stmp_abs_second_min;
                    else    
                        Magnitude = Stmp_abs_min;
                    end

                    Signum = Stmp_sign_prod * Stmp_sign_1(v);
                    
                    Rcv(j,idx(v)) = Signum*Magnitude;
                    LLR_D2(idx(v)) = S(v) + Rcv(j,idx(v));
                    % --variable node & check node computed
                    NumC = NumC + 1;
                    NumV = NumV + 1;
                end
            end
        otherwise
            error('Unknown TxRx.Decoder.LDPC.Type.')
    end
    
end

% -- compute binary-valued estimates
bit_output = 0.5*(1-mysign(LLR_D2(1:LDPC.inf_bits)));

clear Rcv

return

function s = mysign(inp)
s = 2*double(inp>0)-1;
return