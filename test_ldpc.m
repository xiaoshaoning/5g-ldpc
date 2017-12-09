% 5g ldpc encoding and decoding test
% reference: 3GPP TS 38.212 section 5.3.2
% author: Xiao, Shaoning 萧少宁
% license: MIT
% currently the supported cases are K = 8448 and K = 3840

function test_ldpc(base_graph_index)

if nargin == 0
  base_graph_index = 1;    
end

SNR_list = 15:15;

trial_number = 1;

BLER = zeros(1, length(SNR_list));

TxRx.Decoder.LDPC.Iterations = 10;
% TxRx.Decoder.LDPC.Type = 'OMS';
% TxRx.Decoder.LDPC.Type = 'MPA';

if base_graph_index == 1
%   load base_graph_1_check_node_list
%   base_graph_check_node_list = base_graph_1_check_node_list;
  LDPC.inf_bits = 8448;
elseif base_graph_index == 2
%   load base_graph_2_check_node_list
%   base_graph_check_node_list = base_graph_2_check_node_list;
  LDPC.inf_bits = 3840;
else
  error('wrong base graph index.');
end

for SNR_list_index = 1:length(SNR_list)
    
    sigma_square = 10^(-SNR_list(SNR_list_index)/10);
    
    for trial = 1:trial_number
                
        tx_bits = randi([0, 1], LDPC.inf_bits, 1);
        
        [encoded_bits, LDPC.H] = ldpc_encode(tx_bits, base_graph_index);
        
        [LDPC.par_bits, LDPC.tot_bits] = size(LDPC.H);
        
        symbols = 1 - 2 * encoded_bits;
        
        noise = randn(size(symbols));
        
        waveform = symbols + noise * sqrt(sigma_square);
        
        LLR_received = 2 * waveform / sigma_square;
        
%         rx_bits = decLDPC_layered(TxRx, LDPC, LLR_received.', base_graph_check_node_list);
        
        f0=1./(1+exp(-2*LLR_received));
        f1 = 1 - f0;

        rx_bits_prime = ldpc_decode(f0,f1, full(LDPC.H), TxRx.Decoder.LDPC.Iterations);

%         rx_bits = rx_bits(:);
        
%         if ~isequal(tx_bits, rx_bits)
%           BLER(SNR_list_index) = BLER(SNR_list_index) + 1;    
%         end
        
        if ~isequal(rx_bits_prime(1:LDPC.inf_bits), tx_bits)
          BLER(SNR_list_index) = BLER(SNR_list_index) + 1;
        end
        
    end
    
end

fprintf('BLER is %f\n', BLER/trial_number);

end
