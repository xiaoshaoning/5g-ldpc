% 5g ldpc encoding and decoding test
% reference: 3GPP TS 38.212 section 5.3.2
% author: Xiao, Shaoning 萧少宁
% license: MIT
% currently the only supported case is K = 8448

function test_ldpc

SNR_list = 15:15;

trial_number = 1;

BLER = zeros(1, length(SNR_list));

TxRx.Decoder.LDPC.Iterations = 10;
TxRx.Decoder.LDPC.Type = 'OMS';

load H
LDPC.H = H;
[LDPC.par_bits,LDPC.tot_bits] = size(LDPC.H);
LDPC.inf_bits = LDPC.tot_bits - LDPC.par_bits;

for SNR_list_index = 1:length(SNR_list)
    
    sigma_square = 10^(-SNR_list(SNR_list_index)/10);
    
    for trial = 1:trial_number
                
        tx_bits = randi([0, 1], LDPC.inf_bits, 1);
        
        encoded_bits = ldpc_encode(tx_bits);
        
        symbols = 1 - 2 * encoded_bits;
        
        noise = randn(1, length(symbols));
        
        waveform = symbols + noise * sqrt(sigma_square);
        
        LLR_received = 2 * waveform / sigma_square;
        
        rx_bits = decLDPC_layered(TxRx, LDPC, LLR_received);
        
        rx_bits = rx_bits(:);
        
        if ~isequal(tx_bits, rx_bits)
          BLER(SNR_list_index) = BLER(SNR_list_index) + 1;    
        end
        
    end
    
end

fprintf('BLER is %f\n', BLER/trial_number);

end
