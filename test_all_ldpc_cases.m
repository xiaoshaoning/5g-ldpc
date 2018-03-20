base_graph_list = 1:2;

Z_c_list = [2, 4, 8, 16, 32, 64, 128, 256, ...
3, 6, 12, 24, 48, 96, 192, 384, ...
5, 10, 20, 40, 80, 160, 320, ...
7, 14, 28, 56, 112, 224, ...
9, 18, 36, 72, 144, 288, ...
11, 22, 44, 88, 176, 352, ...
13, 26, 52, 104, 208, ...
15, 30, 60, 120, 240];

pass_or_failure = 1;
for base_graph = base_graph_list
    for Z_c = Z_c_list
        pass_flag = test_ldpc(base_graph, Z_c);
        if ~pass_flag
          pass_or_failure = 0;
          disp(['The case base_graph = ', num2str(base_graph), ' Z_c= ', num2str(Z_c), ' passed.']);
          break;
        else
          disp(['case base_graph = ', num2str(base_graph), ' Z_c= ', num2str(Z_c), ' passed.']);  
        end
    end
end

if pass_or_failure == 1
   disp('test passed.');
else
   disp('test failed.');
end
