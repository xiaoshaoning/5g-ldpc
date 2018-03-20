5g ldpc codes
=============

To test the LDPC encoding and decoding functions, run the following function under matlab 
```
test_all_ldpc_cases
```

The LDPC decoding function decLDPC_layered.m is from [Simulator for LDPC decoding in IEEE 802.11n](http://www.csl.cornell.edu/~studer/software_ldpc.html) and the author is Christoph Studer. I made some small modifications on it to accelerate its execution.

The LDPC decoding function ldpc_decode.m is from https://github.com/willemolding/FastPythonLDPC .

The two included excel files are from 3gpp.
The matlab codes other than decLDPC_layered.m and ldpc_decode.m follow the MIT license.
