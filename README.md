5g ldpc codes
=============

To test the LDPC encoding and decoding functions, run the following function under matlab 
```
test_ldpc
```
To test the LDPC encoding function, execute the below function under matlab,
```
test_ldpc_encode
```

The encoding function was only tested for K = 8448 and K = 3840.
More functions such as rate matching, interleaving will be added in the future.

The LDPC decoding function decLDPC_layered.m is from [Simulator for LDPC decoding in IEEE 802.11n](http://www.csl.cornell.edu/~studer/software_ldpc.html) and the author is Christoph Studer. I only made a small modification on the decoding function.

The two included excel files are from 3gpp.
The matlab codes other than decLDPC_layered.m follow the MIT license.
