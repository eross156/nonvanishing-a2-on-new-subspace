
# Nonvanishing Second Coefficient on the New Subspace

In the paper [Nonvanishing of Second Coefficients of Hecke Polynomials on the New Subspace](https://arxiv.org/abs/2407.11694), we showed that `a2^new(2,N,k)` and  `a2^new(4,N,k)` are nonvanishing for sufficiently large `N+k`.
This repository contains code to check the finitely many remaining cases.

- The computation for `a2^new(2,N,k)` was split into 120 batches: `A0` - `A119`.
  - `Ai` computed all `i*10^8 <= N < (i+1)*10^8`.
- The computation for `a2^new(4,N,k)` was split into 61 batches: `B1` - `B61`.
  - `Bi` computed `N=i` for `1 <= i <= 60`.
  - `B61` computed `N >= 61`.
- `Table A` gives the complete list of triples `(2,N,k)` for which `a2^new(2,N,k) >= 0`.
- `Table B` gives the complete list of triples `(4,N,k)` for which `a2^new(4,N,k) <= 0`.


Here is an example of how the code was run.
```
$ # verify that our a2 formula matches what is computed by Sage
$ sage compute_a2.sage test
$ # check batch A0
$ sage compute_a2.sage A0 > out_A0.txt
...
$ # double check that our computations give the same list 
$ # from just checking all the small values of N and k.
$ sage compute_dim.sage search_small_N_k > out_small_N_k.txt
```
