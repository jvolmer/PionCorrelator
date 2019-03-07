# Pion Correlator

Compute the neutral Pion Correlator for already computed configurations (configurations are too large to put here, you have to create your own to use this code. I used configurations computed via the tmLQCD-code, <github.com/etmc/tmLQCD>).

The computation is done in two steps:

## 1. Invert Dirac matrix
Using e.g. 4 stochastic sources
```cpp
./run_dis.x p 4
```
The result, the propagator, is saved in output/propagator/. The used stochastic sources are stored in output/source/.

## 2. Combine the inverted matrices to the correlator 
```cpp
./run_dis.x c 4
```
The correlator is saved in output/correlator/.