# Gonorrhea-modeling-1
Nneed download Armadillo C++ library from [Armadillo](http://arma.sourceforge.net/download.html).

 
Complie the source code as follows:

```
g++ .\Model_optimistic.cpp -o .\Model_optimistic.exe -std=c++14
```

There are six inputs for the program.
* running id
* recommended level of testing (F/T)
* contact tracing of regular partners (F/T)
* contact tracing of casual partners (F/T)
* probability to test a regular partner
* probability to test a casual partner


Run the EXE file as follows:
```
.\Model_optimistic.exe 123 F F F 0.0 0.0
```
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5002980.svg)](https://doi.org/10.5281/zenodo.5002980)
