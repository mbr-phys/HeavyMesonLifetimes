# HeavyMesonLifetimes
Matrix Elements for Heavy Meson Lifetimes from Lattice QCD

Software source code based on [Grid](https://github.com/paboyle/Grid) and [Hadrons](https://github.com/aportelli/Hadrons), which are both free software under GPLv2. Both Grid & Hadrons need to be installed as dependencies, and production binaries were using Grid up to commit XXXX and Hadrons up to XXXX.

The code can be compiled using the sequence of commands below
```
./bootstrap.sh
mkdir build
cd build
../configure --with-grid=<Grid prefix> --with-hadrons=<Hadrons prefix> CXX=<compiler used for Grid>
make
```
