# Hermite RBF Quasi-interpolation (HRBFQI)

:exclamation: Note: Please note this is a fork of the [official repository](https://github.com/GCVGroup/HRBFQI). 

This is a fast surface reconstruction method from Hermite points, presented in the paper

> Shengjun Liu, Charlie C.L. Wang, Guido Brunnett, and Jun Wang, "A closed-form formulation of HRBF-based surface reconstruction by approximate solution", Computer-Aided Design, Special Issue of 2016 Symposium on Solid and Physical Modeling, June 20-24, 2016, Berlin, Germany, vol.78, pp.147-157, September 2016.

<img src="https://github.com/zishun/HRBFQI/raw/master/Bin/archer.png" width="350"/>


## Build

### C++ Executable

* Windows + MSVS: use ```RBFRecon.sln```. ![](https://img.shields.io/badge/build-passing-brightgreen)
* Linux + Qt Creator: use ```RBFRecon.pro```. ![](https://img.shields.io/badge/build-passing-brightgreen)
* Linux + CMake: ![](https://img.shields.io/badge/build-passing-brightgreen)
    ```sh
    mkdir build && cd build
    cmake .. && make -j9
    ```
    
### Python Binding
* Linux + CMake: requires pybind11 ![](https://img.shields.io/badge/build-passing-brightgreen)
    ```sh
    cd pyHRBFQI
    mkdir build && cd build
    cmake .. && make -j9
    python test_pyHRBFQI.py
    ```
* Windows: ![](https://img.shields.io/badge/build-to_do-yellow)

## Usage
Please follow [ReadMe.txt](https://github.com/zishun/HRBFQI/blob/master/ReadMe.txt)
