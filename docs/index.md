<div align="center">
  <img src="images/logo.png" alt="LetzElPhC Logo" style="width:40%;">
  <p style="font-size:0.8em; color: gray; margin-top: 5px;">
    Logo designed by Henry Fried
  </p>
</div>

# Welcome to LetzElPhC Documentation

This is the official documentation for the LetzElPhC project.
Here you will find installation instructions, usage examples.

## About the Code

LetzElPhC is a C code designed to compute electron-phonon coupling matrix elements from the outputs of standard Density Functional Theory (DFT) and Density Functional Perturbation Theory (DFPT) calculations.  
Currently, it only supports Quantum Espresso, with plans for Abinit.  
The main objective is to facilitate electron-phonon calculations within Yambo (version 5.2+).  
The code is released under the MIT license and hosted on GitHub: [link](https://github.com/muralidhar-nalabothula/LetzElPhC).

### Main Features

- Utilizes full crystal symmetries, ensuring compatibility with Yambo code
- Multiple levels of parallelization (OpenMP, plane-wave, k-point, q-point)
- Fully parallel I/O via NetCDF-4/HDF5
- Highly portable across CPU architectures and OS
