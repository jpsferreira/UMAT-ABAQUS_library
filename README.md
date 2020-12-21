
UMAT-ABAQUS_library | A library with user-defined material models for soft biological tissues and filamentous networks. | A framework to test material response under multiple loading conditions.

# Description

This library provides subroutines for 3D implementation of large deformation constitutive behavior using continuum mechanics. The code is written in either fixed-form or modern Fortran. You can use it as a primary subroutine interface for ABAQUS (Dassault Systèmes) as a UMAT or an interface for your own finite element code. 

Each material model comprises: source code and a one-element ABAQUS example. 

This repository also has a framework to test the material response under multiple loading conditions, in particular:
* Load type: Monotonic and cyclic loading 
* Homogeneous deformation mode: Uniaxial, biaxial, simple shear and pure shear

This repository is oriented for experienced researchers in biomechanics and continuum mechanics. You are highly welcome to use these open-resource codes for your non-commercial usage; please cite the below papers. Pull requests are welcome. 

# Compiling

Simple bash scripts at each directory are provided for building with gfortran. The code can also be compiled with the Intel Fortran Compiler (and presumably any other Fortran compiler that supports modern standards). 

# Future Releases

*Makefile
*Benchmarking
### License:
  * GNU GPL 3.0.  


# References

J P S Ferreira et al. "On the mechanical response of the actomyosin cortex during cell indentations". English. In: Biomechanics and Modeling in Mechanobiology 130.10 (Apr. 2020), pp. 202–19. doi: [(10.1007/s10237-020-01324-5)](https://link.springer.com/article/10.1007%2Fs10237-020-01324-5).

J P S Ferreira et al. "Altered mechanics of vaginal smooth muscle cells due to the lysyl oxidase-like1 knockout." English. In: Acta Biomaterialia (Apr. 2020). doi: [(10.1016/j.actbio.2020.03.046)](https://doi.org/10.1016/j.actbio.2020.03.046).

J A López-Campos et al. "Characterization of hyperelastic and damage behavior of tendons." English. In: Computer Methods in Biomechanics and Biomedical Engineering 81.2 (Jan. 2020), pp. 1–11. doi: [(10.1080/10255842.2019.1710742)](https://doi.org/10.1080/10255842.2019.1710742).

M C P Vila Pouca et al. "Viscous effects in pelvic floor muscles during childbirth: A numerical study." English. In: International Journal for Numerical Methods in Biomedical Engineering 34.3 (Mar. 2018), e2927. doi: [(10.1002/cnm.2927)](https://doi.org/10.1002/cnm.2927). 

J P S Ferreira, M P L Parente, and R M Natal Jorge. "Continuum mechanical model for cross-linked actin networks with contractile bundles". English. In: Journal of the Mechanics and Physics of Solids (2017). doi: [(10.1016/j.jmps.2017.09.009)](https://doi.org/10.1016/j.jmps.2017.09.009).

J P S Ferreira et al. "A general framework for the numerical implementation of anisotropic hyperelastic material models including non-local damage". In: Biomechanics and Modeling in Mechanobiology 44.18–19 (2017), pp. 5894–1140. doi: [(10.1007/s10237-017- 0875-9)](https://link.springer.com/article/10.1007/s10237-017-0875-9).