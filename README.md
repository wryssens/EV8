# EV8
**v1 by P. Bonche, H. Flocard, P.-H. Heenen**   
**v2 by W. Ryssens, V. Hellemans, M. Bender and P.-H. Heenen**

The EV8 code solves the nuclear Skyrme-Hartree-Fock+BCS problem, representing the single-particle wavefunctions on a three-dimensional coordinate mesh. The [original](https://www.sciencedirect.com/science/article/abs/pii/S0010465505002821]) version (v1) as published in 2005, the [second version](https://www.sciencedirect.com/science/article/abs/pii/S0010465514003361)(v2) in 2015.

This github repository provides public access to a copy of the code. Contrary to the static repository at [Mendeley Data](https://data.mendeley.com/datasets/pbv7bz59rj/1), this repository can evolve as we find and fix errors in both the text and code.

Compared to the published version, this repository offers in addition

* fixes for small issues discovered by ourselves and other users,
* a list of known errors in the paper,
* this expanded readme file.

### Files 
---
Contained in this repository, you will find

* **README.md**: the file you are reading.
* **ev8_errors.pdf**: List of known errors (and their corrections) in the published paper of EV8v2.
* **Codes/**
   * **ev8.f**:  the EV8 source code.
   * **nil8.f**: Auxiliary program to generate starting wavefunctions for EV8.
   * **int8.f**: Auxiliary program to interpolate wavefunctions to different box parameters.
   * **den8.f**:  Auxiliary program to extract the nuclear density to a format that is readable by humans and non-FORTRAN programs. 

* **Examples/**
	* **Pb208.sh**: example demonstrating the useage of NIL8 and EV8 for the 208Pb for different sizes of the simulation box. Reference outpus of this script are
		* **Pb208.ev8.Sly4.SmallBox.out**
		* **Pb208.ev8.Sly4.BigBox.out**
	* **Zr84.sh**: example demonstrating NIL8 and EV8 for 84Zr, both without and with a constraint on the quadrupole deformation. Reference outpus of this script are
		* **zr84.out**
		* **z84.ev8.q200.out**
		
Warning: the example scripts assume the program codes are in the same directory as the scripts. Therefore, before running, either copy the program files from the "Codes" directory to the Examples directory or edit the script files to add the path, "../Codes/".

### Compilation
---
We do not provide a Makefile with the code. In order to compile an executable, you should
1. Create a param8.h datafile containing the numerical parameters of the mesh.
2. Compile ev8.f with your preferred compiler in the same directory of param8.h.

Please refer to the paper for more details on the contents of param8.h.

Caveats for the compilation process:
* The code is "fixed format FORTRAN". Some compilers requires you to specify a command-line option in order to correctly compile source code of this type.
* The code uses more continuation lines than allowed by the standard FORTRAN. The f95 compiler for example doesn't allow this by default. In that particular case the compiler option '-maxcontin 500' (or similar) might help you.
*  When dealing with large boxes or small stepsize, the code requires a large memory, which is not always handled well by all compilers. On the gfortran and portland compilers, for example, this results in severe compilation errors. You might want to add the '-mcmodel=medium' option to the compilation command in that case. When dealing with extremely large boxes or extremely small stepsizes, the previous option might not be enough. Again, for the gfortran and portland group compilers, you can probably use the '-mcmodel=large -fno-optimize-sibling-calls'.Note that the '-fno-optimize-sibling-calls' is not available for older versions of these compilers.
*  For recent versions of the gfortran compilers (starting with the 4.8 series) you will see a large number of compiler warnings. These are of no consequence however, and simply represent a new error-checking feature of gfortran. They are related to the way certain loops are coded, of which gfortran cannot internally check their validity.

### Running the code (in- and output)
----
Execution of the code is typically achieved by
```
./ev8.exe < ev8.data > output.out
```

where ev8.data is input for the code and output.out will contain all output pushed to STDOUT by the code. Furthermore, the code requires the presence of a file named (exactly) **fort.12** file and will produce a file named (exactly) **fort.13**. 

* **ev8.data**: Input datafile, containing **all** required input in the required format. Please refer to the paper and the examples for a detailed explanation.
* **fort.12**: An input wavefunction file that provides a starting point for the code. Can be produced with NIL8 or with previous runs of EV8. It should match the compilation parameters of your executables.
* **fort.13**: An output wavefunction file containing the endpoint of the iterative process. Can be used as a starting point for further EV8 runs. 


### Citing our work
-----

If you use EV8 in research work, please consider citing us. 

> P. Bonche, H. Flocard, P.-H. Heenen,  
> *Solution of the Skyrme-HF+BCS equation on a 3D mesh, II: A new version of the Ev8 code.*, Comp. Phys. Comm. 171, 49-62 (2005).
> [https://doi.org/10.1016/j.cpc.2005.05.001](https://doi.org/10.1016/j.cpc.2005.05.001)
>
> W. Ryssens, V. Hellemans, M. Bender and P.-H. Heenen,  
> *Solution of the Skyrme-HF+BCS equation on a 3D mesh, II: A new version of the Ev8 code*, Comp. Phys. Comm. 187, 175-194 (2015).
> [https://doi.org/10.1016/j.cpc.2014.10.001](https://doi.org/10.1016/j.cpc.2014.10.001)
>
>  W. Ryssens, V. Hellemans, M. Bender and P.-H. Heenen,
>  *Erratum: Solution of the Skyrme HF+BCS equation on a 3D mesh II. A new version of the Ev8 code (Computer Physics Communications (2015) 187:2 (175-194))*, Comp. Phys. Comm. 190, 231 (2015).
> [https://doi.org/10.1016/j.cpc.2015.01.011](https://doi.org/10.1016/j.cpc.2014.10.001)
  
