# Second order metric perturbation of Kerr black holes

A Fortran ('08) and python code that solves Teukolsky equation for
the linearly perturbed Newman-Penrose scalar Psi4 about a Kerr black hole.
The code also directly reconstructs
the linear spacetime metric in outgoing raditation gauge from
the linearized Newman-Penrose scalar Psi\_4, and then
solves the equations of motion for the second order Psi\_4.
The code evolves fields in the time domain, and spatial derivatives
are evaluated using pseudo-spectral methods. 
For more information about the code and the formalism we use
see the papers listed under `Citation`.

Look under **Releases** for the latest stable version of this code.

Runtime parameters are configured in the `setup.py` file.

## Libraries

* mpmath: 
	http://mpmath.org/

* Fastest Fourier Transform in the West (FFTW): 
	http://fftw.org

* OpenMP (this is optional): 
	https://www.openmp.org/

I have successfully compiled the code with
gfortran (version 9) and ifort (version 17).

## Running on clusters
See `slurm` directory for some template slurm scripts 

## Derivation of equations of motion in coordinate form

A Mathematica notebook that contains the equations of motion
(as described in the `code paper` listed under `Citation`) in coordinate
form can be found [here](https://github.com/JLRipley314/2nd-order-teuk-derivations).

## Visualization

I use pyqtgraph-graph derived software
(see [here](https://github.com/JLRipley314/sci-vis))
to visualize the data, which are saved as csv files. 

## Citation

If you use this code please cite
our [paper](https://inspirehep.net/literature/1820630)
which describes and uses this code. Also consider citing our 
[paper](https://inspirehep.net/literature/1813628)
which describes our formalism in more detail.

## Contact

For questions please contact
Justin Ripley: lloydripley [at] gmail [dot] com
