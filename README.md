# Track Analysis
This app analyses the track data of multi-particle tracking microrheology and calculate mean-square displacements (MSDs), the linear viscoelastic moduli (*G*' and *G*'') by the generalized Stokes-Einstein relation (GSER), and the non-Gaussian parameters, *etc*. 

Features:
* supports analysis upon multiple track files
* supports manual optimization for GSER calculation
* supports calculation of the distribution of diffusivity p(D) by Lucy's algorithm (see: Soft Matter 2018, 14, 3694-3703)

### Installation

Drag the .mlappinstall file to the command window of MATLAB. Then it will appear in the app list.

### Data preparation

The app can open any .mat file that contain a N x 4 matrix named tr.

