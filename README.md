# Track Analysis
### Introduction

This app analyses the track data of multi-particle tracking microrheology and calculates the mean-square displacements (MSDs), the linear viscoelastic moduli (*G*&prime; and *G*&Prime;) by the generalized Stokes--Einstein relation (GSER), and the non-Gaussian parameters, *etc*. [^*]

Features:
* supports analysis upon multiple track files
* supports manual optimization for GSER calculation
* supports calculation of the distribution of diffusivity *p*(*D*) by Lucy's algorithm[^†][^‡]

### Installation

Unzip and drag the `.mlappinstall` file to the command window of MATLAB. After the installation it will appear in the app list. 

### Data preparation

Consider a time-lapse video of *M* frames. The resolution of the video is *X* pixels &times; *Y* pixels. Each frame is labeled by a positive integer, called "frame ID". Suppose that previous image analysis has provided a result of a number of tracks from this video. Each track is labeled by a positive integer, called "track ID". The information of a track (bearing the same track ID) is provided by the data of *x* and *y* positions (in the unit of pixels), each of which corresponds to a distinct frame_ID's. So, each row of a track data consists of the ordered quadruple (x, y, frame ID, track ID). In general, different tracks would start at different frame ID's and last for different number of frames. Also, tracks that skipping some of the frames (i.e. under the same track ID, the frame ID's are not a set of continuous natural numbers) can be handled by this app. An *N* &times; 4 matrix each row of which is the quadruple (x, y, frame ID, track ID) is sufficient to completely describe a set of multiple tracks.

The app can open any `.mat` file that contains an *N* &times; 4 matrix named `tr`. The 4 columns of `tr` are *x* positions, *y* positions, frame ID, and tracks ID in order. See `tr_example.mat` for an example.

Another app, [Particle Tracking](https://github.com/andrewx101/particle_tracking), is designed for generating the `tr` data from time-lapse series.

Within the same `tr` data, the ordered pairs (frame ID, track ID) must be pair-wise distinct, since no one particle (the same track ID) can appear at two or more positions in the same moment (frame ID).

The app can process multiple '.mat' files. For example, they may be the results from different videos but on the same sample. In this case the track ID and frame ID from the `tr` of different imported `.mat` files do not need to be distinct video-wise. In other words, in each of these `tr` variables, the track ID's and frame ID's can, without exception, be (re-)starting from 1, regardless of the others.

### Load data

Click `Add...` to select `.mat` files that contain the matrix `tr`. Multiple selection is supported. Click `Add...` many times to append more tracks. The files are listed in the table next to the button. Click `Remove` to remove the currently selected file. The number of tracks contains in each file is also listed in this table.

### Setting parameters

The parameters in the `Parameters` panel should be set before calculation. Their meanings are listed below.

#### Resolution (&micro;m/pixel)

The conversion factor of the unit of *x*/*y* position from pixel to micron. This value is determined by independent calibration of the microscope with a digital camera.

#### Temperature (&deg;C)

The temperature of the experiment.

#### Particle radius (&mu;m)

The radius of the probe particles.

#### Shutter time (s)

The lag time between adjacent steps of each tracks. This should be the same for all loaded tracks.

### Calculate

Click `Calculate` to start the calculation upon the currently loaded tracks. Status that indicating the currently processing quantity will be temporarily shown below the axes. Before finishing the calculation a window will pop up for the user to adjust the parameter for GSER calculation (will be introduced in later sections).  After closing this window the calculation finishes. The status label read `Done!` and the mean-square displacement (MSD) is plotted on the axes. 

### The GSER pop-up window

During the first calculation or after clicking the `GSER` button the GSER pop-up window is shown for tuning. GSER calculation can be optimized by a) cutting unwanted portion of the MSD and b) smoothing the MSD curve and its derivative successively by the same smoothing algorithm and setting. By default the MSD is uncut; it and its derivative are smoothed by the moving-average algorithm. In the left figure, the original MSD curve is plotted as blue dots, and the cut and smoothed MSD curve is plotted as a red solid line. In the right figure the resulted *G*' and *G*'' are plotted. Change the two slider the set the minimum and maximum &Delta;*t* to cut out unwanted portion of the MSD. The resulted upper and lower bounds are indicated as black dashed line in the left figure. Change smoothing algorithm by selecting the corresponding radio button in the `Smooth algorithm` panel. Change the `Span` (for moving-average and robust lowess algorithms), and `Knots` (for SLM engine) by the slider next to the panel. After any change both figures will update. Tune until the results are good. Click `OK` or directly close the GSER pop-up windows will quit this session and confirm the GSER results.

### Select plotting data

The `Plot data` panel allow selection of the following observables to plot.

#### MSD

The multi-particle mean-square displacements <&Delta;*x*<sup>2</sup>>, <&Delta;*y*<sup>2</sup>>, and <&Delta;*r*<sup>2</sup>>, where &Delta;*r*<sup>2</sup> = &Delta;*x*<sup>2</sup> + &Delta;*y*<sup>2</sup>. A black dash line indicates the slope of unity.

#### Modulus

The storage and loss moduli *G*&prime; and *G*&Prime; calculated by the `GSER` pop-up window.

#### alpha_2

The mean single-particle and multi-particle non-Gaussian parameters[^†]. The mean single-particle parameter is plotted with its standard error among all tracks.

#### van Hove

The distribution density of &Delta;*x* (van Hove function, scattered) and the corresponding Gaussian distribution (solid) estimated from the MSD for the lag time &Delta;*t* adjusted by the `Lag time` slider.

#### SP-MSD

Single-particle MSD <&Delta;*x*<sup>2</sup>><sub>single</sub> (solid) and the corresponding multi-particle counterpart (dash).

### Lucy

For each lag time &Delta;*t*, adjusted by the `Lag time` slider, the van Hove function can be deconvolve into a distribution of diffusivity *P*(*D*) using Lucy's algorithm[^‡]. First select the lag time under which you want to calculate *P*(*D*), then click the `Lucy` button to start the calculation. After the calculation the result will replace the content in the figure. *P*(*D*) is ploted as blue patches, and the corresponding cumulated distribution is plotted as line. Each click on the `Lucy` button will replace the old result of *P*(*D*).

### Saving data

By clicking the `Save...` button a folder selection dialogue popup and the current calculated results will be saved into a MATLAB data file named `msd_data.mat` in the selected folder.  If a Lucy calculation has performed an addition MATLAB data file named `PD_data_xxx.mat` in the same folder, where xxx is the lag time under which the *P*(*D*) was calculated. If no Lucy calculation has performed, a dialog stating this will appear after saving.

The `msd_data.mat` contains the following 6 vairables.

#### `msd`

*n* &times; 4 matrix of the multi-particle MSDs, where *n* is the number of lag times. The 4 columns are: &Delta;*t* (s),  <&Delta;*x*<sup>2</sup>> (&micro;m<sup>2</sup>), <&Delta;*y*<sup>2</sup>> (&micro;m<sup>2</sup>), and <&Delta;*r*<sup>2</sup>> (&micro;m<sup>2</sup>), respectively. 

#### `msd_i`

*m* &times; 1 cell array of the single-particle MSD, where *m* is the total number of tracks processed. Each cell is an *n* &times 4 matrix with the same column meaning as `msd`.

#### `G`

3-column matrix of the GSER result. The columns are angular frequency *&omega;* (rad/s), *G*&prime; (Pa), *G*&Prime; (Pa), respectively.

#### `alpha_2` and `alpha_2_i`

Multi-particle and mean single-particle non-Gaussian parameters, both of which are 3-column matrices. The 3 columns of `alpha_2_i` are &Delta;*t* (s), the mean *&alpha;*<sub>2</sub> among multiple tracks, and the standard error of *&alpha;*<sub>2</sub> among multiple tracks. The 3 columns of `alpha_2` are &Delta;*t* (s), the multi-particle *&alpha;*<sub>2</sub>, and the finite-sample variance of the multi-particle *&alpha;*<sub>2</sub>.

#### `van_Hove`

The displacement &Delta;*x* distribution density *p*(&Delta;*x*), an *n* &times 2 cell array.  The 2 cells in each row are the distribution densities and the bins edges, respectively.

The `PD_data_xxx.mat` contains the following 3 variables

#### `p_D` and `P_D`

The distribution density and distribution function of diffusivity, respectively. Both are 2-column matrices: *D* (m<sup>2</sup>/s), *p*(*D*) (s m<sup>&#8722;2</sup>) or *P*(*D*) (-).

#### `van_Hove_rec`

A reconstructed van_Hove function from the calculated *P*(*D*) data. It is a 2-column matrix: &Delta;*x* (micro;m), distribution density *p*(&Delta;*x*) (You may check the validity of the latter by comparing `van_Hove_rec` with the van_Hove data of the same lag time in the `msd_data.mat` file.

[![View track_analysis on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://ww2.mathworks.cn/matlabcentral/fileexchange/118310-track_analysis)

[^*]: This `readme.md` and also the app are for users familiar with the field of probe microrheology. See for example: Waigh, T. A. "Microrheology of Complex Fluids." *Reports on Progress in Physics* 68, no. 3 (2005): 685-742.
[^†]: See for example: Hong, W., *et al.* &ldquo;Colloidal Probe Dynamics in Gelatin Solution During the Sol–Gel Transition.&rdquo; *Soft Matter* 14, no. 19 (2018): 3694-703. 
[^‡]: Lucy, L. B. &ldquo;An Iterative Technique for the Rectification of Observed Distributions.&rdquo; *The Astronomical Journal* 79, no. 6 (1974): 745-54. 
