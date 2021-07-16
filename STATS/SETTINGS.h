//DEFINE ALL THE PARAMETERS NEEDED TO CALCULATE THE STATISTICS FOR THE GIOMETTO CASE

//GENERAL SETTINGS

//------------------------------------------------------
//Clear: delete all the environment, INCLUDING the STATISTICS OUTPUT. This is the very first operation of this program, and it can be followed by the "First time run" option, if selected.
bool clear_bool=true;
//First time run: create the environment (copy data files, create the subdirectories etc). This should be done for a first-time run or if the previous "Clear" option has been selected.
bool environment_bool=false;
//------------------------------------------------------

//SIMULATION SETTINGS: geometrical and physical parameters defining the current simulation

//------------------------------------------------------
//NOTE: all the following parameters are intended to be dimensionful, the MKS standard unit system is adopted.
//Grid stencil: number of cells, cell width in x and y direction
int Nx=154;
int Ny=154;
int Nz=410;
double Lx=0.0095183;
double Ly=0.0095183;
//Geometrical and dynamical parameters
double alpha=3.14159/6.0;
double N=std::sqrt(0.164414138);
double Pr=1.0;
double nu=1.5e-5;
double Gr=1.0/(1.5e-5*1.5e-5*pow(N, 6));
//Parameters to normalize data: if no normalization is needed, set to 1.0
double b_norm=1.0;
double U_norm=1.0/N;
double T_norm=1.0/N;
double Z_norm=1.0/(N*N);
//Indexes of the 5 selected heights. These are the heights where the autocorrelation profiles will be extracted from for visualization. The most suitable heights should be chosen according to the current simulation (check the Umean profiles first, for example) and the associated vertical index should be found in the z.dat file. For later convenience, it is very useful to indicate in a comment on the right, what are the real heights, given in dimensionless units. The following values are just an example: z_i= 1.1e-3; 1.0e-2; 0.051; 0.1; 0.151. 
int k1=10, k2=67, k3=182, k4=248, k5=296; //respectively, 1.1e-3; 1.0e-2; 0.051; 0.1; 0.151; ( <= dimless)
//------------------------------------------------------

//SETTINGS FOR THE STATISTICS UTILITIES

//------------------------------------------------------
//NOTE: in case of second order statistics or autocorrelations, please check if you have the first order statistics, otherwise calculate them.
// 1) First order statistics: mean vertical profiles of velocity (U), buoyancy (b) and nu_t/nu ratio.
bool firstOrder_bool=false;
// 2) Second order statistics: rms of b and U, covariances UxUz and bU, turbulent kinetic balance.
bool secondOrder_bool=false;
// 3) Autocorrelation profiles at selected heights: along x and y directions, by averaging over the time and the other horizontal direction.
bool autocorrelation_bool=false;
//------------------------------------------------------

