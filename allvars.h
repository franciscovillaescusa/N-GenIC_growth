#include <srfftw_mpi.h>

#define  PI          3.14159265358979323846 
#define  GRAVITY     6.672e-8
#define  HUBBLE      3.2407789e-18   /* in h/sec */
#define  C           2.9979e10


#define  T_CMB0      2.725	/* present-day CMB temperature */

/* Note there is a slight correction from 4/11
 * due to the neutrinos being slightly coupled at e+- annihilation.
 * See Mangano et al 2005 (hep-ph/0506164)
 * The correction is (3.046/3)^(1/4), for N_eff = 3.046 */
#define TNU     (T_CMB0*pow(4/11.,1/3.)*1.00381)              /* Neutrino + antineutrino background temperature in Kelvin */

/*With slightly relativistic massive neutrinos, for consistency we need to include radiation.
 * A note on normalisation (as of 08/02/2012):
 * CAMB appears to set Omega_Lambda + Omega_Matter+Omega_K = 1,
 * calculating Omega_K in the code and specifying Omega_Lambda and Omega_Matter in the paramfile.
 * This means that Omega_tot = 1+ Omega_r + Omega_g, effectively
 * making h0 slightly larger than specified.
 * CLASS seems to set Omega_tot = 1, calculating Omega_Lambda (by default) and setting everything
 * else in the paramfile.
 * CLASS is clearly the saner option, but for the moment we will follow CAMB.*/
/*Stefan-Boltzmann constant in cgs units*/
#define STEFAN_BOLTZMANN 5.670373e-5
/* Omega_g = 4 \sigma_B T_{CMB}^4 8 \pi G / (3 c^3 H^2)*/
#define OMEGAG (4*STEFAN_BOLTZMANN*8*M_PI*GRAVITY/(3*C*C*C*HUBBLE*HUBBLE)*pow(T_CMB0,4)/HubbleParam/HubbleParam)

/*Neutrinos are included in the radiation*/
/*For massless neutrinos, rho_nu/rho_g = 7/8 (T_nu/T_cmb)^4 *N_eff, but we absorbed N_eff into T_nu above*/
#define OMEGANU (OMEGAG*7/8.*pow(TNU/T_CMB0,4)*3)
/*With massless neutrinos only, add the neutrinos to the radiation*/

/* For convenience define OMEGAK so we can easily switch between CLASS and CAMB normalisations*/
#define OMEGAK (1-Omega - OmegaLambda)
/*CLASS normalisation; needs a slightly lower OmegaLambda*/
//#define OMEGAK (1-All.Omega0 - All.OmegaLambda - OMEGAR)




double PowerSpec(double kmag);
double GrowthFactor(double astart, double aend);
double F_Omega(double a);
int    read_parameter_file(char *fname);
double PowerSpec_EH(double k);
double PowerSpec_Efstathiou(double k);


#ifdef T3E
typedef short int int4byte;	/* Note: int has 8 Bytes on the T3E ! */
typedef unsigned short int uint4byte;	/* Note: int has 8 Bytes on the T3E ! */
#else
typedef int int4byte;
typedef unsigned int uint4byte;
#endif



extern struct io_header_1
{
  uint4byte npart[6];      /*!< npart[1] gives the number of particles in the present file, other particle types are ignored */
  double mass[6];          /*!< mass[1] gives the particle mass */
  double time;             /*!< time (=cosmological scale factor) of snapshot */
  double redshift;         /*!< redshift of snapshot */
  int4byte flag_sfr;       /*!< flags whether star formation is used (not available in L-Gadget2) */
  int4byte flag_feedback;  /*!< flags whether feedback from star formation is included */
  uint4byte npartTotal[6]; /*!< npart[1] gives the total number of particles in the run. If this number exceeds 2^32, the npartTotal[2] stores
                                the result of a division of the particle number by 2^32, while npartTotal[1] holds the remainder. */
  int4byte flag_cooling;   /*!< flags whether radiative cooling is included */
  int4byte num_files;      /*!< determines the number of files that are used for a snapshot */
  double BoxSize;          /*!< Simulation box size (in code units) */
  double Omega0;           /*!< matter density */
  double OmegaLambda;      /*!< vacuum energy density */
  double HubbleParam;      /*!< little 'h' */
  int4byte flag_stellarage;     /*!< flags whether the age of newly formed stars is recorded and saved */
  int4byte flag_metals;         /*!< flags whether metal enrichment is included */
  int4byte hashtabsize;         /*!< gives the size of the hashtable belonging to this snapshot file */
  char fill[84];		/*!< fills to 256 Bytes */
}
header, header1;


extern int      Nglass;
extern int      *Local_nx_table;
extern int      WhichSpectrum;


extern FILE     *FdTmp, *FdTmpInput;

extern int      Nmesh, Nsample;

extern int      SphereMode;

extern long long IDStart;


extern char     GlassFile[500]; 
extern char     FileWithInputSpectrum[500];
extern char     FileWithTransfer[500];
extern char     FileWithBGrowth[500];
extern char     FileWithCDMGrowth[500];
extern char     FileWithNUGrowth[500];
extern int      GlassTileFac; 

extern double   Box;
extern int RayleighSampling;
extern int Phase_flip;
extern int Seed;


extern long long TotNumPart;

extern int      NumPart;

extern int      NTaskWithN;


extern int      *Slab_to_task;


extern struct part_data 
{
  float Pos[3];
  float displ[3]; //Particle positions are updated only at the end
                  //displ keeps the displacement field
  float Vel[3];
#ifdef  MULTICOMPONENTGLASSFILE                      
  int   Type;
#endif
  long long ID;
} *P;


extern double InitTime;
extern double Redshift;
extern double MassTable[6];


extern char OutputDir[100], FileBase[100];
extern int  NumFilesWrittenInParallel;


extern int      ThisTask, NTask;

extern int      Local_nx, Local_x_start;

extern int  IdStart;

extern rfftwnd_mpi_plan Inverse_plan;
extern fftw_real        *Disp, *Workspace, *Disp2, *Workspace2;
extern fftw_complex     *Cdata, *Cdata2;


extern double UnitTime_in_s, UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
extern double InputSpectrum_UnitLength_in_cm;
extern double G, Hubble;
extern double RhoCrit;

extern double Omega, OmegaLambda, OmegaDM_2ndSpecies, Sigma8;
extern double OmegaBaryon, HubbleParam,Hzi;
extern double PrimordialIndex;
extern double ShapeGamma;

extern double Dplus; /* growth factor */


#ifdef DIFFERENT_TRANSFER_FUNC
extern int Type, MinType, MaxType;
#endif

extern int    ReNormalizeInputSpectrum;

extern int    WDM_On;
extern int    WDM_Vtherm_On;
extern double WDM_PartMass_in_kev;

extern int    NU_On;
extern int    NU_Vtherm_On;
extern double NU_PartMass_in_ev;

#ifdef OUTPUT_DF
// define the pointers for the modes coordinates, amplitudes and phases
extern long long *coord;
extern float     *amplitudes, *phases;
#endif
