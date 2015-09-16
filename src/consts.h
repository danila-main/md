#ifndef MD_CONSTS_H 
#define MD_CONSTS_H 

#define Ne 1000
#define Np 1000
#define N (Ne + Np)

extern double Lx, Ly, Lz;       /*cell size in each dimension*/ 
extern double ne;               /*initial free electrons density [cm^3]*/
extern double dt;               /*time step [s]*/
extern double H;                /*magnetic field induction [Gauss]*/ 
extern double initTe;           /*initial temperature of electrons [K]*/
extern double initTp;           /*initial temperature of proton [K]*/
extern double t_max;            /*max time [s]*/ 

extern double initK;
extern double initU;
extern double initE;

extern double K;
extern double KelXY;
extern double KelZ;
extern double KpXY;
extern double KpZ;
extern double U;
extern double E;

/****  Fundamental physical constants  ******/                                                                          
extern const double Pi;                                                                                                 
extern const double Kb;         /* Boltzamann's constant  */
extern const double me;         /* mass of electron */
extern const double mp;         /* mass of proton*/
extern const double a0;         /* Bohr's radius*/
extern const double e;          /* electron's charge */
extern const double c;          /* speed of light*/

extern double rx[];
extern double ry[];
extern double rz[];

extern double vx[];
extern double vy[];
extern double vz[];

extern double Fx[];
extern double Fy[];
extern double Fz[];

/*auxiliary arrays*/
extern double dist[];
extern double q[];
extern double m[];
extern int type[];

extern double rank_tmp[N][N];
extern int rank[];
extern int rank_max;
extern double Imax;

#define ELECTRON 1
#define PROTON 2
#define ALL 0



#endif
