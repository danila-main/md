#include "consts.h"

double Lx = 0.0;
double Ly = 0.0;
double Lz = 0.0;
double ne = 0.0;
double dt = 0.0;
double H  = 0.0;
double initTe = 0.0;
double initTp = 0.0;
double t_max  = 0.0;

double initK = 0.0;
double initU = 0.0;
double initE = 0.0;

double K = 0.0;
double KelXY = 0.0;
double KelZ = 0.0;
double KpXY = 0.0;
double KpZ = 0.0;
double U = 0.0;
double E = 0.0;

const double Pi = 3.14159265;                                                                                           
const double Kb = 1.38e-16;         /* Boltzamann's constant  */ 
const double me = 9.1e-28;          /* mass of electron */ 
const double mp = 1.67e-24;         /* mass of proton*/ 
const double a0 = 5.29e-9;          /* Bohr's radius*/ 
const double e = 4.8e-10;           /* electron's charge */ 
const double c = 29979245800.0;     /* light speed */           

double rx[N];
double ry[N];
double rz[N];

double vx[N];
double vy[N];
double vz[N];

double Fx[N];
double Fy[N];
double Fz[N];

double dist[N/2*(N-1)];
double q[N];
double m[N];
int type[N];

double rank_tmp[N][N];
int rank[N];
int rank_max = 0;
double Imax = 0.0;
