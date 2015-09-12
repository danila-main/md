#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "consts.h"
#include "func.h"

double sgn(double x);
double rnd();
double normal_rnd();

void update_forces()
{
    int i, j, l;
    double dx, dy, dz;
    double dr, dr2;
    double k, f, fx, fy, fz;
    U = 0.0;
    for(i = 0; i < N; ++i){
        Fx[i] = Fy[i] = Fz[i] = 0.0;}

    if(1)
    {
        l = 0;
        for(i = 0; i < N-1; ++i)
        {
            for(j = i+1; j < N; ++j)
            {
                get_distance_vector(i, j, &dx, &dy, &dz);
                dr = sqrt(dx*dx + dy*dy + dz*dz);
                dist[l] = dr;
                k = (-1.0) * q[i] * q[j];
                f = k * (1.0/(dr * dr));
                fx = f * dx/dr;
                fy = f * dy/dr;
                fz = f * dz/dr;

                Fx[i] += fx; Fx[j] -= fx;
                Fy[i] += fy; Fy[j] -= fy;
                Fz[i] += fz; Fz[j] -= fz;

                U += (-k*(1.0/dr));
                ++l;
            }
        }
    }
}

void set_particles()
{
    int i;
    double sigma;
    sigma = sqrt(Kb*initTe/me);
    for(i = 0; i < Ne; ++i) 
    {  
        type[i] = ELECTRON; 
           q[i] = -e;    
           m[i] = me;   
          rx[i] = rnd() * Lx; 
          ry[i] = rnd() * Ly;
          rz[i] = rnd() * Lz;
          vx[i] = normal_rnd() * sigma; 
          vy[i] = normal_rnd() * sigma;
          vz[i] = normal_rnd() * sigma;
    } 

    sigma = sqrt(Kb*initTp/mp);
    for(i = Ne; i < Ne+Np; ++i) 
    {  
        type[i] = PROTON; 
           q[i] = e;    
           m[i] = mp;   
          rx[i] = rnd() * Lx; 
          ry[i] = rnd() * Ly;
          rz[i] = rnd() * Lz;
          vx[i] = normal_rnd() * sigma; 
          vy[i] = normal_rnd() * sigma;
          vz[i] = normal_rnd() * sigma;
    }
}

void get_kinetic_energy(int i)
{
    if(type[i] == ELECTRON)                                                                                        
    {                                                                                                                   
        KelXY += (0.5 * (vx[i] * vx[i]) * m[i] +                                                     
                  0.5 * (vy[i] * vy[i]) * m[i]);                                                      
        KelZ +=   0.5 * (vz[i] * vz[i]) * m[i];                                                      
    }                                                                                                                   
    else if (type[i] == PROTON)                                                                                    
    {                                                                                                                   
        KpXY += (0.5 * (vx[i] * vx[i]) * m[i] +                                                     
                 0.5 * (vy[i] * vy[i]) * m[i]);                                                      
        KpZ +=   0.5 * (vz[i] * vz[i]) * m[i];                                                      
    }
    K += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]) * m[i];
}

void get_distance_vector(int i, int j, double* dx, double* dy, double* dz)
{
    *dx = rx[j] - rx[i];
    *dy = ry[j] - ry[i];
    *dz = rz[j] - rz[i];

    if(fabs(*dx) > 0.5 * Lx)
        *dx -= sgn(*dx) * Lx;
    if(fabs(*dy) > 0.5 * Ly)
        *dy -= sgn(*dy) * Ly;
    if(fabs(*dz) > 0.5 * Lz)
        *dz -= sgn(*dz) * Lz;
}

double sgn(double x)
{
    return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}

double rnd()
{
    return (((double)rand())/((double)RAND_MAX));
}

double normal_rnd()
{
    int i;
    int iter = 12;
    double sum = 0.0;
    for(i = 0; i < iter; i++){
        sum += rnd();}
    return sqrt(12.0/((double)iter)) * (sum - 0.5*iter);
}
