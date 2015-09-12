#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>

#include "consts.h"
#include "func.h"
#include "step.h"

int main()
{
    int i;
    ne = 1e9;
    dt = 1e-13;
    H = 0.0;
    initTe = 10.0;
    initTp = 10.0;
    t_max = 500.0*200.0*dt;
    Lx = Ly = Lz = pow((double)Ne/ne,0.333333);

    assert(N % 2 == 0 );
    srand(time(NULL));

    set_particles();

    FILE* f_exp = fopen("exp_dist","w");
    FILE* f_th = fopen("th_dist", "w");
    for(i = 0; i < Ne; i++)
    {
        fprintf(f_exp,"%g\t%g\t%g\t%g\n",vx[i], vy[i], vz[i],
                                         sqrt(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]));
    }
    double Vx, V;
    double Nx, Ntot;
    for(i = 0; i < Ne; i++)
    {
        Vx = vx[i];
        Nx = Ne*sqrt(me/(2.0*Pi*Kb*initTe))*exp(-0.5*me*Vx*Vx/(Kb*initTe))*0.1e6;
        V = sqrt(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
        Ntot = 4.0*Pi*Ne*V*V*pow((me/(2.0*Pi*Kb*initTe)),1.5)*exp(-0.5*me*V*V/(Kb*initTp))*0.1e6;
        fprintf(f_th,"%g\t%g\t%g\t%g\n",Vx,Nx,V,Ntot);
    }

    return EXIT_SUCCESS;
}
