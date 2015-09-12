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
    dt = 1e-11;
    H = 5e4;
    t_max = 500.0*10*dt;
    Lx = Ly = Lz = pow(50.0/1e10,0.333333);
    double T = 2.0 * Pi * me * c/ (e * H);
    fprintf(stderr,"T = %g\t dt = %g\n", T, dt/T);

    assert(N % 2 == 0 );
    srand(time(NULL));

    rx[0] = 0.333 * Lx + 0.01*Lx;
    vx[0] = -1.0*sqrt(10.0 * Kb / me);
    ry[0] = 0.5 * Ly;  
    vy[0] = 0.0;
    rz[0] = 0.5 * Lz;  
    vz[0] = 0.0;
    m[0] = me; q[0] = -e; type[0] = 1;

    rx[1] = 0.666 * Lx - 0.01*Lx;
    vx[1] = sqrt(10.0 * Kb / mp);
    ry[1] = 0.5 * Ly;  
    vy[1] = 0.0;
    rz[1] = 0.5 * Lz;  
    vz[1] = 0.0;
    m[1] = mp; q[1] = e; type[1] = 2;

    U = K = KelXY = KelZ = KpXY = KpZ = 0.0;
    for(i = 0; i < N; ++i) {get_kinetic_energy(i);}
    update_forces();
    E = K + U;

    initK = K;
    initU = U;
    initE = E;

    int st;
    double t;
    double dx, dy, dz;
    FILE* f_en = fopen("energy","w");
    FILE* f_pos = fopen("part_pos", "w");
    st = 0;
    fprintf(f_pos,"%d\t%g\t%g\t%g\t%g\n",st, rx[0]/Lx, ry[0]/Ly, rx[1]/Lx, ry[1]/Ly);
    fprintf(f_en,"%d\t%g\t%g\t%g\n", st, E, initE, (1.0 - fabs(E/initE)));
    //FILE* f_pos = fopen("part_pos", "w");
    for(t = 0.0, st = 0; t < t_max; t += dt)
    {
        U = K = KelXY = KelZ = KpXY = KpZ = 0.0;
        //VV_step(dt);
        B_step(dt);
        E = U + K;
        st++;

        fprintf(f_pos,"%d\t%g\t%g\t%g\t%g\n",st, rx[0]/Lx, ry[0]/Ly, rx[1]/Lx, ry[1]/Ly);
        fprintf(f_en,"%d\t%g\t%g\t%g\n", st, E, initE, (1.0 - fabs(E/initE)));
        if(initE != 0.0 && fabs(1.0 - E/initE) > 0.1)
        {
            fprintf(stderr,"!\n");
            return EXIT_FAILURE;
        }
    }
    
    return EXIT_SUCCESS;
}
