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
    dt = 1e-13;
    H = 5e4;
    t_max = 500.0*100.*dt;
    Lx = Ly = Lz = pow(50.0/1e10,0.333333);

    assert(N % 2 == 0 );
    srand(time(NULL));

    rx[0] = 0.5 * Lx + 100.0 * 100.0 * a0;
    vx[0] = 0.0;
    ry[0] = 0.5 * Ly;  
    vy[0] = sqrt(e*e / (me * 100.0 * 100.0 * a0));
    rz[0] = 0.5 * Lz;  
    vz[0] = 0.0;
    m[0] = me; q[0] = -e; type[0] = 1;

    rx[1] = 0.5 * Lx;
    vx[1] = 0.0;
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

    int st = 0;
    double t;
    double dx, dy, dz;
    FILE* f_en = fopen("energy","w");
    FILE* f_pos = fopen("part_pos", "w");

    get_distance_vector(0, 1, &dx, &dy, &dz);
    fprintf(f_pos,"%d\t%g\t%g\t%g\t%g\t%g\n",st, rx[0]/Lx, ry[0]/Ly, rx[1]/Lx, ry[1]/Ly, sqrt(dx*dx + dy*dy + dz*dz));
    fprintf(f_en,"%d\t%g\t%g\t%g\n", st, E, initE, (1.0 - fabs(E/initE)));
    for(t = 0.0; t < t_max; t += dt)
    {
        U = K = KelXY = KelZ = KpXY = KpZ = 0.0;
        /*VV_step(dt);*/
        B_step(dt);
        E = U + K;
        st++;
        
        get_distance_vector(0, 1, &dx, &dy, &dz);
        fprintf(f_pos,"%d\t%g\t%g\t%g\t%g\t%g\n",st, rx[0]/Lx, ry[0]/Ly, rx[1]/Lx, ry[1]/Ly, sqrt(dx*dx + dy*dy + dz*dz));
        fprintf(f_en,"%d\t%g\t%g\t%g\n", st, E, initE, (1.0 - fabs(E/initE)));
        if(initE != 0.0 && fabs(initE - E) > 0.1 * fabs(initE))
        {
            fprintf(stderr,"!\n");
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}
