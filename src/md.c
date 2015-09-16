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
    dt = 1e-12;
    H = 0.0;
    initTe = 10.0;
    initTp = 10.0;
    t_max = 500.0*200.0*10.0*dt;
    Lx = Ly = Lz = pow((double)Ne/ne,0.333333);

    assert(N % 2 == 0 );
    srand(time(NULL));

    set_particles();
    U = K = KelXY = KelZ = KpXY = KpZ = 0.0;
    for(i = 0; i < N; ++i) {get_kinetic_energy(i);}
    update_forces();
    E = K + U;

    initK = K;
    initU = U;
    initE = E;
    Imax = 0.0015 * fabs(initE) / t_max;

    FILE* f = fopen("energy", "w");
    double t;
    for(t = 0.0, i=0; t < t_max; t += dt, i++)
    {
        
        if(i % 1000 == 0)
            fprintf(f,"%.5g\t%.5g\t%.5g\t%.5g\n", t, KelXY, KelZ, E);
        
        U = K = KelXY = KelZ = KpXY = KpZ = 0.0;
        B_step(dt);
        E = U + K;

        if(initE != 0.0 && fabs(1.0 - E/initE) > 0.1)
        {
            fprintf(stderr,"!\n");
            FILE* d = fopen("dist", "w");
            double d_min = 1e6;
            double d_max = -1e6;
            double d_av  = 0.0;
            for(i = 0; i< N/2*(N-1); ++i)
            {
                if(dist[i] < d_min) {d_min = dist[i];}
                if(dist[i] > d_max) {d_max = dist[i];}
                d_av += dist[i];
            }
            d_av /= (N/2*(N-1));
            fprintf(d,"%.5g\t%.5g\t%.5g\t%.5g\n", t, d_av, d_min/d_av, d_max/d_av);
            return EXIT_FAILURE;
        }
    }
    
    return EXIT_SUCCESS;
}
