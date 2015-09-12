#ifndef MD_FUNC_H
#define MD_FUNC_H

void update_forces();
void get_kinetic_energy(int i);

void set_particles();
void get_distance_vector(int i, int j, double* dx, double* dy, double* dz);

#endif
