#ifndef MD_FUNC_H
#define MD_FUNC_H

void update_rank();
void update_forces();
void update_rank_forces(int mesh_node);
void get_kinetic_energy(int i);

void set_particles();
void get_distance_vector(int i, int j, double* dx, double* dy, double* dz);
int need_force(int pat_num, int mesh_node);
#endif
