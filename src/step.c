#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "consts.h"
#include "step.h"
#include "func.h"

/*helpers*/
double Ax, Ay, Az;
double Mx, My, Mz;
double Vx, Vy, Vz;
double Ex, Ey, Ez;
double avEx, avEy, avEz;
double pEx[N];
double pEy[N];
double pEz[N];
double omega;
double I, II, III;

void B_update_position(int i, double dt);
void B_update_velocity(int i, double dt);

void B_step(double dt)
{
    int i, m;
    int steps_num;
    U = K = KelXY = KelZ = KpXY = KpZ = 0.0;

    update_rank();
    if (rank_max > 10)
        fprintf(stderr, "%d\n", rank_max);
    steps_num = (int)pow(2.0,(double)rank_max);
    double tau = dt / steps_num;
    for(m = 0; m < steps_num; ++m)
    {
        for(i = 0; i < N; ++i)
        {
            B_update_position(i, tau);
        }

        if((m+1) != steps_num)
            update_rank_forces(m+1);
        else
        {
            U = K = KelXY = KelZ = KpXY = KpZ = 0.0;
            update_forces();
        }

        for(i = 0; i < N; ++i)
        {
            if ( need_force(i, m+1) )
            {
                avEx = 0.5 * (Fx[i]/q[i] + pEx[i]); 
                avEy = 0.5 * (Fy[i]/q[i] + pEy[i]); 
                avEz = 0.5 * (Fz[i]/q[i] + pEz[i]);
            }
            else
            {
                assert(pEx[i] == Fx[i]/q[i]);
                avEx = pEx[i]; 
                avEy = pEy[i]; 
                avEz = pEz[i];
            }
            B_update_velocity(i, tau);
            get_kinetic_energy(i);
        }
    }
	return;
}

void B_update_position(int i, double dt)
{
    omega = (H * q[i])/(c * m[i]);
    Ex = Fx[i]/q[i]; Ey = Fy[i]/q[i]; Ez = Fz[i]/q[i]; 

    /*now save E*/
    pEx[i] = Ex; pEy[i] = Ey; pEz[i] = Ez;

    /*Calculate position at dt*/
    if(H != 0.0)
    {
        I = vx[i] * (sin(omega*dt) / omega);
        II = (q[i]/m[i])*Ey * ((omega*dt - sin(omega*dt)) / (omega*omega));
        III = (vy[i]*omega + (q[i]/m[i])*Ex) * ((1.0-cos(omega*dt))/(omega*omega));
    }
    else
    {
        I = vx[i]*dt;
        II = 0.0;
        III = (vy[i]*omega + (q[i]/m[i])*Ex) * (0.5*dt*dt);
    }
    rx[i] += I + II + III;

    if(H != 0.0)
    {
        I = vy[i] * (sin(omega*dt) / omega);
        II = -1.0 * (q[i]/m[i])*Ex * ((omega*dt - sin(omega*dt)) / (omega*omega));
        III = -1.0 * (vx[i]*omega - (q[i]/m[i])*Ey) * ((1.0-cos(omega*dt))/(omega*omega));
    }
    else
    {
        I = vy[i]*dt;
        II = 0.0;
        III = -1.0 * (vx[i]*omega - (q[i]/m[i])*Ey) * (0.5*dt*dt);
    }
    ry[i] += I + II + III;

    rz[i] += vz[i]*dt + (q[i]/m[i])*Ez * 0.5*dt*dt;

    if( rx[i] > Lx) rx[i] -= Lx;
    if( rx[i] < 0)  rx[i] += Lx;

    if( ry[i] > Ly) ry[i] -= Ly;
    if( ry[i] < 0)  ry[i] += Ly;

    if( rz[i] > Lz) rz[i] -= Lz;
    if( rz[i] < 0)  rz[i] += Lz;

    return;
}

void B_update_velocity(int i, double dt)
{
    omega = (H * q[i])/(c * m[i]);

    Vx = vx[i]; Vy = vy[i]; Vz = vz[i];

    if(H != 0.0)
    {
        I = Vx * cos(omega*dt);
        II = q[i]/(m[i]*omega)*avEy*(1 - cos(omega*dt));
        III = (Vy*omega + (q[i]/m[i])*avEx) * sin(omega*dt)/omega;
    }
    else
    {
        I = Vx;
        II = 0.0;
        III = (Vy*omega + (q[i]/m[i])*avEx) * dt;
    }
    vx[i] =  I + II + III;

    if(H != 0.0)
    {
        I = Vy * cos(omega*dt);
        II = -1.0 * q[i]/(m[i]*omega)*avEx*(1 - cos(omega*dt));
        III = -1.0 * (Vx*omega - (q[i]/m[i])*avEy) * sin(omega*dt)/omega;
    }
    else
    {
        I = Vy;
        II = 0.0;
        III = -1.0 * (Vx*omega - (q[i]/m[i])*avEy) * dt;
    }
    vy[i] = I + II + III;

    vz[i] = Vz + (q[i]/m[i])*avEz*dt;
    return;
}






















































void VV_update_position(int i, double dt);
void VV_update_velocity(int i, double dt);

void VV_step(double dt)
{
    int i;
    U = K = KelXY = KelZ = KpXY = KpZ = 0.0;
    for(i = 0; i < N; i++)
    {
        VV_update_velocity(i, 0.5*dt);
        VV_update_position(i, dt);
    }

    update_forces();

    for(i = 0; i < N; i++)
    {
        VV_update_velocity(i, 0.5*dt);
        get_kinetic_energy(i);
    }
}

void VV_update_velocity(int i, double dt)
{
    omega = (H * q[i])/(c * m[i]);
    Mx = omega * vy[i];
    My = -1.0 * omega * vx[i];
    Mz = 0.0;
    Ax = Fx[i]/m[i];
    Ay = Fy[i]/m[i];
    Az = Fz[i]/m[i];

    vx[i] += (Mx + Ax) *dt;
    vy[i] += (My + Ay) *dt;
    vz[i] += (Mz + Az) *dt;
}

void VV_update_position(int i, double dt)
{
    rx[i] += vx[i]*dt;
    ry[i] += vy[i]*dt;
    rz[i] += vz[i]*dt;
    
    if( rx[i] > Lx) rx[i] -= Lx;
    if( rx[i] < 0)  rx[i] += Lx;

    if( ry[i] > Ly) ry[i] -= Ly;
    if( ry[i] < 0)  ry[i] += Ly;

    if( rz[i] > Lz) rz[i] -= Lz;
    if( rz[i] < 0)  rz[i] += Lz;
}


void update_position(int i, double dt);
void update_velocity(int i, double dt);

void Step(double dt)
{
    int i;
    U = K = KelXY = KelZ = KpXY = KpZ = 0.0;
    for(i = 0; i < N; i++)
    {
        update_velocity(i, 0.5*dt);
        update_position(i, dt);
    }

    update_forces();

    for(i = 0; i < N; i++)
    {
        update_velocity(i, 0.5*dt);
        get_kinetic_energy(i);
    }
}


void update_velocity(int i, double dt)
{
	omega = (H * q[i])/(c * m[i]);
    Ex = Fx[i]/q[i]; Ey = Fy[i]/q[i]; Ez = Fz[i]/q[i]; 
	Vx = vx[i]; Vy = vy[i]; Vz = vz[i];

    if(H != 0.0)
    {
        I = Vx * cos(omega*dt);
        II = q[i]/(m[i]*omega)*Ey*(1 - cos(omega*dt));
        III = (Vy*omega + (q[i]/m[i])*Ex) * sin(omega*dt)/omega;
    }
    else
    {
        I = Vx;
        II = 0.0;
        III = (Vy*omega + (q[i]/m[i])*Ex) * dt;
    }
    vx[i] =  I + II + III;

    if(H != 0.0)
    {
        I = Vy * cos(omega*dt);
        II = -1.0 * q[i]/(m[i]*omega)*Ex*(1 - cos(omega*dt));
        III = -1.0 * (Vx*omega - (q[i]/m[i])*Ey) * sin(omega*dt)/omega;
    }
    else
    {
        I = Vy;
        II = 0.0;
        III = -1.0 * (Vx*omega - (q[i]/m[i])*Ey) * dt;
    }
    vy[i] = I + II + III;

    vz[i] = Vz + (q[i]/m[i])*Ez*dt;
    return;
}

void update_position(int i, double dt)
{
    omega = (H * q[i])/(c * m[i]);
    Ex = Fx[i]/q[i]; Ey = Fy[i]/q[i]; Ez = Fz[i]/q[i]; 

    if(H != 0.0)
    {
        I = vx[i] * (sin(omega*dt) / omega);
        II = (q[i]/m[i])*Ey * ((omega*dt - sin(omega*dt)) / (omega*omega));
        III = (vy[i]*omega + (q[i]/m[i])*Ex) * ((1.0-cos(omega*dt))/(omega*omega));
    }
    else
    {
        I = vx[i]*dt;
        II = 0.0;
        III = (vy[i]*omega + (q[i]/m[i])*Ex) * (0.5*dt*dt);
    }
    rx[i] += I + II + III;

    if(H != 0.0)
    {
        I = vy[i] * (sin(omega*dt) / omega);
        II = -1.0 * (q[i]/m[i])*Ex * ((omega*dt - sin(omega*dt)) / (omega*omega));
        III = -1.0 * (vx[i]*omega - (q[i]/m[i])*Ey) * ((1.0-cos(omega*dt))/(omega*omega));
    }
    else
    {
        I = vy[i]*dt;
        II = 0.0;
        III = -1.0 * (vx[i]*omega - (q[i]/m[i])*Ey) * (0.5*dt*dt);
    }
    ry[i] += I + II + III;

    rz[i] += vz[i]*dt + (q[i]/m[i])*Ez * 0.5*dt*dt;

    if( rx[i] > Lx) rx[i] -= Lx;
    if( rx[i] < 0)  rx[i] += Lx;

    if( ry[i] > Ly) ry[i] -= Ly;
    if( ry[i] < 0)  ry[i] += Ly;

    if( rz[i] > Lz) rz[i] -= Lz;
    if( rz[i] < 0)  rz[i] += Lz;
    return;
}


