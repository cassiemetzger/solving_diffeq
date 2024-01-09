#ifndef LAX_WENDROFF_H
#define LAX_WENDROFF_H
#include <vector>
#include <iostream>
#include <cmath>

std::vector<double> lax_wendroff(const std::vector<double> &u, double v, double dx, double dt){
    int Nx = u.size(); 
    std::vector<double> u_new(Nx, 0.0); //fill vector u_new with 0s while giving it length of u 
    //ghost cells 
    u_new[0] = u[0]; 
    u_new[0] = u_new[Nx-1];
    double u_half_before = 0.5 * (u[0] + u[1]) - v*dt/(2.0*dx)*(u[1] - u[0]);
    double u_half_after;
    for(int j = 1; j < Nx-1; j++){
        u_half_after = 0.5 *(u[j] + u[j+1])- v*dt/(2.0*dx)*(u[j+1] - u[j]);
        u_new[j] = u[j] - v*dt/dx*(u_half_after - u_half_before);
        u_half_before = u_half_after;
    }
    //apply boundary conditions
    u_new[0] = u_new[Nx-2];
    u_new[Nx-1] = u_new[1];

    return u_new;
}
#endif
                