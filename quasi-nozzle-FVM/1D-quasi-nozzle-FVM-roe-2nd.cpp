//***********************************************************************************************
// quasi-nozzle-FVM.cpp: Solve the Euler equations for the flow in a quasi-1D nozzle 
//                       using the finite volume method.
// Notes: central approximation for the fluxes; Euler explicit time integration method;
//        2nd and 4th order JST damping for stability
// © 2023 Xu-Hui Zhou All Rights Reserved.
//***********************************************************************************************

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <fstream>
#define EPSILON 1e-8
using namespace std;

// function declaration
void getGeo(int N, double x_min, double x_max, vector<double>& x_cell, vector<double>& x_face, vector<double>& area_face, vector<double>& vol_cell);
double func(double M, double gamma, double Abar);
double derivFunc(double M, double gamma, double Abar);
void getExactSol(double gamma, double Ru, double M_air, double p0, double T0, vector<double>& x_cell, vector<double>& M0_exact, vector<vector<double> >& V0_exact, vector<vector<double> >& U0_exact);
void getInitCons(double M, double gamma, double Ru, double M_air, double p0, double T0, vector<vector<double> >& U0);
// void consToFlux(double gamma, vector<vector<double> >& U0, vector<vector<double> >& F0);
void consToPrim(double gamma, vector<vector<double> >& U0, vector<vector<double> >& V0);
void primToLimiter(vector<vector<double> >& V0, vector<vector<double> >& phiplus, vector<vector<double> >& phiminus);
void getLeftRightStates(vector<vector<double> >& V0, vector<vector<double> >& phiplus, vector<vector<double> >& phiminus, vector<vector<double> >& VL, vector<vector<double> >& VR, double kappa, double epsilon);
void primToFlux(vector<vector<double> >& VL, vector<vector<double> >& VR, vector<vector<double> >& F0, double gamma);
// void getDamp(double gamma, double K2, double K4, vector<vector<double> >& U0, vector<vector<double> >& d0);
void getSource(double gamma, double delta_x, vector<vector<double> >& U0, vector<vector<double> >& S0, vector<double>& area_face);
void setInBoun(vector<vector<double> >& U0, double gamma, double p_stag, double T_stag, double Ru, double M_air);
void setOutBoun(vector<vector<double> >& U0, bool isPackPres, double pb, double gamma);
void getResidual(vector<vector<double> >& F0, vector<vector<double> >& S0, vector<vector<double> >& R0, vector<double>& area_face, double delta_x);
void getDeltaT(vector<vector<double> >& U0, vector<double>& t0, double CFL, double delta_x, double gamma);
void doUpdate(vector<vector<double> >& U0, vector<vector<double> >& U1, vector<vector<double> >& R0, vector<double>& vol_cell, double delta_t, double delta_x);
void getL2Norm(vector<vector<double> >& R0, vector<double>& current_resi_norm);
void getMach(double gamma, vector<vector<double> >& U0, vector<double>& M0);
void transpose2DVec(vector<vector<double> >& U0, vector<vector<double> >& U0T);


int main() {

    // ******** Geometry ********
    int num_cell; // needs to be even number
    double x_min, x_max, delta_x;
    num_cell = 50;
    x_min = -1.0;
    x_max = 1.0;
    delta_x = (x_max - x_min)/double(num_cell);

    vector<double> x_cell(num_cell), x_face(num_cell+1), area_face(num_cell+1), vol_cell(num_cell);
    getGeo(num_cell, x_min, x_max, x_cell, x_face, area_face, vol_cell);

    // ********************************* Computation *********************************
    vector<vector<double> > U0(3, vector<double>(num_cell+4)), U1(3, vector<double>(num_cell+4));
    vector<vector<double> > V0(3, vector<double>(num_cell+4));
    vector<vector<double> > F0(3, vector<double>(num_cell+1));
    vector<vector<double> > phiplus(3, vector<double>(num_cell+2)), phiminus(3, vector<double>(num_cell+2));
    vector<vector<double> > VL(3, vector<double>(num_cell+1)), VR(3, vector<double>(num_cell+1));
    vector<vector<double> > S0(3, vector<double>(num_cell));
    vector<vector<double> > R0(3, vector<double>(num_cell));
    vector<double> M0(num_cell+4);
    vector<double> t0(num_cell);
    
    double gamma, Ru, M_air, T0, p0, pb;
    double kappa, epsilon;
    bool isPackPres;
    string kappa_name;
    gamma = 1.4;
    Ru = 8314.0;
    M_air = 28.96;
    T0 = 600.0; // stagnation temperature
    p0 = 300000.0; // stagnation pressure
    isPackPres = true;
    pb = 120000.0; // back pressure
    epsilon = 1.0; // 0 for 1st order and 1 for 2nd order
    kappa = 0.0; // -1 for fully upwind and 0 for upwind biased
    kappa_name = "0p0";

    // ******** Exact Solutions (for isentropic case only) ********
    vector<double> M0_exact(num_cell);
    vector<vector<double> > V0_exact(3, vector<double>(num_cell)), U0_exact(3, vector<double>(num_cell));
    getExactSol(gamma, Ru, M_air, p0, T0, x_cell, M0_exact, V0_exact, U0_exact);

    // ******** Initial Conditions ********
    double M_init;
    M_init = 1.5;
    getInitCons(M_init, gamma, Ru, M_air, p0, T0, U0); // calculate initial conservative variable
    consToPrim(gamma, U0, V0); // get initial primitive variable
    primToLimiter(V0, phiplus, phiminus); // calcuate limiter phi
    getLeftRightStates(V0, phiplus, phiminus, VL, VR, kappa, epsilon); // calculate left and right states for primitive variables
    primToFlux(VL, VR, F0, gamma); // left and right primitives to flux
    getSource(gamma, delta_x, U0, S0, area_face);

    // ******** Solve Equations Iteratively ********
    int num_iter;
    num_iter = 200000;
    double CFL, delta_t;
    CFL = 0.1;
    vector<vector<double> > resi_norm(num_iter, vector<double>(3, 0.0));
    
    for (int i = 0; i < num_iter; i++){
        getResidual(F0, S0, R0, area_face, delta_x);
        getL2Norm(R0,resi_norm[i]);
        // cout << resi_norm[i][0] << " " << resi_norm[i][1] << " " << resi_norm[i][2] << endl;
        if ((resi_norm[i][0]/resi_norm[0][0] < EPSILON && resi_norm[i][1]/resi_norm[0][1] < EPSILON) && (resi_norm[i][2]/resi_norm[0][2] < EPSILON)) {
            cout << resi_norm[i][0] << " " << resi_norm[i][1] << " " << resi_norm[i][2] << endl;
            cout << resi_norm[i][0]/resi_norm[0][0] << " " << resi_norm[i][1]/resi_norm[0][1] << " " << resi_norm[i][2]/resi_norm[0][2] << endl;
            cout << "The number of steps for convergence is " << i << endl;
            break;
        }
        getDeltaT(U0, t0, CFL, delta_x, gamma);
        delta_t = *min_element(t0.begin(), t0.end()); // global minimal time step
        doUpdate(U0, U1, R0, vol_cell, delta_t, delta_x);
        setInBoun(U0, gamma, p0, T0, Ru, M_air);
        setOutBoun(U0, isPackPres, pb, gamma);
        consToPrim(gamma, U0, V0);
        primToLimiter(V0, phiplus, phiminus);
        getLeftRightStates(V0, phiplus, phiminus, VL, VR, kappa, epsilon);
        primToFlux(VL, VR, F0, gamma);
        getSource(gamma, delta_x, U0, S0, area_face);
    }
    // consToPrim(gamma, U0, V0);
    getMach(gamma, V0, M0);

    // ******** Write results into files ********
    // string fileResi = "residual_" + to_string(num_cell) + ".txt";
    // ofstream writeResi("2nd-order-shock/Roe/"+fileResi); // write residual with iteration
    // for (const auto& row : resi_norm) {
    //     for (const auto& element : row) {
    //         writeResi << element << ' ';
    //     }
    //     writeResi << '\n';
    // }
    // writeResi.close();

    // vector<vector<double> > V0T(num_cell+4, vector<double>(3)); // write primitive variable
    // transpose2DVec(V0,V0T);
    // string fileV = "V_" + to_string(num_cell) + ".txt";
    // ofstream writePrim("2nd-order-shock/Roe/"+fileV);
    // for (const auto& row : V0T) {
    //     for (const auto& element : row) {
    //         writePrim << element << ' ';
    //     }
    //     writePrim << '\n';
    // }
    // writePrim.close();

    // vector<vector<double> > U0T(num_cell+4, vector<double>(3)); // write primitive variable
    // transpose2DVec(U0,U0T);
    // string fileU = "U_" + to_string(num_cell) + ".txt";
    // ofstream writeCons("2nd-order-shock/Roe/"+fileU);
    // for (const auto& row : U0T) {
    //     for (const auto& element : row) {
    //         writeCons << element << ' ';
    //     }
    //     writeCons << '\n';
    // }
    // writeCons.close();

    string fileMa = "Ma_" + to_string(num_cell) + "_" + kappa_name + ".txt";
    ofstream writeMach("2nd-order-shock-diff-kappa/Roe/"+fileMa); // write Mach number
    for (const auto& element : M0) {
        writeMach << element << '\n';
    }
    writeMach.close();

    string fileXCell = "xcell_" + to_string(num_cell) + ".txt";
    ofstream writeX("2nd-order-shock-diff-kappa/Roe/"+fileXCell); // write x cell
    for (const auto& element : x_cell) {
        writeX << element << '\n';
    }
    writeX.close();

    // vector<vector<double> > V0_exact_T(num_cell, vector<double>(3)); // write exact primitive variable
    // transpose2DVec(V0_exact,V0_exact_T);
    // string fileVExact = "V_exact_" + to_string(num_cell) + ".txt";
    // ofstream writeExactPrim("2nd-order-shock/Roe/"+fileVExact);
    // for (const auto& row : V0_exact_T) {
    //     for (const auto& element : row) {
    //         writeExactPrim << element << ' ';
    //     }
    //     writeExactPrim << '\n';
    // }
    // writeExactPrim.close();

    // vector<vector<double> > U0_exact_T(num_cell, vector<double>(3)); // write exact primitive variable
    // transpose2DVec(U0_exact,U0_exact_T);
    // string fileUExact = "U_exact_" + to_string(num_cell) + ".txt";
    // ofstream writeExactCons("2nd-order-shock/Roe/"+fileUExact);
    // for (const auto& row : U0_exact_T) {
    //     for (const auto& element : row) {
    //         writeExactCons << element << ' ';
    //     }
    //     writeExactCons << '\n';
    // }
    // writeExactCons.close();

    // string fileMaExact = "Ma_exact_" + to_string(num_cell) + ".txt";
    // ofstream writeExactMach("2nd-order-shock/Roe/"+fileMaExact); // write Mach number
    // for (const auto& element : M0_exact) {
    //     writeExactMach << element << '\n';
    // }
    // writeExactMach.close();
    
    return 0;
}

// Function: To get cell center coordinates, face coordinates, face areas, and cell volumes
void getGeo(int N, double x_min, double x_max, vector<double>& x_cell, vector<double>& x_face, vector<double>& area_face, vector<double>& vol_cell)
{
    for(int i = 0; i < N; i++){
        x_cell[i] = x_min + (double(i)+0.5)/double(N) * (x_max - x_min);
    }
    double PI = 2.0 * acos(0.0);
    for(int i = 0; i < N + 1; i++){
        x_face[i] = x_min + double(i)/double(N) * (x_max - x_min);
        area_face[i] = 0.2 + 0.4 * (1 + sin(PI*(x_face[i]-0.5)));  // A(x) = 0.2 + 0.4 * [1 + sin(pi{x-0.5})]
    }
    for (int i = 0; i < N; i++){
        vol_cell[i] = 0.5 * (area_face[i]+area_face[i+1]) * (x_max-x_min)/double(N);
    }
}

// Function: To get initial conservative variables (rho, rho*u, rho*et)
void getInitCons(double M, double gamma, double Ru, double M_air, double p0, double T0, vector<vector<double> >& U0)
{    
    double phi, T, p, rho, u;
    phi = 1.0+(gamma-1.0)/2.0*pow(M,2.0);
    T = T0 / phi;
    p = p0 / pow(phi, gamma/(gamma-1.0));
    rho = p / (Ru/M_air * T);
    u = M * sqrt(gamma * Ru/M_air * T);
    
    for (int j = 0; j < U0[0].size(); j++){
        U0[0][j] = rho;
        U0[1][j] = rho * u;
        U0[2][j] = 1.0/(gamma-1.0)*p + 0.5*rho*pow(u,2.0);
    }
}

// Function: To convert conservative variable to primitive variable
void consToPrim(double gamma, vector<vector<double> >& U0, vector<vector<double> >& V0)
{
    for (int j = 0; j < V0[0].size(); j++){
        V0[0][j] = U0[0][j];
        V0[1][j] = U0[1][j]/(U0[0][j]+1e-16);
        V0[2][j] = (U0[2][j] - 0.5*pow(U0[1][j],2.0)/(U0[0][j]+1e-16)) * (gamma-1.0);
    }
    int N;
    N = V0[0].size();
    V0[0][0] = max(V0[0][0], 0.0001);
    V0[1][0] = max(V0[1][0], 10.0);
    V0[2][0] = max(V0[2][0], 500.0);
    V0[0][1] = max(V0[0][1], 0.0001);
    V0[1][1] = max(V0[1][1], 10.0);
    V0[2][1] = max(V0[2][1], 500.0);
    V0[0][N-2] = max(V0[0][N-2], 0.0001);
    V0[1][N-2] = max(V0[1][N-2], 10.0);
    V0[2][N-2] = max(V0[2][N-2], 500.0);
    V0[0][N-1] = max(V0[0][N-1], 0.0001);
    V0[1][N-1] = max(V0[1][N-1], 10.0);
    V0[2][N-1] = max(V0[2][N-1], 500.0);
}

// Function: To get Mach number from primitive variable
void getMach(double gamma, vector<vector<double> >& V0, vector<double>& M0)
{
    for (int j = 0; j < M0.size(); j++){
        M0[j] = fabs(V0[1][j])/(sqrt(gamma*V0[2][j]/(V0[0][j]+1e-16))+1e-16);
    }
}

// Function: To calculate source term from conservative variable
void getSource(double gamma, double delta_x, vector<vector<double> >& U0, vector<vector<double> >& S0, vector<double>& area_face){
    for (int j = 0; j < S0[0].size(); j++){
        S0[0][j] = 0.0;
        S0[1][j] = (U0[2][j+2]-0.5*pow(U0[1][j+2],2.0)/(U0[0][j+2]+1e-16)) * (gamma - 1.0) * (area_face[j+1]-area_face[j])/delta_x;
        S0[2][j] = 0.0;
    }
}

// Function: To set inflow boundary conditions
void setInBoun(vector<vector<double> >& U0, double gamma, double p_stag, double T_stag, double Ru, double M_air){
    double u2, p2, a2, M2, u3, p3, a3, M3, M0, M1, phi0, phi1, T0, T1, p0, p1, rho0, rho1, u0, u1;
    // inflow boundary conditions
    u2 = U0[1][2]/(U0[0][2]+1e-16);
    p2 = (U0[2][2]-0.5*pow(U0[1][2],2.0)/(U0[0][2]+1e-16)) * (gamma - 1.0);
    a2 = sqrt(gamma*p2/(U0[0][2]+1e-16));
    M2 = fabs(u2)/(a2+1e-16);
    u3 = U0[1][3]/(U0[0][3]+1e-16);
    p3 = (U0[2][3]-0.5*pow(U0[1][3],2.0)/(U0[0][3]+1e-16)) * (gamma - 1.0);
    a3 = sqrt(gamma*p3/(U0[0][3]+1e-16));
    M3 = fabs(u3)/(a3+1e-16);
    M0 = max(3.0*M2 - 2.0*M3,0.005);
    M1 = max(2.0*M2 - M3,0.005);
    phi0 = 1.0 + (gamma-1.0)*pow(M0,2.0)/2.0;
    phi1 = 1.0 + (gamma-1.0)*pow(M1,2.0)/2.0;
    T0 = T_stag / (phi0+1e-16);
    T1 = T_stag / (phi1+1e-16);
    p0 = p_stag / (pow(phi0,gamma/(gamma-1.0))+1e-16);
    p1 = p_stag / (pow(phi1,gamma/(gamma-1.0))+1e-16);
    rho0 = p0 / (Ru/M_air * T0 + 1e-16);
    rho1 = p1 / (Ru/M_air * T1 + 1e-16);
    u0 = M0 * sqrt(gamma * Ru/M_air * T0);
    u1 = M1 * sqrt(gamma * Ru/M_air * T1);
    U0[0][0] = rho0;
    U0[1][0] = rho0 * u0;
    U0[2][0] = 1.0/(gamma-1.0)*p0 + 0.5*rho0*pow(u0,2.0);
    U0[0][1] = rho1;
    U0[1][1] = rho1 * u1;
    U0[2][1] = 1.0/(gamma-1.0)*p1 + 0.5*rho1*pow(u1,2.0);
}

// Function: To set outflow boundary conditions
void setOutBoun(vector<vector<double> >& U0, bool isPackPres, double pb, double gamma){
    int N;
    double um4, um3, um2, um1, pm3, pm2, pm1;
    N = U0[0].size();
    if (isPackPres == false) {
        for (int i = 0; i < U0.size(); i++){
            U0[i][N-2] = 2.0*U0[i][N-3] - U0[i][N-4];
            U0[i][N-1] = 3.0*U0[i][N-3] - 2.0*U0[i][N-4];
        }
    } else {
        U0[0][N-2] = 2.0*U0[0][N-3] - U0[0][N-4];
        U0[0][N-1] = 3.0*U0[0][N-3] - 2.0*U0[0][N-4];
        um4 = U0[1][N-4]/(U0[0][N-4]+1e-16);
        um3 = U0[1][N-3]/(U0[0][N-3]+1e-16);
        um2 = 2.0*um3 - um4;
        um1 = 3.0*um3 - 2.0*um4;
        U0[1][N-2] = U0[0][N-2]*um2;
        U0[1][N-1] = U0[0][N-1]*um1;
        pm3 = (U0[2][N-3]-0.5*pow(U0[1][N-3],2.0)/(U0[0][N-3]+1e-16)) * (gamma - 1.0);
        pm2 = 2.0*pb - pm3;
        pm1 = 4.0*pb - 3.0*pm3;
        U0[2][N-2] = 1.0/(gamma-1.0)*pm2 + 0.5*pow(U0[1][N-2],2.0)/(U0[0][N-2]+1e-16);
        U0[2][N-1] = 1.0/(gamma-1.0)*pm1 + 0.5*pow(U0[1][N-1],2.0)/(U0[0][N-1]+1e-16);
    }
}

// Function: To get residual
void getResidual(vector<vector<double> >& F0, vector<vector<double> >& S0, vector<vector<double> >& R0, vector<double>& area_face, double delta_x){
    for (int i = 0; i < R0.size(); i++){
        for (int j = 0; j < R0[0].size(); j++){
            R0[i][j] = F0[i][j+1]*area_face[j+1] - F0[i][j]*area_face[j] - S0[i][j]*delta_x;
        }
    }
}

// Function: To get the L2 norm of residual
void getL2Norm(vector<vector<double> >& R0, vector<double>& current_resi_norm) {
    vector<double> norm(3, 0.0);
    for (int i = 0; i < R0.size(); i++){
        for (int j = 0; j < R0[0].size(); j++){
            norm[i] += pow(R0[i][j],2.0);
        }
    }
    for (int i = 0; i < current_resi_norm.size(); i++){
        current_resi_norm[i] = sqrt(norm[i]/double(R0[0].size()));
    }
}

// Function: To get local time step
void getDeltaT(vector<vector<double> >& U0, vector<double>& t0, double CFL, double delta_x, double gamma){
    double u, p, a;
    for (int j = 0; j < t0.size(); j++){
        u = U0[1][j+2]/(U0[0][j+2]+1e-16);
        p = (U0[2][j+2]-0.5*pow(U0[1][j+2],2.0)/(U0[0][j+2]+1e-16)) * (gamma - 1.0);
        a = sqrt(gamma*p/(U0[0][j+2]+1e-16));
        t0[j] = CFL*delta_x/(fabs(u)+a+1e-16);
    }
}

// Function: To do iteration
void doUpdate(vector<vector<double> >& U0, vector<vector<double> >& U1, vector<vector<double> >& R0, vector<double>& vol_cell, double delta_t, double delta_x){
    for (int i = 0; i < U0.size(); i++){
        for (int j = 2; j < U0[0].size()-2; j++){
            U1[i][j] = U0[i][j] - (delta_t/vol_cell[j-2])*R0[i][j-2];
        }
    }
    U0 = U1;
}

// Function: Transpose a 2D vector
void transpose2DVec(vector<vector<double> >& U0, vector<vector<double> >& U0T){
    for (int i = 0; i < U0T.size(); i++){
        for (int j = 0; j < U0T[0].size(); j++){
            U0T[i][j] = U0[j][i];
        }
    }
}

// Functions: To get exact solution for isentropic case
double func(double M, double gamma, double Abar) {
    return pow(2.0/(gamma+1.0)*(1.0+(gamma-1.0)/2.0*pow(M,2.0)),(gamma+1.0)/(gamma-1.0)) - pow(Abar,2.0)*pow(M,2.0);
}

double derivFunc(double M, double gamma, double Abar) {
    return 2.0 * M * (pow(2.0/(gamma+1.0)*(1.0+(gamma-1.0)/2.0*pow(M,2)),2.0/(gamma-1.0)) - pow(Abar,2.0));
}

void getExactSol(double gamma, double Ru, double M_air, double p0, double T0, vector<double>& x_cell, vector<double>& M0_exact, vector<vector<double> >& V0_exact, vector<vector<double> >& U0_exact)
{
    double PI, Abar, h;
    PI = 2 * acos(0.0);
    for (int i = 0; i < x_cell.size(); i++){
        // Give initial Mach number
        if (x_cell[i] < 0.0){
            M0_exact[i] = 0.5;
        } else {
            M0_exact[i] = 3.0;
        }
        Abar = (0.2 + 0.4 * (1.0 + sin(PI*(x_cell[i]-0.5))))/0.2; // Calculate A/A*
        // Newton Iteration
        h = func(M0_exact[i], gamma, Abar)/(derivFunc(M0_exact[i], gamma, Abar)+1e-16);
        while (fabs(h) >= 1e-11)
        {
            M0_exact[i] = M0_exact[i] - h;
            h = func(M0_exact[i], gamma, Abar)/(derivFunc(M0_exact[i], gamma, Abar)+1e-16);
        }
        // get exact solution of primitive variable from solved Mach number
        double phi, T_curr;
        phi = 1.0 + (gamma-1.0)/2.0 * pow(M0_exact[i],2.0);
        T_curr = T0/(phi+1e-16);
        V0_exact[2][i] = p0/(pow(phi,gamma/(gamma-1.0))+1e-16);
        V0_exact[0][i] = V0_exact[2][i]/((Ru/M_air*T_curr)+1e-16);
        V0_exact[1][i] = M0_exact[i] * sqrt(gamma*Ru/M_air*T_curr);
        // convert V0_exact to U0_exact
        U0_exact[0][i] = V0_exact[0][i];
        U0_exact[1][i] = V0_exact[0][i]*V0_exact[1][i];
        U0_exact[2][i] = 1.0/(gamma-1.0)*V0_exact[2][i] + 0.5*V0_exact[0][i]*pow(V0_exact[1][i],2.0);
    }
}

// // Function: To convert primitive variable to flux
// void primToFlux(vector<vector<double> >& V0, vector<vector<double> >& F0, double gamma)
// {
//     double R, rhobar, ubar, ht_L, ht_R, htbar, abar;
//     double lambda1bar, lambda2bar, lambda3bar, lambda1abs, lambda2abs, lambda3abs;
//     double deltarho, deltau, deltap, deltaw1, deltaw2, deltaw3;
//     vector<double> r1bar(3), r2bar(3), r3bar(3), fL(3), fR(3);

//     for (int j = 0; j < F0[0].size(); j++){
//         R = sqrt(V0[0][j+2]/(V0[0][j+1]+1e-16));
//         rhobar = R*V0[0][j+1];
//         ubar = (R*V0[1][j+2]+V0[1][j+1])/(R+1.0);
//         ht_L = gamma/(gamma-1.0)*V0[2][j+1]/(V0[0][j+1]+1e-16) + 0.5*pow(V0[1][j+1],2.0);
//         ht_R = gamma/(gamma-1.0)*V0[2][j+2]/(V0[0][j+2]+1e-16) + 0.5*pow(V0[1][j+2],2.0);
//         htbar = (R*ht_R+ht_L)/(R+1.0);
//         abar = sqrt((gamma-1.0)*(htbar-0.5*pow(ubar,2.0)));

//         lambda1bar = ubar;
//         lambda2bar = ubar + abar;
//         lambda3bar = ubar - abar;
//         if (fabs(lambda1bar) >= 2.0*0.1*abar) {
//             lambda1abs = fabs(lambda1bar);
//         } else {
//             lambda1abs = pow(lambda1bar,2.0)/(4.0*0.1*abar+1e-16) + 0.1*abar;
//         }
//         if (fabs(lambda2bar) >= 2.0*0.1*abar) {
//             lambda2abs = fabs(lambda2bar);
//         } else {
//             lambda2abs = pow(lambda2bar,2.0)/(4.0*0.1*abar+1e-16) + 0.1*abar;
//         }
//         if (fabs(lambda3bar) >= 2.0*0.1*abar) {
//             lambda3abs = fabs(lambda3bar);
//         } else {
//             lambda3abs = pow(lambda3bar,2.0)/(4.0*0.1*abar+1e-16) + 0.1*abar;
//         }

//         deltarho = V0[0][j+2] - V0[0][j+1];
//         deltau = V0[1][j+2] - V0[1][j+1];
//         deltap = V0[2][j+2] - V0[2][j+1];
//         deltaw1 = deltarho - deltap/(pow(abar,2.0)+1e-16);
//         deltaw2 = deltau + deltap/(rhobar*abar+1e-16);
//         deltaw3 = deltau - deltap/(rhobar*abar+1e-16);

//         r1bar[0] = 1.0;
//         r1bar[1] = ubar;
//         r1bar[2] = 0.5*pow(ubar,2.0);
//         r2bar[0] = rhobar/(2.0*abar+1e-16);
//         r2bar[1] = rhobar/(2.0*abar+1e-16)*(ubar+abar);
//         r2bar[2] = rhobar/(2.0*abar+1e-16)*(htbar+ubar*abar);
//         r3bar[0] = -rhobar/(2.0*abar+1e-16);
//         r3bar[1] = -rhobar/(2.0*abar+1e-16)*(ubar-abar);
//         r3bar[2] = -rhobar/(2.0*abar+1e-16)*(htbar-ubar*abar);

//         fL[0] = V0[0][j+1]*V0[1][j+1];
//         fL[1] = V0[0][j+1]*pow(V0[1][j+1],2.0)+V0[2][j+1];
//         fL[2] = V0[0][j+1]*V0[1][j+1]*ht_L;
//         fR[0] = V0[0][j+2]*V0[1][j+2];
//         fR[1] = V0[0][j+2]*pow(V0[1][j+2],2.0)+V0[2][j+2];
//         fR[2] = V0[0][j+2]*V0[1][j+2]*ht_R;

//         for (int i = 0; i < 3; i++){
//             F0[i][j] = 0.5*(fL[i]+fR[i]) - 0.5*(lambda1abs*deltaw1*r1bar[i]+lambda2abs*deltaw2*r2bar[i]+lambda3abs*deltaw3*r3bar[i]);
//         }
//     }
// }

// Function: To calculate the limiting function phiplus and phiminus from primitive variables
void primToLimiter(vector<vector<double> >& V0, vector<vector<double> >& phiplus, vector<vector<double> >& phiminus) 
{
    vector<double> rplus(3), rminus(3);
    for (int j = 0; j < phiplus[0].size(); j++){
        for (int i = 0; i < 3; i++){
            rplus[i] = (V0[i][j+2]-V0[i][j+1])/(copysign(1.0,V0[i][j+1]-V0[i][j])*max(fabs(V0[i][j+1]-V0[i][j]),1e-6));
            phiplus[i][j] = (rplus[i]+fabs(rplus[i]))/(1.0+rplus[i]+1e-16);

            rminus[i] = (V0[i][j+1]-V0[i][j])/(copysign(1.0,V0[i][j+2]-V0[i][j+1])*max(fabs(V0[i][j+2]-V0[i][j+1]),1e-6));
            phiminus[i][j] = (rminus[i]+fabs(rminus[i]))/(1.0+rminus[i]+1e-16);
        }
    }
}

// Function: To calculate left and right states
void getLeftRightStates(vector<vector<double> >& V0, vector<vector<double> >& phiplus, vector<vector<double> >& phiminus, vector<vector<double> >& VL, vector<vector<double> >& VR, double kappa, double epsilon)
{
    for (int i = 0; i < VL.size(); i++){
        for (int j = 0; j < VL[0].size(); j++){
            VL[i][j] = V0[i][j+1] + epsilon/4.0*((1.0-kappa)*phiplus[i][j]*(V0[i][j+1]-V0[i][j])+(1.0+kappa)*phiminus[i][j]*(V0[i][j+2]-V0[i][j+1]));
            VR[i][j] = V0[i][j+2] - epsilon/4.0*((1.0-kappa)*phiminus[i][j+1]*(V0[i][j+3]-V0[i][j+2])+(1.0+kappa)*phiplus[i][j+1]*(V0[i][j+2]-V0[i][j+1]));
        }
    }
}

void primToFlux(vector<vector<double> >& VL, vector<vector<double> >& VR, vector<vector<double> >& F0, double gamma)
{
    double R, rhobar, ubar, ht_L, ht_R, htbar, abar;
    double lambda1bar, lambda2bar, lambda3bar, lambda1abs, lambda2abs, lambda3abs;
    double deltarho, deltau, deltap, deltaw1, deltaw2, deltaw3;
    vector<double> r1bar(3), r2bar(3), r3bar(3), fL(3), fR(3);

    for (int j = 0; j < F0[0].size(); j++){
        R = sqrt(VR[0][j]/(VL[0][j]+1e-16));
        rhobar = R*VL[0][j];
        ubar = (R*VR[1][j]+VL[1][j])/(R+1.0);
        ht_L = gamma/(gamma-1.0)*VL[2][j]/(VL[0][j]+1e-16) + 0.5*pow(VL[1][j],2.0);
        ht_R = gamma/(gamma-1.0)*VR[2][j]/(VR[0][j]+1e-16) + 0.5*pow(VR[1][j],2.0);
        htbar = (R*ht_R+ht_L)/(R+1.0);
        abar = sqrt((gamma-1.0)*(htbar-0.5*pow(ubar,2.0)));

        lambda1bar = ubar;
        lambda2bar = ubar + abar;
        lambda3bar = ubar - abar;
        if (fabs(lambda1bar) >= 2.0*0.1*abar) {
            lambda1abs = fabs(lambda1bar);
        } else {
            lambda1abs = pow(lambda1bar,2.0)/(4.0*0.1*abar+1e-16) + 0.1*abar;
        }
        if (fabs(lambda2bar) >= 2.0*0.1*abar) {
            lambda2abs = fabs(lambda2bar);
        } else {
            lambda2abs = pow(lambda2bar,2.0)/(4.0*0.1*abar+1e-16) + 0.1*abar;
        }
        if (fabs(lambda3bar) >= 2.0*0.1*abar) {
            lambda3abs = fabs(lambda3bar);
        } else {
            lambda3abs = pow(lambda3bar,2.0)/(4.0*0.1*abar+1e-16) + 0.1*abar;
        }

        deltarho = VR[0][j] - VL[0][j];
        deltau = VR[1][j] - VL[1][j];
        deltap = VR[2][j] - VL[2][j];
        deltaw1 = deltarho - deltap/(pow(abar,2.0)+1e-16);
        deltaw2 = deltau + deltap/(rhobar*abar+1e-16);
        deltaw3 = deltau - deltap/(rhobar*abar+1e-16);

        r1bar[0] = 1.0;
        r1bar[1] = ubar;
        r1bar[2] = 0.5*pow(ubar,2.0);
        r2bar[0] = rhobar/(2.0*abar+1e-16);
        r2bar[1] = rhobar/(2.0*abar+1e-16)*(ubar+abar);
        r2bar[2] = rhobar/(2.0*abar+1e-16)*(htbar+ubar*abar);
        r3bar[0] = -rhobar/(2.0*abar+1e-16);
        r3bar[1] = -rhobar/(2.0*abar+1e-16)*(ubar-abar);
        r3bar[2] = -rhobar/(2.0*abar+1e-16)*(htbar-ubar*abar);

        fL[0] = VL[0][j]*VL[1][j];
        fL[1] = VL[0][j]*pow(VL[1][j],2.0)+VL[2][j];
        fL[2] = VL[0][j]*VL[1][j]*ht_L;
        fR[0] = VR[0][j]*VR[1][j];
        fR[1] = VR[0][j]*pow(VR[1][j],2.0)+VR[2][j];
        fR[2] = VR[0][j]*VR[1][j]*ht_R;

        for (int i = 0; i < 3; i++){
            F0[i][j] = 0.5*(fL[i]+fR[i]) - 0.5*(lambda1abs*deltaw1*r1bar[i]+lambda2abs*deltaw2*r2bar[i]+lambda3abs*deltaw3*r3bar[i]);
        }
    }
}

