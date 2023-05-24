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
void getGeo(vector<vector<double> >& x_node, vector<vector<double> >& y_node, vector<vector<double> >& x_cell, vector<vector<double> >& y_cell, 
vector<vector<double> >& Area_x, vector<vector<double> >& Area_y, vector<vector<double> >& Vol);
void getManufSol(vector<vector<vector<double> > >& VE, vector<vector<vector<double> > >& SE, vector<vector<double> >& VLeft, vector<vector<double> >& VRight,
vector<vector<double> >& VTop, vector<vector<double> >& VBottom, vector<vector<double> >& x_cell, vector<vector<double> >& y_cell,
vector<vector<double> >& x_node, vector<vector<double> >& y_node);
void getInitCons(double gamma, vector<vector<vector<double> > >& U0, vector<vector<vector<double> > >& V0);
void setBounCond(double gamma, vector<vector<vector<double> > >& V0, vector<vector<double> >& VLeft, vector<vector<double> >& VRight, vector<vector<double> >& VBottom, vector<vector<double> >& VTop);
void consToPrim(double gamma, vector<vector<vector<double> > >& U0, vector<vector<vector<double> > >& V0);
// void limitGlobal(vector<vector<vector<double> > >& V0);
void primToLimiter(vector<vector<vector<double> > >& V0, vector<vector<vector<double> > >& phiplus_x, vector<vector<vector<double> > >& phiminus_x,
vector<vector<vector<double> > >& phiplus_y, vector<vector<vector<double> > >& phiminus_y, bool islimiter);
void getLeftRightStates(vector<vector<vector<double> > >& V0, vector<vector<vector<double> > >& phiplus_x, vector<vector<vector<double> > >& phiminus_x, 
vector<vector<vector<double> > >& VL_x, vector<vector<vector<double> > >& VR_x,
vector<vector<vector<double> > >& phiplus_y, vector<vector<vector<double> > >& phiminus_y, 
vector<vector<vector<double> > >& VL_y, vector<vector<vector<double> > >& VR_y, double kappa, double epsilon);
void primToFlux(vector<vector<vector<double> > >& VL_x, vector<vector<vector<double> > >& VR_x, vector<vector<vector<double> > >& VL_y, vector<vector<vector<double> > >& VR_y, 
vector<vector<vector<double> > >& Fx, vector<vector<vector<double> > >& Fy, vector<vector<double> >& x_node, vector<vector<double> >& y_node, double gamma);
void getResidual(vector<vector<vector<double> > >& Fx, vector<vector<vector<double> > >& Fy, vector<vector<vector<double> > >& R0, vector<vector<double> >& Area_x, vector<vector<double> >& Area_y,
vector<vector<vector<double> > >& SE, vector<vector<double> >& Vol);
void getL2Norm(vector<vector<vector<double> > >& R0, vector<double>& current_resi_norm);
void getDeltaT(vector<vector<vector<double> > >& V0, vector<vector<double> >& Vol, vector<vector<double> >& x_node, vector<vector<double> >& y_node, 
vector<vector<double> >& t0, double & delta_t, double CFL, double gamma);
void doUpdate1(vector<vector<vector<double> > >& U0, vector<vector<vector<double> > >& U1, vector<vector<vector<double> > >& R0, vector<vector<double> >& Vol, double alpha1, double delta_t);
void doUpdate2(vector<vector<vector<double> > >& U0, vector<vector<vector<double> > >& U1, vector<vector<vector<double> > >& R0, vector<vector<double> >& Vol, double alpha2, double delta_t);
void doUpdate(vector<vector<vector<double> > >& U0, vector<vector<vector<double> > >& U1, vector<vector<vector<double> > >& R0, vector<vector<double> >& Vol, double delta_t);


int main() {

    // ************************************************ Geometry ************************************************

    // Declare variables
    int nzones, imax, jmax, kmax;
    int nx_cell, ny_cell;
    string cell_name;

    cell_name = "2d129";

    // Open file
    // std::ifstream inFile("../Project_Files/Grids/curviliniear-grids/curv2d257.grd");
    std::ifstream inFile("../Project_Files/Grids/curviliniear-grids/curv"+cell_name+".grd");

    // Read in data
    inFile >> nzones;
    inFile >> imax >> jmax >> kmax; // ith column & jth row !!!
    // cout << nzones << " " << imax << " " << jmax << " " << kmax << endl;

    // calculate cell numbers in x and y directions
    nx_cell = imax - 1;
    ny_cell = jmax - 1;

    // Read in x-coordinate and y-coordinate
    vector<vector<double> > x_node(imax, vector<double>(jmax));
    vector<vector<double> > y_node(imax, vector<double>(jmax));
    vector<vector<double> > Area_x(imax, vector<double>(jmax-1));
    vector<vector<double> > Area_y(imax-1, vector<double>(jmax));
    vector<vector<double> > x_cell(nx_cell, vector<double>(ny_cell));
    vector<vector<double> > y_cell(nx_cell, vector<double>(ny_cell));
    vector<vector<double> > Vol(nx_cell, vector<double>(ny_cell));

    for (int k = 0; k < kmax; k++) {
        for (int j = 0; j < jmax; j++) {
            for (int i = 0; i < imax; i++) {
                inFile >> x_node[imax-1-i][jmax-1-j];
            }
        }
    }

    for (int k = 0; k < kmax; k++) {
        for (int j = 0; j < jmax; j++) {
            for (int i = 0; i < imax; i++) {
                inFile >> y_node[imax-1-i][jmax-1-j];
            }
        }
    }

    // cout << x_node[0][0] << " " << x_node[2][0] << " " << x_node[4][0] << " " << x_node[6][0] << " " << x_node[8][0] << endl;
    // cout << y_node[0][0] << " " << y_node[2][0] << " " << y_node[4][0] << " " << y_node[6][0] << " " << y_node[8][0] << endl;
    // cout << x_node[0][8] << " " << x_node[2][8] << " " << x_node[4][8] << " " << x_node[6][8] << " " << x_node[8][8] << endl;
    // cout << y_node[0][8] << " " << y_node[2][8] << " " << y_node[4][8] << " " << y_node[6][8] << " " << y_node[8][8] << endl;

    // Close file
    inFile.close();

    // calculate face areas and cell volumes
    getGeo(x_node, y_node, x_cell, y_cell, Area_x, Area_y, Vol); // Have checked correctness!

    // ********************************* Manufactured Solution *********************************
    vector<vector<vector<double> > > VE(4, vector<vector<double> >(nx_cell, vector<double>(ny_cell))); // primitive variable at cell centers
    vector<vector<vector<double> > > SE(4, vector<vector<double> >(nx_cell, vector<double>(ny_cell))); // source terms at cell centers
    vector<vector<double> > VLeft(4, vector<double>(ny_cell)); // primitive variable at face centers on the left boundary
    vector<vector<double> > VRight(4, vector<double>(ny_cell)); // primitive variable at face centers on the right boundary
    vector<vector<double> > VBottom(4, vector<double>(nx_cell)); // primitive variable at face centers on the bottom boundary
    vector<vector<double> > VTop(4, vector<double>(nx_cell)); // primitive variable at face centers on the top boundary
    getManufSol(VE, SE, VLeft, VRight, VTop, VBottom, x_cell, y_cell, x_node, y_node);

    // write manufactured solution to tecplot file
    ofstream outfile("supersonic-vanleer-result/MMS-"+cell_name+".dat");
    outfile << "TITLE = \"Manufactured Solution\"\n";
    outfile << "variables=\"x(m)\" \"y(m)\" \"rho(kg/m^3)\" \"u(m/s)\" \"v(m/s)\" \"press(N/m^2)\" \"mass\" \"xmom\" \"ymom\" \"energy\"\n";

    outfile << "zone T=\"" << 1 << "\"\n";
    outfile << "I=" << nx_cell << " J=" << ny_cell << "\n";
    outfile << "DATAPACKING=POINT\n";

    for (int j = 0; j < x_cell[0].size(); j++) {
        for (int i = 0; i < x_cell.size(); i++) {
            outfile << x_cell[i][j] << " " << y_cell[i][j] << " ";
            outfile << VE[0][i][j] << " " << VE[1][i][j] << " ";
            outfile << VE[2][i][j] << " " << VE[3][i][j] << " ";
            outfile << SE[0][i][j] << " " << SE[1][i][j] << " ";
            outfile << SE[2][i][j] << " " << SE[3][i][j] << "\n";
        }
    }
    outfile.close();

    
    // ********************************* Computation *********************************

    vector<vector<vector<double> > > U0(4, vector<vector<double> >(nx_cell+4, vector<double>(ny_cell+4)));
    vector<vector<vector<double> > > U1(4, vector<vector<double> >(nx_cell+4, vector<double>(ny_cell+4)));
    vector<vector<vector<double> > > V0(4, vector<vector<double> >(nx_cell+4, vector<double>(ny_cell+4)));
    vector<vector<vector<double> > > V1(4, vector<vector<double> >(nx_cell+4, vector<double>(ny_cell+4)));
    vector<vector<vector<double> > > Fx(4, vector<vector<double> >(nx_cell+1, vector<double>(ny_cell)));
    vector<vector<vector<double> > > Fy(4, vector<vector<double> >(nx_cell, vector<double>(ny_cell+1)));
    vector<vector<vector<double> > > phiplus_x(4, vector<vector<double> >(nx_cell+2, vector<double>(ny_cell)));
    vector<vector<vector<double> > > phiminus_x(4, vector<vector<double> >(nx_cell+2, vector<double>(ny_cell)));
    vector<vector<vector<double> > > phiplus_y(4, vector<vector<double> >(nx_cell, vector<double>(ny_cell+2)));
    vector<vector<vector<double> > > phiminus_y(4, vector<vector<double> >(nx_cell, vector<double>(ny_cell+2)));
    vector<vector<vector<double> > > VL_x(4, vector<vector<double> >(nx_cell+1, vector<double>(ny_cell)));
    vector<vector<vector<double> > > VR_x(4, vector<vector<double> >(nx_cell+1, vector<double>(ny_cell)));
    vector<vector<vector<double> > > VL_y(4, vector<vector<double> >(nx_cell, vector<double>(ny_cell+1)));
    vector<vector<vector<double> > > VR_y(4, vector<vector<double> >(nx_cell, vector<double>(ny_cell+1)));
    vector<vector<vector<double> > > R0(4, vector<vector<double> >(nx_cell, vector<double>(ny_cell)));
    vector<vector<double> > t0(nx_cell, vector<double>(ny_cell));

    double gamma;
    double epsilon, kappa;
    bool islimiter;

    gamma = 1.4;
    epsilon = 1.0; // 0 for 1st order and 1 for 2nd order
    kappa = -1.0; // -1 for fully upwind and 0 for upwind biased
    islimiter = true;

    // ******** Initial Conditions ********
    getInitCons(gamma, U0, V0); // calculate initial primitive & conservative variables; Have checked correctness!
    setBounCond(gamma, V0, VLeft, VRight, VBottom, VTop); // set boundary conditions; Have checked correctness!
    primToLimiter(V0, phiplus_x, phiminus_x, phiplus_y, phiminus_y, islimiter); // calcuate limiter phi; Have checked correctness!
    getLeftRightStates(V0, phiplus_x, phiminus_x, VL_x, VR_x, phiplus_y, phiminus_y, VL_y, VR_y, kappa, epsilon); // calculate left and right states for primitive variables, Have checked correctness!
    primToFlux(VL_x, VR_x, VL_y, VR_y, Fx, Fy, x_node, y_node, gamma); // left and right primitives to flux

    // ******** Solve Equations Iteratively (Runge-Kutta 2 stages) ********
    int num_iter;;
    double CFL, delta_t;
    double alpha1, alpha2;
    num_iter = 5000;
    CFL = 1.0;
    alpha1 = 0.5;
    alpha2 = 1.0;
    vector<vector<double> > resi_norm(num_iter, vector<double>(4, 0.0));

    getResidual(Fx, Fy, R0, Area_x, Area_y, SE, Vol);

    for (int i = 0; i < num_iter; i++){
        if (i < 1000) {
            islimiter = true;
        } else {
            islimiter = true;
        }

        getResidual(Fx, Fy, R0, Area_x, Area_y, SE, Vol);
        getL2Norm(R0, resi_norm[i]);
        cout << i << " " << resi_norm[i][0]/resi_norm[0][0] << " " << resi_norm[i][1]/resi_norm[0][1] << " " << resi_norm[i][2]/resi_norm[0][2] << " " << resi_norm[i][3]/resi_norm[0][3] << endl;
        if (resi_norm[i][0]/resi_norm[0][0] < EPSILON && resi_norm[i][1]/resi_norm[0][1] < EPSILON && resi_norm[i][2]/resi_norm[0][2] < EPSILON && resi_norm[i][3]/resi_norm[0][3] < EPSILON) {
            cout << resi_norm[i][0] << " " << resi_norm[i][1] << " " << resi_norm[i][2] << " " << resi_norm[i][3] << endl;
            cout << resi_norm[i][0]/resi_norm[0][0] << " " << resi_norm[i][1]/resi_norm[0][1] << " " << resi_norm[i][2]/resi_norm[0][2] << " " << resi_norm[i][3]/resi_norm[0][3] << endl;
            cout << "The number of steps for convergence is " << i << endl;
            break;
        }
        // calculate delta t
        getDeltaT(V0, Vol, x_node, y_node, t0, delta_t, CFL, gamma);
        // // do update
        // doUpdate(U0, U1, R0, Vol, delta_t);
        // consToPrim(gamma, U0, V0);
        // setBounCond(gamma, V0, VLeft, VRight, VBottom, VTop);
        // primToLimiter(V0, phiplus_x, phiminus_x, phiplus_y, phiminus_y); 
        // getLeftRightStates(V0, phiplus_x, phiminus_x, VL_x, VR_x, phiplus_y, phiminus_y, VL_y, VR_y, kappa, epsilon);
        // primToFlux(VL_x, VR_x, VL_y, VR_y, Fx, Fy, x_node, y_node, gamma); // update flux (Fx & Fy)

        // do update step 1
        doUpdate1(U0, U1, R0, Vol, alpha1, delta_t);
        consToPrim(gamma, U1, V1);
        setBounCond(gamma, V1, VLeft, VRight, VBottom, VTop); // get primitive at ghost cell (including limit range of primitive variables and update conservative)
        primToLimiter(V1, phiplus_x, phiminus_x, phiplus_y, phiminus_y, islimiter); 
        getLeftRightStates(V1, phiplus_x, phiminus_x, VL_x, VR_x, phiplus_y, phiminus_y, VL_y, VR_y, kappa, epsilon);
        primToFlux(VL_x, VR_x, VL_y, VR_y, Fx, Fy, x_node, y_node, gamma); // update flux (Fx & Fy)
        getResidual(Fx, Fy, R0, Area_x, Area_y, SE, Vol); // calculate R(1)
        // do update step 2
        doUpdate2(U0, U1, R0, Vol, alpha2, delta_t);
        consToPrim(gamma, U0, V0);
        setBounCond(gamma, V0, VLeft, VRight, VBottom, VTop);
        primToLimiter(V0, phiplus_x, phiminus_x, phiplus_y, phiminus_y, islimiter); 
        getLeftRightStates(V0, phiplus_x, phiminus_x, VL_x, VR_x, phiplus_y, phiminus_y, VL_y, VR_y, kappa, epsilon);
        primToFlux(VL_x, VR_x, VL_y, VR_y, Fx, Fy, x_node, y_node, gamma); // update flux (Fx & Fy)
        
    }
    consToPrim(gamma,U0,V0);

    ofstream outFile("supersonic-vanleer-result/My-solution-"+cell_name+".dat");
    outFile << "TITLE = \"Manufactured Solution\"\n";
    outFile << "variables=\"x(m)\" \"y(m)\" \"rho(kg/m^3)\" \"u(m/s)\" \"v(m/s)\" \"press(N/m^2)\"\n";

    outFile << "zone T=\"" << 1 << "\"\n";
    outFile << "I=" << nx_cell << " J=" << ny_cell << "\n";
    outFile << "DATAPACKING=POINT\n";

    for (int j = 0; j < x_cell[0].size(); j++) {
        for (int i = 0; i < x_cell.size(); i++) {
            outFile << x_cell[i][j] << " " << y_cell[i][j] << " ";
            outFile << V0[0][i+2][j+2] << " " << V0[1][i+2][j+2] << " ";
            outFile << V0[2][i+2][j+2] << " " << V0[3][i+2][j+2] << "\n";
        }
    }
    outFile.close();

    string fileResi = "residual_" + cell_name + ".txt";
    ofstream writeResi("supersonic-vanleer-result/"+fileResi); // write residual with iteration
    for (const auto& row : resi_norm) {
        for (const auto& element : row) {
            writeResi << element << ' ';
        }
        writeResi << '\n';
    }
    writeResi.close();
    
    return 0;
}

// Function: To get initial conservative variables (rho, rho*u, rho*et)
void getInitCons(double gamma, vector<vector<vector<double> > >& U0, vector<vector<vector<double> > >& V0)
{    
    double rho0, uvel0, vvel0, press0;
    rho0 = 1.0;
    uvel0 = 800.0;
    vvel0 = 800.0;
    press0 = 1.0*1e5;
    
    for (int i = 0; i < U0[0].size(); i++){
        for (int j = 0; j < U0[0][i].size(); j++){
            // primitive variable
            V0[0][i][j] = rho0;
            V0[1][i][j] = uvel0;
            V0[2][i][j] = vvel0;
            V0[3][i][j] = press0;
            // conservative variable
            U0[0][i][j] = rho0;
            U0[1][i][j] = rho0 * uvel0;
            U0[2][i][j] = rho0 * vvel0;
            U0[3][i][j] = 1.0/(gamma-1.0)*press0 + 0.5*rho0*(pow(uvel0,2.0)+pow(vvel0,2.0));
        }
    }
}

// Function: To convert conservative variables to primitive variables
void consToPrim(double gamma, vector<vector<vector<double> > >& U0, vector<vector<vector<double> > >& V0)
{
    for (int i = 0; i < V0[0].size(); i++){
        for (int j = 0; j < V0[0][i].size(); j++){
            V0[0][i][j] = U0[0][i][j];
            V0[1][i][j] = U0[1][i][j]/(U0[0][i][j]+1e-16);
            V0[2][i][j] = U0[2][i][j]/(U0[0][i][j]+1e-16);
            V0[3][i][j] = (U0[3][i][j]-0.5*(pow(U0[1][i][j],2.0)+pow(U0[2][i][j],2.0))/(U0[0][i][j]+1e-16))*(gamma-1.0);
        }
    }
}

// Function: To set boundary conditions to get ghost cell values
void setBounCond(double gamma, vector<vector<vector<double> > >& V0, vector<vector<double> >& VLeft, vector<vector<double> >& VRight, vector<vector<double> >& VBottom, vector<vector<double> >& VTop){
    int Nx, Ny;
    Nx = V0[0].size();
    Ny = V0[0][0].size();
    // left and right boundary
    for (int d = 0; d < V0.size(); d++){
        for (int j = 2; j < V0[0][0].size()-2; j++){
            V0[d][1][j] = 2.0*VLeft[d][j-2] - V0[d][2][j];
            V0[d][0][j] = 2.0*V0[d][1][j] - V0[d][2][j];
            V0[d][Nx-2][j] = 2.0*V0[d][Nx-3][j] - V0[d][Nx-4][j];
            V0[d][Nx-1][j] = 2.0*V0[d][Nx-2][j] - V0[d][Nx-3][j];
        }
    }
    // lower and upper boundary
    for (int d = 0; d < V0.size(); d++){
        for (int i = 2; i < V0[0].size()-2; i++){
            V0[d][i][1] = 2.0*VBottom[d][i-2] - V0[d][i][2];
            V0[d][i][0] = 2.0*V0[d][i][1] - V0[d][i][2];
            V0[d][i][Ny-2] = 2.0*V0[d][i][Ny-3] - V0[d][i][Ny-4];
            V0[d][i][Ny-1] = 2.0*V0[d][i][Ny-2] - V0[d][i][Ny-3];
        }
    }
    // set global limitation for primitive variable
    for (int i = 0; i < V0[0].size(); i++){
        for (int j = 0; j < V0[0][i].size(); j++){
            V0[0][i][j] = max(V0[0][i][j],0.0001);
            V0[3][i][j] = max(V0[3][i][j],500.0);
        }
    }
    // // update conservative variable
    // for (int i = 0; i < U0[0].size(); i++){
    //     for (int j = 0; j < U0[0][i].size(); j++){
    //         U0[0][i][j] = V0[0][i][j];
    //         U0[1][i][j] = V0[0][i][j] * V0[1][i][j];
    //         U0[2][i][j] = V0[0][i][j] * V0[2][i][j];
    //         U0[3][i][j] = 1.0/(gamma-1.0)*V0[3][i][j] + 0.5*V0[0][i][j]*(pow(V0[1][i][j],2.0)+pow(V0[2][i][j],2.0));
    //     }
    // }
}

// Function: To set limitation for rho and P (globally)
void limitGlobal(vector<vector<vector<double> > >& V0){
    for (int i = 0; i < V0[0].size(); i++){
        for (int j = 0; j < V0[0][i].size(); j++){
            V0[0][i][j] = max(V0[0][i][j],0.0001);
            V0[3][i][j] = max(V0[3][i][j],500.0);
        }
    }
}

// Function: To calculate the limiting function phiplus and phiminus from primitive variables
void primToLimiter(vector<vector<vector<double> > >& V0, vector<vector<vector<double> > >& phiplus_x, vector<vector<vector<double> > >& phiminus_x,
vector<vector<vector<double> > >& phiplus_y, vector<vector<vector<double> > >& phiminus_y, bool islimiter)
{
    vector<double> rplus_x(4,0.0), rminus_x(4,0.0), rplus_y(4,0.0), rminus_y(4,0.0);
    if (islimiter == true) {
        for (int d = 0; d < 4; d++){
            for (int i = 0; i < phiplus_x[d].size(); i++){
                for (int j = 0; j < phiplus_x[d][i].size(); j++){
                    rplus_x[d] = (V0[d][i+2][j+2]-V0[d][i+1][j+2])/(copysign(1.0,V0[d][i+1][j+2]-V0[d][i][j+2])*max(fabs(V0[d][i+1][j+2]-V0[d][i][j+2]),1e-6));
                    phiplus_x[d][i][j] = (rplus_x[d]+fabs(rplus_x[d]))/(1.0+rplus_x[d]+1e-16);

                    rminus_x[d] = (V0[d][i+1][j+2]-V0[d][i][j+2])/(copysign(1.0,V0[d][i+2][j+2]-V0[d][i+1][j+2])*max(fabs(V0[d][i+2][j+2]-V0[d][i+1][j+2]),1e-6));
                    phiminus_x[d][i][j] = (rminus_x[d]+fabs(rminus_x[d]))/(1.0+rminus_x[d]+1e-16);
                }
            }
        }
        for (int d = 0; d < 4; d++){
            for (int i = 0; i < phiplus_y[d].size(); i++){
                for (int j = 0; j < phiplus_y[d][i].size(); j++){
                    rplus_y[d] = (V0[d][i+2][j+2]-V0[d][i+2][j+1])/(copysign(1.0,V0[d][i+2][j+1]-V0[d][i+2][j])*max(fabs(V0[d][i+2][j+1]-V0[d][i+2][j]),1e-6));
                    phiplus_y[d][i][j] = (rplus_y[d]+fabs(rplus_y[d]))/(1.0+rplus_y[d]+1e-16);

                    rminus_y[d] = (V0[d][i+2][j+1]-V0[d][i+2][j])/(copysign(1.0,V0[d][i+2][j+2]-V0[d][i+2][j+1])*max(fabs(V0[d][i+2][j+2]-V0[d][i+2][j+1]),1e-6));
                    phiminus_y[d][i][j] = (rminus_y[d]+fabs(rminus_y[d]))/(1.0+rminus_y[d]+1e-16);
                }
            }
        }
    } else {
        for (int d = 0; d < 4; d++){
            for (int i = 0; i < phiplus_x[d].size(); i++){
                for (int j = 0; j < phiplus_x[d][i].size(); j++){
                    phiplus_x[d][i][j] = 1.0;
                    phiminus_x[d][i][j] = 1.0;
                }
            }
        }
        for (int d = 0; d < 4; d++){
            for (int i = 0; i < phiplus_y[d].size(); i++){
                for (int j = 0; j < phiplus_y[d][i].size(); j++){
                    phiplus_y[d][i][j] = 1.0;
                    phiminus_y[d][i][j] = 1.0;
                }
            }
        }
    }
}

// Function: To calculate left and right states 
void getLeftRightStates(vector<vector<vector<double> > >& V0, vector<vector<vector<double> > >& phiplus_x, vector<vector<vector<double> > >& phiminus_x, 
vector<vector<vector<double> > >& VL_x, vector<vector<vector<double> > >& VR_x,
vector<vector<vector<double> > >& phiplus_y, vector<vector<vector<double> > >& phiminus_y, 
vector<vector<vector<double> > >& VL_y, vector<vector<vector<double> > >& VR_y, double kappa, double epsilon)
{
    for (int d = 0; d < 4; d++){
        for (int i = 0; i < VL_x[d].size(); i++){
            for (int j = 0; j < VL_x[d][i].size(); j++){
                VL_x[d][i][j] = V0[d][i+1][j+2] + epsilon/4.0*((1.0-kappa)*phiplus_x[d][i][j]*(V0[d][i+1][j+2]-V0[d][i][j+2])+(1.0+kappa)*phiminus_x[d][i][j]*(V0[d][i+2][j+2]-V0[d][i+1][j+2]));
                VR_x[d][i][j] = V0[d][i+2][j+2] - epsilon/4.0*((1.0-kappa)*phiminus_x[d][i+1][j]*(V0[d][i+3][j+2]-V0[d][i+2][j+2])+(1.0+kappa)*phiplus_x[d][i+1][j]*(V0[d][i+2][j+2]-V0[d][i+1][j+2]));
            }
        }
    }

    for (int d = 0; d < 4; d++){
        for (int i = 0; i < VL_y[d].size(); i++){
            for (int j = 0; j < VL_y[d][i].size(); j++){
                VL_y[d][i][j] = V0[d][i+2][j+1] + epsilon/4.0*((1.0-kappa)*phiplus_y[d][i][j]*(V0[d][i+2][j+1]-V0[d][i+2][j])+(1.0+kappa)*phiminus_y[d][i][j]*(V0[d][i+2][j+2]-V0[d][i+2][j+1]));
                VR_y[d][i][j] = V0[d][i+2][j+2] - epsilon/4.0*((1.0-kappa)*phiminus_y[d][i][j+1]*(V0[d][i+2][j+3]-V0[d][i+2][j+2])+(1.0+kappa)*phiplus_y[d][i][j+1]*(V0[d][i+2][j+2]-V0[d][i+2][j+1]));
            }
        }
    }
}

// Function: To convert left and right states to flux
void primToFlux(vector<vector<vector<double> > >& VL_x, vector<vector<vector<double> > >& VR_x, vector<vector<vector<double> > >& VL_y, vector<vector<vector<double> > >& VR_y, 
vector<vector<vector<double> > >& Fx, vector<vector<vector<double> > >& Fy, vector<vector<double> >& x_node, vector<vector<double> >& y_node, double gamma)
{
    double aL, aR, htL, htR, UhatL, UhatR, ML, MR, Mplus, Mminus, betaL, betaR, alphaplus, alphaminus, cplus, cminus, pbarplus, pbarminus, dplus, dminus;
    double nx, ny, A;
    // calculate Fx
    for (int i = 0; i < Fx[0].size(); i++){
        for (int j = 0; j < Fx[0][i].size(); j++){
            aL = sqrt(gamma*VL_x[3][i][j]/(VL_x[0][i][j]+1e-16));
            aR = sqrt(gamma*VR_x[3][i][j]/(VR_x[0][i][j]+1e-16));
            htL = gamma/(gamma-1.0)*VL_x[3][i][j]/(VL_x[0][i][j]+1e-16) + 0.5*(pow(VL_x[1][i][j],2.0)+pow(VL_x[2][i][j],2.0));
            htR = gamma/(gamma-1.0)*VR_x[3][i][j]/(VR_x[0][i][j]+1e-16) + 0.5*(pow(VR_x[1][i][j],2.0)+pow(VR_x[2][i][j],2.0));
            A = sqrt(pow(x_node[i][j+1]-x_node[i][j],2.0)+pow(y_node[i][j+1]-y_node[i][j],2.0));
            nx = (y_node[i][j+1]-y_node[i][j])/(A+1e-16);
            ny = -1.0*(x_node[i][j+1]-x_node[i][j])/(A+1e-16);
            UhatL = VL_x[1][i][j]*nx + VL_x[2][i][j]*ny;
            UhatR = VR_x[1][i][j]*nx + VR_x[2][i][j]*ny;
            ML = UhatL/(aL+1e-16);
            MR = UhatR/(aR+1e-16);
            Mplus = 0.25*pow((ML+1.0),2.0);
            Mminus = -0.25*pow((MR-1.0),2.0);
            betaL = -1.0 * max(0.0, 1.0-trunc(fabs(ML)));
            betaR = -1.0 * max(0.0, 1.0-trunc(fabs(MR)));
            alphaplus = 0.5*(1.0+copysign(1.0,ML));
            alphaminus = 0.5*(1.0-copysign(1.0,MR));
            cplus = alphaplus*(1.0+betaL)*ML - betaL*Mplus;
            cminus = alphaminus*(1.0+betaR)*MR - betaR*Mminus;
            pbarplus = Mplus*(-1.0*ML+2.0);
            pbarminus = Mminus*(-1.0*MR-2.0);
            dplus = alphaplus*(1+betaL) - betaL*pbarplus;
            dminus = alphaminus*(1+betaR) - betaR*pbarminus;

            Fx[0][i][j] = VL_x[0][i][j]*aL*cplus + VR_x[0][i][j]*aR*cminus;
            Fx[1][i][j] = VL_x[0][i][j]*aL*cplus*VL_x[1][i][j] + VR_x[0][i][j]*aR*cminus*VR_x[1][i][j] + dplus*nx*VL_x[3][i][j] + dminus*nx*VR_x[3][i][j];
            Fx[2][i][j] = VL_x[0][i][j]*aL*cplus*VL_x[2][i][j] + VR_x[0][i][j]*aR*cminus*VR_x[2][i][j] + dplus*ny*VL_x[3][i][j] + dminus*ny*VR_x[3][i][j];
            Fx[3][i][j] = VL_x[0][i][j]*aL*cplus*htL + VR_x[0][i][j]*aR*cminus*htR;
        }
    }
    // calculate Fy
    for (int i = 0; i < Fy[0].size(); i++){
        for (int j = 0; j < Fy[0][i].size(); j++){
            aL = sqrt(gamma*VL_y[3][i][j]/(VL_y[0][i][j]+1e-16));
            aR = sqrt(gamma*VR_y[3][i][j]/(VR_y[0][i][j]+1e-16));
            htL = gamma/(gamma-1.0)*VL_y[3][i][j]/(VL_y[0][i][j]+1e-16) + 0.5*(pow(VL_y[1][i][j],2.0)+pow(VL_y[2][i][j],2.0));
            htR = gamma/(gamma-1.0)*VR_y[3][i][j]/(VR_y[0][i][j]+1e-16) + 0.5*(pow(VR_y[1][i][j],2.0)+pow(VR_y[2][i][j],2.0));
            A = sqrt(pow(x_node[i+1][j]-x_node[i][j],2.0)+pow(y_node[i+1][j]-y_node[i][j],2.0));
            nx = -1.0*(y_node[i+1][j]-y_node[i][j])/(A+1e-16);
            ny = (x_node[i+1][j]-x_node[i][j])/(A+1e-16);
            UhatL = VL_y[1][i][j]*nx + VL_y[2][i][j]*ny;
            UhatR = VR_y[1][i][j]*nx + VR_y[2][i][j]*ny;
            ML = UhatL/(aL+1e-16);
            MR = UhatR/(aR+1e-16);
            Mplus = 0.25*pow((ML+1.0),2.0);
            Mminus = -0.25*pow((MR-1.0),2.0);
            betaL = -1.0 * max(0.0, 1.0-trunc(fabs(ML)));
            betaR = -1.0 * max(0.0, 1.0-trunc(fabs(MR)));
            alphaplus = 0.5*(1.0+copysign(1.0,ML));
            alphaminus = 0.5*(1.0-copysign(1.0,MR));
            cplus = alphaplus*(1.0+betaL)*ML - betaL*Mplus;
            cminus = alphaminus*(1.0+betaR)*MR - betaR*Mminus;
            pbarplus = Mplus*(-1.0*ML+2.0);
            pbarminus = Mminus*(-1.0*MR-2.0);
            dplus = alphaplus*(1+betaL) - betaL*pbarplus;
            dminus = alphaminus*(1+betaR) - betaR*pbarminus;

            Fy[0][i][j] = VL_y[0][i][j]*aL*cplus + VR_y[0][i][j]*aR*cminus;
            Fy[1][i][j] = VL_y[0][i][j]*aL*cplus*VL_y[1][i][j] + VR_y[0][i][j]*aR*cminus*VR_y[1][i][j] + dplus*nx*VL_y[3][i][j] + dminus*nx*VR_y[3][i][j];
            Fy[2][i][j] = VL_y[0][i][j]*aL*cplus*VL_y[2][i][j] + VR_y[0][i][j]*aR*cminus*VR_y[2][i][j] + dplus*ny*VL_y[3][i][j] + dminus*ny*VR_y[3][i][j];
            Fy[3][i][j] = VL_y[0][i][j]*aL*cplus*htL + VR_y[0][i][j]*aR*cminus*htR;
        }
    }
}

// Function: To get face area and cell volume
void getGeo(vector<vector<double> >& x_node, vector<vector<double> >& y_node, vector<vector<double> >& x_cell, vector<vector<double> >& y_cell, 
vector<vector<double> >& Area_x, vector<vector<double> >& Area_y, vector<vector<double> >& Vol)
{
    double dx1, dy1, dx2, dy2;
    // Area at | | | |
    for (int i = 0; i < Area_x.size(); i++){
        for (int j = 0; j < Area_x[i].size(); j++){
            Area_x[i][j] = sqrt(pow(x_node[i][j+1]-x_node[i][j],2.0)+pow(y_node[i][j+1]-y_node[i][j],2.0));
        }
    }
    // Area at _ _ _
    for (int i = 0; i < Area_y.size(); i++){
        for (int j = 0; j < Area_y[i].size(); j++){
            Area_y[i][j] = sqrt(pow(x_node[i+1][j]-x_node[i][j],2.0)+pow(y_node[i+1][j]-y_node[i][j],2.0));
        }
    }
    // cell volume
    for (int i = 0; i < Vol.size(); i++){
        for (int j = 0; j < Vol[i].size(); j++){
            dx1 = x_node[i+1][j+1]-x_node[i][j];
            dy1 = y_node[i+1][j+1]-y_node[i][j];
            dx2 = x_node[i][j+1]-x_node[i+1][j];
            dy2 = y_node[i][j+1]-y_node[i+1][j];
            Vol[i][j] = 0.5*fabs(dx1*dy2 - dx2*dy1);   
        }
    }
    // cell center coordinates
    for (int i = 0; i < x_cell.size(); i++){
        for (int j = 0; j < x_cell[i].size(); j++){
            x_cell[i][j] = 0.25 * (x_node[i][j] + x_node[i][j+1] + x_node[i+1][j] + x_node[i+1][j+1]);
            y_cell[i][j] = 0.25 * (y_node[i][j] + y_node[i][j+1] + y_node[i+1][j] + y_node[i+1][j+1]);
        }
    }
}

// Function: To get residual
void getResidual(vector<vector<vector<double> > >& Fx, vector<vector<vector<double> > >& Fy, vector<vector<vector<double> > >& R0, vector<vector<double> >& Area_x, vector<vector<double> >& Area_y,
vector<vector<vector<double> > >& SE, vector<vector<double> >& Vol)
{
    for (int d = 0; d < R0.size(); d++){
        for (int i = 0; i < R0[d].size(); i++){
            for (int j = 0; j < R0[d][i].size(); j++){
                R0[d][i][j] = -1.0*Fx[d][i][j]*Area_x[i][j] + Fx[d][i+1][j]*Area_x[i+1][j] - 1.0*Fy[d][i][j]*Area_y[i][j] + Fy[d][i][j+1]*Area_y[i][j+1] - SE[d][i][j]*Vol[i][j];
            }
        }
    }
}

// // Function: To get the L1 norm of residual
// void getL1Norm(vector<vector<vector<double> > >& R0, vector<double>& current_resi_norm)
// {
//     vector<double> norm(4, 0.0);
//     for (int d = 0; d < R0.size(); d++){
//         for (int i = 0; i < R0[d].size(); i++){
//             for (int j = 0; j < R0[d][i].size(); j++){
//                 norm[d] += fabs(R0[d][i][j]);
//             }
//         }
//     }
//     for (int i = 0; i < current_resi_norm.size(); i++){
//         current_resi_norm[i] = norm[i]/(R0[0].size()*R0[0][0].size());
//     }
// }

// Function: To get the L2 norm of residual
void getL2Norm(vector<vector<vector<double> > >& R0, vector<double>& current_resi_norm) {
    vector<double> norm(4, 0.0);
    for (int d = 0; d < R0.size(); d++){
        for (int i = 0; i < R0[d].size(); i++){
            for (int j = 0; j < R0[d][i].size(); j++){
                norm[d] += pow(R0[d][i][j],2.0);
            }
        }
    }
    for (int i = 0; i < current_resi_norm.size(); i++){
        current_resi_norm[i] = sqrt(norm[i]/double(R0[0].size()*R0[0][0].size()));
    }
}

// Function: To calculate exact solutions and source terms
void getManufSol(vector<vector<vector<double> > >& VE, vector<vector<vector<double> > >& SE, vector<vector<double> >& VLeft, vector<vector<double> >& VRight,
vector<vector<double> >& VTop, vector<vector<double> >& VBottom, vector<vector<double> >& x_cell, vector<vector<double> >& y_cell,
vector<vector<double> >& x_node, vector<vector<double> >& y_node){
    double Pi, L, gamma;
    double rho0, rhox, rhoy;
    double uvel0, uvelx, uvely, vvel0, vvelx, vvely, wvel0, wvelx, wvely;
    double press0, pressx, pressy;
    double x, y;
    
    Pi = acos(-1.); L = 1.0; gamma = 1.4;
    rho0 = 1.0; rhox = 0.15; rhoy = -0.1;
    uvel0 = 800.0; uvelx = 50.0; uvely = -30.0;
    vvel0 = 800.0; vvelx = -75.0; vvely = 40.0;
    wvel0 = 0.0; wvelx = 0.0; wvely = 0.0;
    press0 = 1.0*1e5; pressx = 0.2*1e5; pressy = 0.5*1e5;
    // cell centers
    for (int i = 0; i < VE[0].size(); i++){
        for (int j = 0; j < VE[0][i].size(); j++){
            x = x_cell[i][j];
            y = y_cell[i][j];
            // primitive variable at cell centers
            VE[0][i][j] = rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L);
            VE[1][i][j] = uvel0 + uvely*cos((3.*Pi*y)/(5.*L)) + uvelx*sin((3.*Pi*x)/(2.*L));
            VE[2][i][j] = vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2.*Pi*y)/(3.*L));
            VE[3][i][j] = press0 + pressx*cos((2.*Pi*x)/L) + pressy*sin((Pi*y)/L);
            // source term at cell centers
            SE[0][i][j] = (3*Pi*uvelx*cos((3*Pi*x)/(2.*L))*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L)))/(2.*L) + 
                        (2*Pi*vvely*cos((2*Pi*y)/(3.*L))*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L)))/(3.*L) + 
                        (Pi*rhox*cos((Pi*x)/L)*(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L))))/L - 
                        (Pi*rhoy*sin((Pi*y)/(2.*L))*(vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/(3.*L))))/(2.*L);
            SE[1][i][j] = (3*Pi*uvelx*cos((3*Pi*x)/(2.*L))*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))*
                        (uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L))))/L + (2*Pi*vvely*cos((2*Pi*y)/(3.*L))*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + 
                        rhox*sin((Pi*x)/L))*(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L))))/(3.*L) + (Pi*rhox*cos((Pi*x)/L)*pow(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + 
                        uvelx*sin((3*Pi*x)/(2.*L)),2.)/L - 2*Pi*pressx*sin((2*Pi*x)/L))/L - (Pi*rhoy*(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + 
                        uvelx*sin((3*Pi*x)/(2.*L)))*sin((Pi*y)/(2.*L))*(vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/(3.*L))))/(2.*L) - 
                        (3*Pi*uvely*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))*sin((3*Pi*y)/(5.*L))*(vvel0 + vvelx*cos((Pi*x)/(2.*L)) + 
                        vvely*sin((2*Pi*y)/(3.*L))))/(5.*L);
            SE[2][i][j] = (Pi*pressy*cos((Pi*y)/L))/L - (Pi*vvelx*sin((Pi*x)/(2.*L))*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))*(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + 
                        uvelx*sin((3*Pi*x)/(2.*L))))/(2.*L) + (3*Pi*uvelx*cos((3*Pi*x)/(2.*L))*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))*(vvel0 + vvelx*cos((Pi*x)/(2.*L)) + 
                        vvely*sin((2*Pi*y)/(3.*L))))/(2.*L) + (4*Pi*vvely*cos((2*Pi*y)/(3.*L))*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))*(vvel0 + vvelx*cos((Pi*x)/(2.*L)) + 
                        vvely*sin((2*Pi*y)/(3.*L))))/(3.*L) + (Pi*rhox*cos((Pi*x)/L)*(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L)))*(vvel0 + vvelx*cos((Pi*x)/(2.*L)) + 
                        vvely*sin((2*Pi*y)/(3.*L))))/L - Pi*rhoy*sin((Pi*y)/(2.*L))*pow(vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/(3.*L)),2.0)/(2.*L);
            SE[3][i][j] = (uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L)))*
                        ((-2*Pi*pressx*sin((2*Pi*x)/L))/L + (rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))*
                        ((-2*Pi*pressx*sin((2*Pi*x)/L))/((-1 + gamma)*L*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))) + 
                        ((3*Pi*uvelx*cos((3*Pi*x)/(2.*L))*(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L))))/L - 
                        (Pi*vvelx*sin((Pi*x)/(2.*L))*(vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/(3.*L))))/L)/2. - 
                        (Pi*rhox*cos((Pi*x)/L)*(press0 + pressx*cos((2*Pi*x)/L) + pressy*sin((Pi*y)/L)))/
                        ((-1 + gamma)*L*pow(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L),2.))) + 
                        (Pi*rhox*cos((Pi*x)/L)*((pow(wvel0,2.) + pow(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L)),2.) + 
                        pow(vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/(3.*L)),2.))/2. + 
                        (press0 + pressx*cos((2*Pi*x)/L) + pressy*sin((Pi*y)/L))/
                        ((-1 + gamma)*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L)))))/L) + 
                        (3*Pi*uvelx*cos((3*Pi*x)/(2.*L))*(press0 + pressx*cos((2*Pi*x)/L) + pressy*sin((Pi*y)/L) + 
                        (rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))*
                        ((pow(wvel0,2.) + pow(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L)),2.) + 
                        pow(vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/(3.*L)),2.))/2. + 
                        (press0 + pressx*cos((2*Pi*x)/L) + pressy*sin((Pi*y)/L))/
                        ((-1 + gamma)*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))))))/(2.*L) + 
                        (2*Pi*vvely*cos((2*Pi*y)/(3.*L))*(press0 + pressx*cos((2*Pi*x)/L) + pressy*sin((Pi*y)/L) + 
                        (rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))*
                        ((pow(wvel0,2.) + pow(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L)),2.) + 
                        pow(vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/(3.*L)),2))/2. + 
                        (press0 + pressx*cos((2*Pi*x)/L) + pressy*sin((Pi*y)/L))/
                        ((-1 + gamma)*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))))))/(3.*L) + 
                        (vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/(3.*L)))*
                        ((Pi*pressy*cos((Pi*y)/L))/L - (Pi*rhoy*sin((Pi*y)/(2.*L))*
                        ((pow(wvel0,2.) + pow(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L)),2.) + 
                        pow(vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/(3.*L)),2.))/2. + 
                        (press0 + pressx*cos((2*Pi*x)/L) + pressy*sin((Pi*y)/L))/
                        ((-1 + gamma)*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L)))))/(2.*L) + 
                        (rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))*
                        ((Pi*pressy*cos((Pi*y)/L))/((-1 + gamma)*L*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))) + 
                        ((-6*Pi*uvely*(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L)))*sin((3*Pi*y)/(5.*L)))/(5.*L) + 
                        (4*Pi*vvely*cos((2*Pi*y)/(3.*L))*(vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/(3.*L))))/(3.*L))/2. + 
                        (Pi*rhoy*sin((Pi*y)/(2.*L))*(press0 + pressx*cos((2*Pi*x)/L) + pressy*sin((Pi*y)/L)))/
                        (2.*(-1 + gamma)*L*pow(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L),2.))));
        }
    }
    // primitive variables at boundary face centers
    int Nx_node, Ny_node;
    Nx_node = x_node.size();
    Ny_node = x_node[0].size();
    // left boundary
    for (int j = 0; j < VLeft[0].size(); j++){
        x = 0.5 * (x_node[0][j] + x_node[0][j+1]);
        y = 0.5 * (y_node[0][j] + y_node[0][j+1]);
        VLeft[0][j] = rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L);
        VLeft[1][j] = uvel0 + uvely*cos((3.*Pi*y)/(5.*L)) + uvelx*sin((3.*Pi*x)/(2.*L));
        VLeft[2][j] = vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2.*Pi*y)/(3.*L));
        VLeft[3][j] = press0 + pressx*cos((2.*Pi*x)/L) + pressy*sin((Pi*y)/L);
    }
    // right boundary
    for (int j = 0; j < VRight[0].size(); j++){
        x = 0.5 * (x_node[Nx_node-1][j] + x_node[Nx_node-1][j+1]);
        y = 0.5 * (y_node[Nx_node-1][j] + y_node[Nx_node-1][j+1]);
        VRight[0][j] = rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L);
        VRight[1][j] = uvel0 + uvely*cos((3.*Pi*y)/(5.*L)) + uvelx*sin((3.*Pi*x)/(2.*L));
        VRight[2][j] = vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2.*Pi*y)/(3.*L));
        VRight[3][j] = press0 + pressx*cos((2.*Pi*x)/L) + pressy*sin((Pi*y)/L);
    }
    // bottom boundary
    for (int i = 0; i < VBottom[0].size(); i++){
        x = 0.5 * (x_node[i][0] + x_node[i+1][0]);
        y = 0.5 * (y_node[i][0] + y_node[i+1][0]);
        VBottom[0][i] = rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L);
        VBottom[1][i] = uvel0 + uvely*cos((3.*Pi*y)/(5.*L)) + uvelx*sin((3.*Pi*x)/(2.*L));
        VBottom[2][i] = vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2.*Pi*y)/(3.*L));
        VBottom[3][i] = press0 + pressx*cos((2.*Pi*x)/L) + pressy*sin((Pi*y)/L);
    }
    // top boundary
    for (int i = 0; i < VTop[0].size(); i++){
        x = 0.5 * (x_node[i][Ny_node-1] + x_node[i+1][Ny_node-1]);
        y = 0.5 * (y_node[i][Ny_node-1] + y_node[i+1][Ny_node-1]);
        VTop[0][i] = rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L);
        VTop[1][i] = uvel0 + uvely*cos((3.*Pi*y)/(5.*L)) + uvelx*sin((3.*Pi*x)/(2.*L));
        VTop[2][i] = vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2.*Pi*y)/(3.*L));
        VTop[3][i] = press0 + pressx*cos((2.*Pi*x)/L) + pressy*sin((Pi*y)/L);
    }
}

// Function: To do iteration with Runge-Kutta 2 stages
void doUpdate1(vector<vector<vector<double> > >& U0, vector<vector<vector<double> > >& U1, vector<vector<vector<double> > >& R0, vector<vector<double> >& Vol, double alpha1, double delta_t){
    for (int d = 0; d < U0.size(); d++){
        for (int i = 2; i < U0[d].size()-2; i++){
            for (int j = 2; j < U0[d][i].size()-2; j++){
                U1[d][i][j] = U0[d][i][j] - alpha1*delta_t/Vol[i-2][j-2]*R0[d][i-2][j-2];
            } 
        }
    }
}

void doUpdate2(vector<vector<vector<double> > >& U0, vector<vector<vector<double> > >& U1, vector<vector<vector<double> > >& R0, vector<vector<double> >& Vol, double alpha2, double delta_t){
    for (int d = 0; d < U0.size(); d++){
        for (int i = 2; i < U0[d].size()-2; i++){
            for (int j = 2; j < U0[d][i].size()-2; j++){
                U1[d][i][j] = U0[d][i][j] - alpha2*delta_t/Vol[i-2][j-2]*R0[d][i-2][j-2];
            } 
        }
    }
    // assign U1 to U0
    for (int d = 0; d < U0.size(); d++){
        for (int i = 0; i < U0[d].size(); i++){
            for (int j = 0; j < U0[d][i].size(); j++){
                U0[d][i][j] = U1[d][i][j]; // U0 = U1
            } 
        }
    }
}

void doUpdate(vector<vector<vector<double> > >& U0, vector<vector<vector<double> > >& U1, vector<vector<vector<double> > >& R0, vector<vector<double> >& Vol, double delta_t){
    for (int d = 0; d < U0.size(); d++){
        for (int i = 2; i < U0[d].size()-2; i++){
            for (int j = 2; j < U0[d][i].size()-2; j++){
                U1[d][i][j] = U0[d][i][j] - delta_t/Vol[i-2][j-2]*R0[d][i-2][j-2];
            } 
        }
    }
    U0 = U1;
}

// Function: To get local time step
void getDeltaT(vector<vector<vector<double> > >& V0, vector<vector<double> >& Vol, vector<vector<double> >& x_node, vector<vector<double> >& y_node, 
vector<vector<double> >& t0, double & delta_t, double CFL, double gamma){
    double nks_x0, nks_y0, nks_x1, nks_y1;
    double neta_x0, neta_y0, neta_x1, neta_y1;
    double nks_x, nks_y, neta_x, neta_y;
    double Aks_0, Aks_1, Aeta_0, Aeta_1;
    double rho, u, v, p, a;
    double lambdaks, lambdaeta;
    for (int i = 0; i < t0.size(); i++){
        for (int j = 0; j < t0[i].size(); j++){
            Aks_0 = sqrt(pow(x_node[i][j+1]-x_node[i][j],2.0)+pow(y_node[i][j+1]-y_node[i][j],2.0));
            Aks_1 = sqrt(pow(x_node[i+1][j+1]-x_node[i+1][j],2.0)+pow(y_node[i+1][j+1]-y_node[i+1][j],2.0));
            Aeta_0 = sqrt(pow(x_node[i+1][j]-x_node[i][j],2.0)+pow(y_node[i+1][j]-y_node[i][j],2.0));
            Aeta_1 = sqrt(pow(x_node[i+1][j+1]-x_node[i][j+1],2.0)+pow(y_node[i+1][j+1]-y_node[i][j+1],2.0));
            nks_x0 = (y_node[i][j+1]-y_node[i][j])/Aks_0;
            nks_y0 = -1.0*(x_node[i][j+1]-x_node[i][j])/Aks_0;
            nks_x1 = (y_node[i+1][j+1]-y_node[i+1][j])/Aks_1;
            nks_y1 = -1.0*(x_node[i+1][j+1]-x_node[i+1][j])/Aks_1;
            neta_x0 = -1.0*(y_node[i+1][j]-y_node[i][j])/Aeta_0;
            neta_y0 = (x_node[i+1][j]-x_node[i][j])/Aeta_0;
            neta_x1 = -1.0*(y_node[i+1][j+1]-y_node[i][j+1])/Aeta_1;
            neta_y1 = (x_node[i+1][j+1]-x_node[i][j+1])/Aeta_1;
            nks_x = 0.5*(nks_x0+nks_x1);
            nks_y = 0.5*(nks_y0+nks_y1);
            neta_x = 0.5*(neta_x0+neta_x1);
            neta_y = 0.5*(neta_y0+neta_y1);
            
            rho = V0[0][i][j];
            u = V0[1][i][j];
            v = V0[2][i][j];
            p = V0[3][i][j];
            a = sqrt(gamma*p/(rho+1e-16)); // speed of sound
            lambdaks = fabs(u*nks_x+v*nks_y)+a;
            lambdaeta = fabs(u*neta_x+v*neta_y)+a;

            t0[i][j] = CFL*Vol[i][j]/(lambdaks*(Aks_0+Aks_1)/2.0+lambdaeta*(Aeta_0+Aeta_1)/2.0+1e-16);
        }
    }
    // get global minimum
    delta_t = 100.0;
    for (int i = 0; i < t0.size(); i++) {
        for (int j = 0; j < t0[i].size(); j++) {
            if (t0[i][j] < delta_t) {
                delta_t = t0[i][j];
            }
        }
    }
}

