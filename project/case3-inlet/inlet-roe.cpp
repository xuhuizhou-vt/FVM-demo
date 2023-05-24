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
void getInitCons(double gamma, vector<vector<vector<double> > >& U0, vector<vector<vector<double> > >& V0);
void setBounCond(double gamma, vector<vector<vector<double> > >& V0, vector<vector<double> >& x_node, vector<vector<double> >& y_node);
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
void getResidual(vector<vector<vector<double> > >& Fx, vector<vector<vector<double> > >& Fy, vector<vector<vector<double> > >& R0, vector<vector<double> >& Area_x, vector<vector<double> >& Area_y);
void getL1Norm(vector<vector<vector<double> > >& R0, vector<double>& current_resi_norm);
void getDeltaT(vector<vector<vector<double> > >& V0, vector<vector<double> >& Vol, vector<vector<double> >& x_node, vector<vector<double> >& y_node, 
vector<vector<double> >& t0, double & delta_t, double CFL, double gamma);
void doUpdate1(vector<vector<vector<double> > >& U0, vector<vector<vector<double> > >& U1, vector<vector<vector<double> > >& R0, vector<vector<double> >& Vol, double alpha1, double delta_t);
void doUpdate2(vector<vector<vector<double> > >& U0, vector<vector<vector<double> > >& U1, vector<vector<vector<double> > >& R0, vector<vector<double> >& Vol, double alpha2, double delta_t);
void doUpdate(vector<vector<vector<double> > >& U0, vector<vector<vector<double> > >& U1, vector<vector<vector<double> > >& R0, vector<vector<double> >& Vol, double delta_t);
void getPressLoss(vector<vector<vector<double> > >& V0, vector<vector<double> >& Area_x, double gamma, double & pressLoss);


int main() {

    // ************************************************ Geometry ************************************************

    // Declare variables
    int nzones, imax, jmax, kmax;
    int nx_cell, ny_cell;
    string cell_name;

    cell_name = "417x129";

    // Open file
    std::ifstream inFile("../Project_Files/Grids/Inlet-Grids/Inlet." + cell_name + ".grd");

    // Read in data
    inFile >> nzones;
    inFile >> imax >> jmax >> kmax; // ith column & jth row !!!
    cout << nzones << " " << imax << " " << jmax << " " << kmax << endl;

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
                inFile >> x_node[i][j];
            }
        }
    }

    for (int k = 0; k < kmax; k++) {
        for (int j = 0; j < jmax; j++) {
            for (int i = 0; i < imax; i++) {
                inFile >> y_node[i][j];
            }
        }
    }

    // cout << x_node[0][0] << " " << x_node[2][0] << " " << x_node[4][0] << " " << x_node[6][0] << " " << x_node[8][0] << endl;
    // cout << y_node[0][0] << " " << y_node[2][0] << " " << y_node[4][0] << " " << y_node[6][0] << " " << y_node[8][0] << endl;
    // cout << x_node[0][8] << " " << x_node[2][8] << " " << x_node[4][8] << " " << x_node[6][8] << " " << x_node[8][8] << endl;
    // cout << y_node[0][8] << " " << y_node[2][8] << " " << y_node[4][8] << " " << y_node[6][8] << " " << y_node[8][8] << endl;

    // Close file
    inFile.close();

    // write manufactured solution to tecplot file
    ofstream outfile("inlet-result/grid" + cell_name + ".dat");
    outfile << "TITLE = \"Grid\"\n";
    outfile << "variables=\"x(m)\" \"y(m)\"\n";

    outfile << "zone T=\"" << 1 << "\"\n";
    outfile << "I=" << imax << " J=" << jmax << "\n";
    outfile << "DATAPACKING=POINT\n";

    for (int j = 0; j < x_node[0].size(); j++) {
        for (int i = 0; i < x_node.size(); i++) {
            outfile << x_node[i][j] << " " << y_node[i][j] << "\n";
        }
    }
    outfile.close();

    // calculate face areas and cell volumes
    getGeo(x_node, y_node, x_cell, y_cell, Area_x, Area_y, Vol); // Have checked correctness!
    
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
    epsilon = 0.0; // 0 for 1st order and 1 for 2nd order
    kappa = -1.0; // -1 for fully upwind and 0 for upwind biased
    islimiter = true;

    // ******** Initial Conditions ********
    getInitCons(gamma, U0, V0); // calculate initial primitive & conservative variables; Have checked correctness!
    setBounCond(gamma, V0, x_node, y_node); // set boundary conditions; Have checked correctness!
    primToLimiter(V0, phiplus_x, phiminus_x, phiplus_y, phiminus_y, islimiter); // calcuate limiter phi; Have checked correctness!
    getLeftRightStates(V0, phiplus_x, phiminus_x, VL_x, VR_x, phiplus_y, phiminus_y, VL_y, VR_y, kappa, epsilon); // calculate left and right states for primitive variables, Have checked correctness!
    primToFlux(VL_x, VR_x, VL_y, VR_y, Fx, Fy, x_node, y_node, gamma); // left and right primitives to flux

    // ******** Solve Equations Iteratively (Runge-Kutta 2 stages) ********
    int num_iter;;
    double CFL, delta_t;
    double alpha1, alpha2;
    double pressLoss;
    num_iter = 500000;
    CFL = 1.2;
    alpha1 = 0.5;
    alpha2 = 1.0;
    pressLoss = 0.0;
    vector<vector<double> > resi_norm(num_iter, vector<double>(4, 0.0));

    // getResidual(Fx, Fy, R0, Area_x, Area_y);

    for (int i = 0; i < num_iter; i++){
        
        if (i < 20000) {
            islimiter = true;
        } else {
            islimiter = true;
        }

        getResidual(Fx, Fy, R0, Area_x, Area_y);
        getL1Norm(R0, resi_norm[i]);
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
        // setBounCond(gamma, V0, x_node, y_node);
        // primToLimiter(V0, phiplus_x, phiminus_x, phiplus_y, phiminus_y); 
        // getLeftRightStates(V0, phiplus_x, phiminus_x, VL_x, VR_x, phiplus_y, phiminus_y, VL_y, VR_y, kappa, epsilon);
        // primToFlux(VL_x, VR_x, VL_y, VR_y, Fx, Fy, x_node, y_node, gamma); // update flux (Fx & Fy)

        // do update step 1
        doUpdate1(U0, U1, R0, Vol, alpha1, delta_t);
        consToPrim(gamma, U1, V1);
        setBounCond(gamma, V1, x_node, y_node); // get primitive at ghost cell (including limit range of primitive variables and update conservative)
        primToLimiter(V1, phiplus_x, phiminus_x, phiplus_y, phiminus_y, islimiter); 
        getLeftRightStates(V1, phiplus_x, phiminus_x, VL_x, VR_x, phiplus_y, phiminus_y, VL_y, VR_y, kappa, epsilon);
        primToFlux(VL_x, VR_x, VL_y, VR_y, Fx, Fy, x_node, y_node, gamma); // update flux (Fx & Fy)
        getResidual(Fx, Fy, R0, Area_x, Area_y); // calculate R(1)
        // do update step 2
        doUpdate2(U0, U1, R0, Vol, alpha2, delta_t);
        consToPrim(gamma, U0, V0);
        setBounCond(gamma, V0, x_node, y_node);
        primToLimiter(V0, phiplus_x, phiminus_x, phiplus_y, phiminus_y, islimiter); 
        getLeftRightStates(V0, phiplus_x, phiminus_x, VL_x, VR_x, phiplus_y, phiminus_y, VL_y, VR_y, kappa, epsilon);
        primToFlux(VL_x, VR_x, VL_y, VR_y, Fx, Fy, x_node, y_node, gamma); // update flux (Fx & Fy)
        
    }
    consToPrim(gamma,U0,V0);
    getPressLoss(V0, Area_x, gamma, pressLoss);
    cout << "Total pressure loss: " << " " << pressLoss << endl;


    ofstream outFile("inlet-result/My-solution" + cell_name + ".dat");
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
    ofstream writeResi("inlet-result/"+fileResi); // write residual with iteration
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
    double Mach, press, temp, R_air;
    double rho, uvel, vvel;
    Mach = 4.0;
    press = 12270.0;
    temp = 217.0;
    R_air = 287.0;
    rho = press/(R_air*temp+1e-16);
    uvel = Mach*sqrt(gamma*R_air*temp);
    vvel = 0.0;

    for (int i = 0; i < U0[0].size(); i++){
        for (int j = 0; j < U0[0][i].size(); j++){
            // primitive variable
            V0[0][i][j] = rho;
            V0[1][i][j] = uvel;
            V0[2][i][j] = vvel;
            V0[3][i][j] = press;
            // conservative variable
            U0[0][i][j] = rho;
            U0[1][i][j] = rho * uvel;
            U0[2][i][j] = rho * vvel;
            U0[3][i][j] = 1.0/(gamma-1.0)*press + 0.5*rho*(pow(uvel,2.0)+pow(vvel,2.0));
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
void setBounCond(double gamma, vector<vector<vector<double> > >& V0, vector<vector<double> >& x_node, vector<vector<double> >& y_node){
    int Nx, Ny, Nxp, Nyp;
    double Mach, press, temp, R_air;
    double rho, uvel, vvel;
    vector<double> s_vec(2);

    Mach = 4.0;
    press = 12270.0;
    temp = 217.0;
    R_air = 287.0;
    rho = press/(R_air*temp+1e-16);
    uvel = Mach*sqrt(gamma*R_air*temp);
    vvel = 0.0;

    Nx = V0[0].size();
    Ny = V0[0][0].size();
    Nxp = x_node.size();
    Nyp = x_node[0].size();

    // top boundary (inlet + upper wall)
    for (int i = 2; i < Nx-2; i++){
        if (x_node[i-1][Nyp-1] < 0.5+1e-5) {
            V0[0][i][Ny-2] = 2.0*rho - V0[0][i][Ny-3];
            V0[0][i][Ny-1] = 2.0*V0[0][i][Ny-2] - V0[0][i][Ny-3];
            V0[1][i][Ny-2] = 2.0*uvel - V0[1][i][Ny-3];
            V0[1][i][Ny-1] = 2.0*V0[1][i][Ny-2] - V0[1][i][Ny-3];
            V0[2][i][Ny-2] = 2.0*vvel - V0[2][i][Ny-3];
            V0[2][i][Ny-1] = 2.0*V0[2][i][Ny-2] - V0[2][i][Ny-3];
            V0[3][i][Ny-2] = 2.0*press - V0[3][i][Ny-3];
            V0[3][i][Ny-1] = 2.0*V0[3][i][Ny-2] - V0[3][i][Ny-3];
        } else {
            s_vec[0] = x_node[i-1][Nyp-1] - x_node[i-2][Nyp-1];
            s_vec[1] = y_node[i-1][Nyp-1] - y_node[i-2][Nyp-1];
            V0[1][i][Ny-2] = (V0[1][i][Ny-3]*(s_vec[0]*s_vec[0]-s_vec[1]*s_vec[1])+2.0*V0[2][i][Ny-3]*s_vec[0]*s_vec[1])/(s_vec[0]*s_vec[0]+s_vec[1]*s_vec[1]+1e-16);
            V0[2][i][Ny-2] = (2.0*V0[1][i][Ny-3]*s_vec[0]*s_vec[1]+V0[2][i][Ny-3]*(s_vec[1]*s_vec[1]-s_vec[0]*s_vec[0]))/(s_vec[0]*s_vec[0]+s_vec[1]*s_vec[1]+1e-16);
            V0[1][i][Ny-1] = 2.0 * V0[1][i][Ny-2] - V0[1][i][Ny-3];
            V0[2][i][Ny-1] = 2.0 * V0[2][i][Ny-2] - V0[2][i][Ny-3];
            V0[3][i][Ny-2] = 2.0 * V0[3][i][Ny-3] - V0[3][i][Ny-4];
            V0[3][i][Ny-1] = 3.0 * V0[3][i][Ny-3] - 2.0 * V0[3][i][Ny-4];
            V0[0][i][Ny-2] = V0[3][i][Ny-2]*V0[0][i][Ny-3]/(V0[3][i][Ny-3]+1e-16);
            V0[0][i][Ny-1] = 2.0 * V0[0][i][Ny-2] - V0[0][i][Ny-3];
        }
    }
    
    // bottom boundary (lower wall - 2)
    for (int i = 2; i < Nx-2; i++){
        s_vec[0] = x_node[i-1][0] - x_node[i-2][0];
        s_vec[1] = y_node[i-1][0] - y_node[i-2][0];
        V0[1][i][1] = (V0[1][i][2]*(s_vec[0]*s_vec[0]-s_vec[1]*s_vec[1])+2.0*V0[2][i][2]*s_vec[0]*s_vec[1])/(s_vec[0]*s_vec[0]+s_vec[1]*s_vec[1]+1e-16);
        V0[2][i][1] = (2.0*V0[1][i][2]*s_vec[0]*s_vec[1]+V0[2][i][2]*(s_vec[1]*s_vec[1]-s_vec[0]*s_vec[0]))/(s_vec[0]*s_vec[0]+s_vec[1]*s_vec[1]+1e-16);
        V0[1][i][0] = 2.0 * V0[1][i][1] - V0[1][i][2];
        V0[2][i][0] = 2.0 * V0[2][i][1] - V0[2][i][2];
        V0[3][i][1] = 2.0 * V0[3][i][2] - V0[3][i][3];
        V0[3][i][0] = 3.0 * V0[3][i][2] - 2.0 * V0[3][i][3];
        V0[0][i][1] = V0[3][i][1]*V0[0][i][2]/(V0[3][i][2]+1e-16);
        V0[0][i][0] = 2.0 * V0[0][i][1] - V0[0][i][2];
    }

    // left boundary (low wall - 1)
    for (int j = 2; j < Ny-2; j++){
        s_vec[0] = x_node[0][j-1] - x_node[0][j-2];
        s_vec[1] = y_node[0][j-1] - y_node[0][j-2];
        V0[1][1][j] = (V0[1][2][j]*(s_vec[0]*s_vec[0]-s_vec[1]*s_vec[1])+2.0*V0[2][2][j]*s_vec[0]*s_vec[1])/(s_vec[0]*s_vec[0]+s_vec[1]*s_vec[1]+1e-16);
        V0[2][1][j] = (2.0*V0[1][2][j]*s_vec[0]*s_vec[1]+V0[2][2][j]*(s_vec[1]*s_vec[1]-s_vec[0]*s_vec[0]))/(s_vec[0]*s_vec[0]+s_vec[1]*s_vec[1]+1e-16);
        V0[1][0][j] = 2.0*V0[1][1][j] - V0[1][2][j];
        V0[2][0][j] = 2.0*V0[2][1][j] - V0[2][2][j];
        V0[3][1][j] = 2.0*V0[3][2][j] - V0[3][3][j];
        V0[3][0][j] = 3.0*V0[3][2][j] - 2.0*V0[3][3][j];
        V0[0][1][j] = V0[3][1][j]*V0[0][2][j]/(V0[3][2][j]+1e-16);
        V0[0][0][j] = 2.0*V0[0][1][j] - V0[0][2][j];
    }

    // right boundary (outlet)
    for (int d = 0; d < V0.size(); d++){
        for (int j = 2; j < V0[0][0].size()-2; j++){
            V0[d][Nx-2][j] = 2.0*V0[d][Nx-3][j] - V0[d][Nx-4][j];
            V0[d][Nx-1][j] = 2.0*V0[d][Nx-2][j] - V0[d][Nx-3][j];
        }
    }

    // set global limitation for primitive variable
    for (int i = 0; i < V0[0].size(); i++){
        for (int j = 0; j < V0[0][i].size(); j++){
            V0[0][i][j] = max(V0[0][i][j],0.001);
            V0[3][i][j] = max(V0[3][i][j],100.0);
        }
    }
}

// Function: To set limitation for rho and P (globally)
void limitGlobal(vector<vector<vector<double> > >& V0){
    for (int i = 0; i < V0[0].size(); i++){
        for (int j = 0; j < V0[0][i].size(); j++){
            V0[0][i][j] = max(V0[0][i][j],0.001);
            V0[3][i][j] = max(V0[3][i][j],100.0);
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

// Function: To convert left and right states to flux with Roe upwind method
void primToFlux(vector<vector<vector<double> > >& VL_x, vector<vector<vector<double> > >& VR_x, vector<vector<vector<double> > >& VL_y, vector<vector<vector<double> > >& VR_y, 
vector<vector<vector<double> > >& Fx, vector<vector<vector<double> > >& Fy, vector<vector<double> >& x_node, vector<vector<double> >& y_node, double gamma)
{
    double R, rhobar, ubar, vbar, ht_L, ht_R, htbar, abar, Uhat;
    double lambda1bar, lambda2bar, lambda3bar, lambda4bar, lambda1abs, lambda2abs, lambda3abs, lambda4abs;
    double deltarho, deltau, deltav, deltap, deltaw1, deltaw2, deltaw3, deltaw4;
    double nx, ny, A;
    vector<double> r1bar(4), r2bar(4), r3bar(4), r4bar(4), fL(4), fR(4);
    // Fx
    for (int i = 0; i < Fx[0].size(); i++){
        for (int j = 0; j < Fx[0][i].size(); j++){
            
            A = sqrt(pow(x_node[i][j+1]-x_node[i][j],2.0)+pow(y_node[i][j+1]-y_node[i][j],2.0));
            nx = (y_node[i][j+1]-y_node[i][j])/(A+1e-16);
            ny = -1.0*(x_node[i][j+1]-x_node[i][j])/(A+1e-16);

            R = sqrt(VR_x[0][i][j]/(VL_x[0][i][j]+1e-16));
            rhobar = R * VL_x[0][i][j];
            ubar = (R*VR_x[1][i][j]+VL_x[1][i][j])/(R+1.0);
            vbar = (R*VR_x[2][i][j]+VL_x[2][i][j])/(R+1.0);
            ht_L  = gamma/(gamma-1.0)*VL_x[3][i][j]/(VL_x[0][i][j]+1e-16) + 0.5*(pow(VL_x[1][i][j],2.0)+pow(VL_x[2][i][j],2.0));
            ht_R  = gamma/(gamma-1.0)*VR_x[3][i][j]/(VR_x[0][i][j]+1e-16) + 0.5*(pow(VR_x[1][i][j],2.0)+pow(VR_x[2][i][j],2.0));
            htbar = (R*ht_R+ht_L)/(R+1.0);
            abar = sqrt((gamma-1.0)*(htbar-0.5*(pow(ubar,2.0)+pow(vbar,2.0))));
            Uhat = ubar*nx + vbar*ny;

            lambda1bar = Uhat;
            lambda2bar = Uhat;
            lambda3bar = Uhat + abar;
            lambda4bar = Uhat - abar;

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
            if (fabs(lambda4bar) >= 2.0*0.1*abar) {
                lambda4abs = fabs(lambda4bar);
            } else {
                lambda4abs = pow(lambda4bar,2.0)/(4.0*0.1*abar+1e-16) + 0.1*abar;
            }

            deltarho = VR_x[0][i][j] - VL_x[0][i][j];
            deltau = VR_x[1][i][j] - VL_x[1][i][j];
            deltav = VR_x[2][i][j] - VL_x[2][i][j];
            deltap = VR_x[3][i][j] - VL_x[3][i][j];
            deltaw1 = deltarho - deltap/(pow(abar,2.0)+1e-16);
            deltaw2 = ny*deltau - nx*deltav;
            deltaw3 = nx*deltau + ny*deltav + deltap/(rhobar*abar+1e-16);
            deltaw4 = nx*deltau + ny*deltav - deltap/(rhobar*abar+1e-16);

            r1bar[0] = 1.0;
            r1bar[1] = ubar;
            r1bar[2] = vbar;
            r1bar[3] = 0.5*(pow(ubar,2.0)+pow(vbar,2.0));
            r2bar[0] = 0.0;
            r2bar[1] = ny*rhobar;
            r2bar[2] = -1.0*nx*rhobar;
            r2bar[3] = rhobar*(ny*ubar-nx*vbar);
            r3bar[0] = rhobar/(2.0*abar+1e-16);
            r3bar[1] = rhobar/(2.0*abar+1e-16)*(ubar+nx*abar);
            r3bar[2] = rhobar/(2.0*abar+1e-16)*(vbar+ny*abar);
            r3bar[3] = rhobar/(2.0*abar+1e-16)*(htbar+Uhat*abar);
            r4bar[0] = -rhobar/(2.0*abar+1e-16);
            r4bar[1] = -rhobar/(2.0*abar+1e-16)*(ubar-nx*abar);
            r4bar[2] = -rhobar/(2.0*abar+1e-16)*(vbar-ny*abar);
            r4bar[3] = -rhobar/(2.0*abar+1e-16)*(htbar-Uhat*abar);
            // conservative variables at left
            fL[0] = VL_x[0][i][j]*(VL_x[1][i][j]*nx+VL_x[2][i][j]*ny);
            fL[1] = VL_x[0][i][j]*VL_x[1][i][j]*(VL_x[1][i][j]*nx+VL_x[2][i][j]*ny)+VL_x[3][i][j]*nx;
            fL[2] = VL_x[0][i][j]*VL_x[2][i][j]*(VL_x[1][i][j]*nx+VL_x[2][i][j]*ny)+VL_x[3][i][j]*ny;
            fL[3] = VL_x[0][i][j]*ht_L*(VL_x[1][i][j]*nx+VL_x[2][i][j]*ny);
            // conservative variables at right
            fR[0] = VR_x[0][i][j]*(VR_x[1][i][j]*nx+VR_x[2][i][j]*ny);
            fR[1] = VR_x[0][i][j]*VR_x[1][i][j]*(VR_x[1][i][j]*nx+VR_x[2][i][j]*ny)+VR_x[3][i][j]*nx;
            fR[2] = VR_x[0][i][j]*VR_x[2][i][j]*(VR_x[1][i][j]*nx+VR_x[2][i][j]*ny)+VR_x[3][i][j]*ny;
            fR[3] = VR_x[0][i][j]*ht_R*(VR_x[1][i][j]*nx+VR_x[2][i][j]*ny);
            // flux at interface
            for (int d = 0; d < 4; d++){
                Fx[d][i][j] = 0.5*(fL[d]+fR[d]) - 0.5*(lambda1abs*deltaw1*r1bar[d]+lambda2abs*deltaw2*r2bar[d]+lambda3abs*deltaw3*r3bar[d]+lambda4abs*deltaw4*r4bar[d]);
            }
        }
    }
    // Fy
    for (int i = 0; i < Fy[0].size(); i++){
        for (int j = 0; j < Fy[0][i].size(); j++){
            
            A = sqrt(pow(x_node[i+1][j]-x_node[i][j],2.0)+pow(y_node[i+1][j]-y_node[i][j],2.0));
            nx = -1.0*(y_node[i+1][j]-y_node[i][j])/(A+1e-16);
            ny = (x_node[i+1][j]-x_node[i][j])/(A+1e-16);

            R = sqrt(VR_y[0][i][j]/(VL_y[0][i][j]+1e-16));
            rhobar = R * VL_y[0][i][j];
            ubar = (R*VR_y[1][i][j]+VL_y[1][i][j])/(R+1.0);
            vbar = (R*VR_y[2][i][j]+VL_y[2][i][j])/(R+1.0);
            ht_L  = gamma/(gamma-1.0)*VL_y[3][i][j]/(VL_y[0][i][j]+1e-16) + 0.5*(pow(VL_y[1][i][j],2.0)+pow(VL_y[2][i][j],2.0));
            ht_R  = gamma/(gamma-1.0)*VR_y[3][i][j]/(VR_y[0][i][j]+1e-16) + 0.5*(pow(VR_y[1][i][j],2.0)+pow(VR_y[2][i][j],2.0));
            htbar = (R*ht_R+ht_L)/(R+1.0);
            abar = sqrt((gamma-1.0)*(htbar-0.5*(pow(ubar,2.0)+pow(vbar,2.0))));
            Uhat = ubar*nx + vbar*ny;

            lambda1bar = Uhat;
            lambda2bar = Uhat;
            lambda3bar = Uhat + abar;
            lambda4bar = Uhat - abar;

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
            if (fabs(lambda4bar) >= 2.0*0.1*abar) {
                lambda4abs = fabs(lambda4bar);
            } else {
                lambda4abs = pow(lambda4bar,2.0)/(4.0*0.1*abar+1e-16) + 0.1*abar;
            }

            deltarho = VR_y[0][i][j] - VL_y[0][i][j];
            deltau = VR_y[1][i][j] - VL_y[1][i][j];
            deltav = VR_y[2][i][j] - VL_y[2][i][j];
            deltap = VR_y[3][i][j] - VL_y[3][i][j];
            deltaw1 = deltarho - deltap/(pow(abar,2.0)+1e-16);
            deltaw2 = ny*deltau - nx*deltav;
            deltaw3 = nx*deltau + ny*deltav + deltap/(rhobar*abar+1e-16);
            deltaw4 = nx*deltau + ny*deltav - deltap/(rhobar*abar+1e-16);

            r1bar[0] = 1.0;
            r1bar[1] = ubar;
            r1bar[2] = vbar;
            r1bar[3] = 0.5*(pow(ubar,2.0)+pow(vbar,2.0));
            r2bar[0] = 0.0;
            r2bar[1] = ny*rhobar;
            r2bar[2] = -1.0*nx*rhobar;
            r2bar[3] = rhobar*(ny*ubar-nx*vbar);
            r3bar[0] = rhobar/(2.0*abar+1e-16);
            r3bar[1] = rhobar/(2.0*abar+1e-16)*(ubar+nx*abar);
            r3bar[2] = rhobar/(2.0*abar+1e-16)*(vbar+ny*abar);
            r3bar[3] = rhobar/(2.0*abar+1e-16)*(htbar+Uhat*abar);
            r4bar[0] = -rhobar/(2.0*abar+1e-16);
            r4bar[1] = -rhobar/(2.0*abar+1e-16)*(ubar-nx*abar);
            r4bar[2] = -rhobar/(2.0*abar+1e-16)*(vbar-ny*abar);
            r4bar[3] = -rhobar/(2.0*abar+1e-16)*(htbar-Uhat*abar);
            // conservative variables at left
            fL[0] = VL_y[0][i][j]*(VL_y[1][i][j]*nx+VL_y[2][i][j]*ny);
            fL[1] = VL_y[0][i][j]*VL_y[1][i][j]*(VL_y[1][i][j]*nx+VL_y[2][i][j]*ny)+VL_y[3][i][j]*nx;
            fL[2] = VL_y[0][i][j]*VL_y[2][i][j]*(VL_y[1][i][j]*nx+VL_y[2][i][j]*ny)+VL_y[3][i][j]*ny;
            fL[3] = VL_y[0][i][j]*ht_L*(VL_y[1][i][j]*nx+VL_y[2][i][j]*ny);
            // conservative variables at right
            fR[0] = VR_y[0][i][j]*(VR_y[1][i][j]*nx+VR_y[2][i][j]*ny);
            fR[1] = VR_y[0][i][j]*VR_y[1][i][j]*(VR_y[1][i][j]*nx+VR_y[2][i][j]*ny)+VR_y[3][i][j]*nx;
            fR[2] = VR_y[0][i][j]*VR_y[2][i][j]*(VR_y[1][i][j]*nx+VR_y[2][i][j]*ny)+VR_y[3][i][j]*ny;
            fR[3] = VR_y[0][i][j]*ht_R*(VR_y[1][i][j]*nx+VR_y[2][i][j]*ny);
            // flux at interface
            for (int d = 0; d < 4; d++){
                Fy[d][i][j] = 0.5*(fL[d]+fR[d]) - 0.5*(lambda1abs*deltaw1*r1bar[d]+lambda2abs*deltaw2*r2bar[d]+lambda3abs*deltaw3*r3bar[d]+lambda4abs*deltaw4*r4bar[d]);
            }
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
void getResidual(vector<vector<vector<double> > >& Fx, vector<vector<vector<double> > >& Fy, vector<vector<vector<double> > >& R0, vector<vector<double> >& Area_x, vector<vector<double> >& Area_y)
{
    for (int d = 0; d < R0.size(); d++){
        for (int i = 0; i < R0[d].size(); i++){
            for (int j = 0; j < R0[d][i].size(); j++){
                R0[d][i][j] = -1.0*Fx[d][i][j]*Area_x[i][j] + Fx[d][i+1][j]*Area_x[i+1][j] - 1.0*Fy[d][i][j]*Area_y[i][j] + Fy[d][i][j+1]*Area_y[i][j+1];
            }
        }
    }
}

// Function: To get the L1 norm of residual
void getL1Norm(vector<vector<vector<double> > >& R0, vector<double>& current_resi_norm)
{
    vector<double> norm(4, 0.0);
    for (int d = 0; d < R0.size(); d++){
        for (int i = 0; i < R0[d].size(); i++){
            for (int j = 0; j < R0[d][i].size(); j++){
                norm[d] += fabs(R0[d][i][j]);
            }
        }
    }
    for (int i = 0; i < current_resi_norm.size(); i++){
        current_resi_norm[i] = norm[i]/(R0[0].size()*R0[0][0].size());
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

// Function: To get total pressure loss
void getPressLoss(vector<vector<vector<double> > >& V0, vector<vector<double> >& Area_x, double gamma, double & pressLoss) {
    int Nx, Ny;
    double p0_infi, H;
    
    Nx = V0[0].size();
    Ny = V0[0][0].size();
    H = 0.4;

    vector<double> V0_exit(4, 0.0);
    vector<double> p0_exit(Ny, 0.0);
    vector<double> Ma_exit(Ny, 0.0);

    p0_infi = 12270.0 * pow(1.0+0.5*(gamma-1.0)*16.0,gamma/(gamma-1.0));
    
    for (int j = 2; j < V0[0][0].size()-2; j++) {
        for (int d = 0; d < V0.size(); d++) {
            V0_exit[d] = 0.5*(3.0*V0[d][Nx-3][j]-V0[d][Nx-4][j]);
        }
        Ma_exit[j] = sqrt(pow(V0_exit[1],2.0)+pow(V0_exit[2],2.0))/(sqrt(gamma*V0_exit[3]/(V0_exit[0]+1e-16))+1e-16);
        p0_exit[j] = V0_exit[3]*pow(1.0+0.5*(gamma-1.0)*pow(Ma_exit[j],2.0), gamma/(gamma-1.0));
        pressLoss += (p0_infi - p0_exit[j]) * Area_x[Nx-4][j-2];
    }
    pressLoss = pressLoss/H;
}

