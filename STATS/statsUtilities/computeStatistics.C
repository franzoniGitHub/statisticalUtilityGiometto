#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <stdio.h>
#include <stdlib.h>
//include the local "settings.h" header file
#include "../SETTINGS.h"
using namespace std;

int main() {
try {
cout << "The programm \"computeStatistics.C\" is running\n";

//CLEAR OPTION
if(clear_bool) {
char answer='n';
cout << "You chose to clear the STATS case and the environment (also the output directories will be deleted): please, confirm (Y/n). Give any other answer to quit:\n";
cin >> answer;
if(answer=='Y') {
system("bash clearAll.sh");
}
else if(answer=='n'){;}
else {
cout << "Closing the program.\n";
}
}

//ENVIRONMENT OPTION
if(environment_bool) {
char answer='n';
cout << "You chose to create the environment: please, confirm (Y/n). Give any other answer to quit:\n";
cin >> answer;
if(answer=='Y') {
system("bash environment.sh"); 
}
else if(answer=='n'){;}
else {
cout << "Closing the program.\n";
} 
}

//Get number of time steps
int nsteps=0;
ifstream nstepsFile("../stats/nSteps.dat");
if(!nstepsFile) {
cout << "Unable to open nSteps.dat in STATS/stats directory. Closing the program.\n";
return 1;
}
nstepsFile >> nsteps;
cout << "The number of time steps is " << nsteps << endl;
nstepsFile.close();

//Build the matrix for the results
vector<vector<double>> matrixOfResults; //first index=rows, second index=columns
for(int i=0; i<Nz; i++) matrixOfResults.push_back({0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}); //18 columns: this is just for initialization.
/*
LEGEND:
Column 0=z values
Column 1=mean b
Column 2=mean Ux
Column 3=mean Uy
Column 4=mean Uz
Column 5=rms b
Column 6=rms Ux
Column 7=rms Uy
Column 8=rms Uz
Column 9=covar UxUz
Column 10=covar Uxb
Column 11=covar Uyb
Column 12=covar Uzb
Column 13=TKE
Column 14=delta/eta -> nut/nu
Column 15=mean dUx/dz
Column 16=mean DIS
Column 17=mean db/dz
*/

//Reading z values: you need to put z.dat file into the STATS directory
//Column 0 plus 14
cout << "Reading the z values from z.dat\n";

ifstream Zfile("../z.dat");
if(Zfile) {
double previousZ=0.0, currentZ=0.0;
Zfile >> previousZ;
for(int i=0; i<Nz; i++) {
Zfile >> currentZ;
matrixOfResults[i][0]=(currentZ+previousZ)/(2.0*Z_norm);
previousZ=currentZ;
}
Zfile.close();
} else {
cout << "Unable to open z.dat in STATS directory. Closing the program.\n";
return 1;
}

//###########################################################################################################################################################
//First order statistics

if(firstOrder_bool) {//1

cout << "Computing the first order statistics\n";

int zcounter=1;
double b=0.0;
vector<double> U={0.0, 0.0, 0.0};

//Columns 1 to 4
for(int i=0; i<nsteps; i++){//2
if(i%20==0) cout << "   current step=" << i << endl;

string filenameb="../stats/b"+to_string(i);
string filenameU="../stats/U"+to_string(i);
string filenameNUT="../stats/nut"+to_string(i);

ifstream bfile(filenameb);
ifstream Ufile(filenameU);
ifstream NUTfile(filenameNUT);

if(bfile && Ufile && NUTfile) {//3

zcounter=1;
for(int j=0; j<22; j++) bfile.ignore(500, '\n');
for(int j=0; j<22; j++) Ufile.ignore(500, '\n');
for(int j=0; j<22; j++) NUTfile.ignore(500, '\n');

for(int j=0; j<Nx*Ny*Nz; j++){//4

bfile >> b;
matrixOfResults[zcounter-1][1]+=b/(b_norm*double(Nx*Ny*nsteps));

NUTfile >> b;
matrixOfResults[zcounter-1][14]+=b/(nu*double(Nx*Ny*nsteps));

Ufile.ignore(1);
Ufile >> U[0] >> U[1] >> U[2];
Ufile.ignore(500, '\n');
matrixOfResults[zcounter-1][2]+=U[0]/(U_norm*double(Nx*Ny*nsteps));
matrixOfResults[zcounter-1][3]+=U[1]/(U_norm*double(Nx*Ny*nsteps));
matrixOfResults[zcounter-1][4]+=U[2]/(U_norm*double(Nx*Ny*nsteps));
if(j==zcounter*Nx*Ny-1) zcounter++;
}//4

bfile.close();
Ufile.close();
NUTfile.close();
} else {
cout << "Unable to open some of the files. Closing the program.\n";
return 1;
}//3
}//2
//******* Write all the output files needed ***********//
ofstream bmean("../output/firstOrder/bmean.dat");
ofstream Umean("../output/firstOrder/Umean.dat");
ofstream nut_nu("../output/firstOrder/nut_nu.dat");

cout << "Writing all the output files required.\n";
/*
LEGEND:
Column 0=z values
Column 1=mean b
Column 2=mean Ux
Column 3=mean Uy
Column 4=mean Uz
Column 5=rms b
Column 6=rms Ux
Column 7=rms Uy
Column 8=rms Uz
Column 9=covar UxUz
Column 10=covar Uxb
Column 11=covar Uyb
Column 12=covar Uzb
Column 13=TKE
Column 14=nut/nu
Column 15=mean dUx/dz
Column 16=mean DIS
Column 17=mean db/dz
*/

for(int i=0; i<Nz; i++) {
bmean      << matrixOfResults[i][0] << "  " << matrixOfResults[i][1]  << endl;
Umean      << matrixOfResults[i][0] << "  " << matrixOfResults[i][2]  << "  " << matrixOfResults[i][3]  << "  " << matrixOfResults[i][4]  << endl;
nut_nu     << matrixOfResults[i][0] << "  " << matrixOfResults[i][14]  << endl;
}

//Building the Umean file for the OpenFOAM-based utility that calculates the dissipation
ofstream OFfile("../Umean");
if(OFfile){
OFfile << "FoamFile{version     2.0;    format      ascii;    class       volVectorField;    object      Umean;}\n";
OFfile << "dimensions      [0 1 -1 0 0 0 0]; internalField   nonuniform List<vector>\n" << Nx*Ny*Nz << endl;
OFfile << "(\n";
zcounter=1;
for(int j=0; j<Nx*Ny*Nz; j++) {
OFfile << "(" << matrixOfResults[zcounter-1][2]*U_norm << "  " << matrixOfResults[zcounter-1][3]*U_norm << "  " << matrixOfResults[zcounter-1][4]*U_norm << ")\n";
if(j==zcounter*Nx*Ny-1) zcounter++;
}
OFfile << ");\n";
OFfile << "boundaryField{    front    {        type            cyclic;    }\n";
OFfile << "    back    {        type            cyclic;    }\n";
OFfile << "    left    {        type            cyclic;    }\n";
OFfile << "    right    {        type            cyclic;    }\n";
OFfile << "    floor    {        type            zeroGradient;    }\n";
OFfile << "    ceiling    {        type            zeroGradient;    }\n";
OFfile << "}";
OFfile.close();
}

}//1

//###########################################################################################################################################################
//Second order statistics

if(secondOrder_bool){//1

cout << "Computing the second order statistics\n";

system("bash forSecondOrder.sh");

//Check that the general number of time steps is equal to that from forSecondOrder.sh
int nsteps2=0;
nstepsFile.open("../stats/nSteps2.dat");
if(!nstepsFile) {
cout << "Unable to open nSteps2.dat in STATS/stats directory. Closing the program.\n";
return 1;
}
nstepsFile >> nsteps2;
if(nsteps2 != nsteps) {
cout << "The general number of time steps is different from that from forSecondOrder.sh. Closing the program.\n";
return 1;
} 
nstepsFile.close();

if(!firstOrder_bool){
ifstream bMeanfile("../output/firstOrder/bmean.dat");
ifstream UMeanfile("../output/firstOrder/Umean.dat");
if(bMeanfile && UMeanfile) {
double useless=0., Ux=0., Uy=0., Uz=0., b=0.;
for(int i=0; i<Nz; i++) {
bMeanfile >> useless >> b;
UMeanfile >> useless >> Ux >> Uy >> Uz;
matrixOfResults[i][1]=b;
matrixOfResults[i][2]=Ux;
matrixOfResults[i][3]=Uy;
matrixOfResults[i][4]=Uz;
}
bMeanfile.close();
UMeanfile.close();
} else {
cout << "Unable to open bmean.dat or Umean.dat in the output directory. Closing the program.\n";
return 1;
}
}

int zcounter=1;
double b=0.0;
double DIS=0.0;
vector<double> U={0.0, 0.0, 0.0};

for(int i=0; i<nsteps; i++){//2
if(i%20==0) cout << "   current step=" << i << endl;

string filenameb="../stats/b"+to_string(i);
string filenameU="../stats/U"+to_string(i);
string filenameDIS="../stats/DIS"+to_string(i);
ifstream bfile(filenameb);
ifstream Ufile(filenameU);
ifstream DISfile(filenameDIS);

if( bfile && Ufile && DISfile) {//3
zcounter=1;
for(int j=0; j<22; j++) bfile.ignore(500, '\n');
for(int j=0; j<22; j++) Ufile.ignore(500, '\n');
for(int j=0; j<22; j++) DISfile.ignore(500, '\n');
for(int j=0; j<Nx*Ny*Nz; j++){//4
bfile >> b;
matrixOfResults[zcounter-1][5]+=(b/b_norm-matrixOfResults[zcounter-1][1])*(b/b_norm-matrixOfResults[zcounter-1][1]);
DISfile >> DIS;
matrixOfResults[zcounter-1][16]+=(DIS*Z_norm)/pow(U_norm, 3.0);
Ufile.ignore(1);
Ufile >> U[0] >> U[1] >> U[2];
U[0]/=U_norm;
U[1]/=U_norm;
U[2]/=U_norm;
Ufile.ignore(500, '\n');
matrixOfResults[zcounter-1][6]+=(U[0]-matrixOfResults[zcounter-1][2])*(U[0]-matrixOfResults[zcounter-1][2]);
matrixOfResults[zcounter-1][7]+=(U[1]-matrixOfResults[zcounter-1][3])*(U[1]-matrixOfResults[zcounter-1][3]);
matrixOfResults[zcounter-1][8]+=(U[2]-matrixOfResults[zcounter-1][4])*(U[2]-matrixOfResults[zcounter-1][4]);

matrixOfResults[zcounter-1][9]+=(U[0]-matrixOfResults[zcounter-1][2])*(U[2]-matrixOfResults[zcounter-1][4]);
matrixOfResults[zcounter-1][10]+=(U[0]-matrixOfResults[zcounter-1][2])*(b/b_norm-matrixOfResults[zcounter-1][1]);
matrixOfResults[zcounter-1][11]+=(U[1]-matrixOfResults[zcounter-1][3])*(b/b_norm-matrixOfResults[zcounter-1][1]);
matrixOfResults[zcounter-1][12]+=(U[2]-matrixOfResults[zcounter-1][4])*(b/b_norm-matrixOfResults[zcounter-1][1]);
if(j==zcounter*Nx*Ny-1) zcounter++;
}//4
bfile.close();
Ufile.close();
DISfile.close();
} else {
cout << "Unable to open some of the files. Closing the program.\n";
return 1;
}//3
}//2

for(int i=0; i<Nz; i++) {//2'
matrixOfResults[i][5]=sqrt(matrixOfResults[i][5]/double(Nx*Ny*nsteps));
matrixOfResults[i][16]/=double(Nx*Ny*nsteps);
matrixOfResults[i][6]=sqrt(matrixOfResults[i][6]/double(Nx*Ny*nsteps));
matrixOfResults[i][7]=sqrt(matrixOfResults[i][7]/double(Nx*Ny*nsteps));
matrixOfResults[i][8]=sqrt(matrixOfResults[i][8]/double(Nx*Ny*nsteps));
matrixOfResults[i][9]/=double(Nx*Ny*nsteps);
matrixOfResults[i][10]/=double(Nx*Ny*nsteps);
matrixOfResults[i][11]/=double(Nx*Ny*nsteps);
matrixOfResults[i][12]/=double(Nx*Ny*nsteps);
matrixOfResults[i][13]=0.5*(pow(matrixOfResults[i][6], 2)+pow(matrixOfResults[i][7], 2)+pow(matrixOfResults[i][8], 2));
}//2'

//Compute the stationary identity
//We use Newton approx on boundaries, central scheme otherwise

cout << "Now computing the stationary identity and writing its output.\n";
{//2''
vector<double> tauXZ;
tauXZ.resize(Nz, 0.0);
tauXZ[0]=(matrixOfResults[1][2]-matrixOfResults[0][2])/((matrixOfResults[1][0]-matrixOfResults[0][0])*sqrt(Gr))-matrixOfResults[0][9];
tauXZ[Nz-1]=(matrixOfResults[Nz-1][2]-matrixOfResults[Nz-2][2])/((matrixOfResults[Nz-1][0]-matrixOfResults[Nz-2][0])*sqrt(Gr))-matrixOfResults[Nz-1][9];

vector<double> taubZ;
taubZ.resize(Nz, 0.0);
taubZ[0]=(matrixOfResults[1][1]-matrixOfResults[0][1])/((matrixOfResults[1][0]-matrixOfResults[0][0])*sqrt(Gr)*Pr)-matrixOfResults[0][12];
taubZ[Nz-1]=(matrixOfResults[Nz-1][1]-matrixOfResults[Nz-2][1])/((matrixOfResults[Nz-1][0]-matrixOfResults[Nz-2][0])*sqrt(Gr)*Pr)-matrixOfResults[Nz-1][12];
for(int i=1; i<Nz-1; i++) {
tauXZ[i]=(matrixOfResults[i+1][2]-matrixOfResults[i-1][2])/((matrixOfResults[i+1][0]-matrixOfResults[i-1][0])*sqrt(Gr))-matrixOfResults[i][9];
taubZ[i]=(matrixOfResults[i+1][1]-matrixOfResults[i-1][1])/((matrixOfResults[i+1][0]-matrixOfResults[i-1][0])*sqrt(Gr)*Pr)-matrixOfResults[i][12];
}

//direct output
ofstream momentumBalance("../output/secondOrder/momentumStationary.dat");
ofstream energyBalance("../output/secondOrder/energyStationary.dat");

momentumBalance<<matrixOfResults[0][1]*sin(alpha)<<"   "<<-(tauXZ[1]-tauXZ[0])/(matrixOfResults[1][0]-matrixOfResults[0][0])<<"   "<<tauXZ[1]-tauXZ[0]<<"   "<<matrixOfResults[1][0]-matrixOfResults[0][0] << endl;
energyBalance<<matrixOfResults[0][2]*sin(alpha)<<"   "<<(taubZ[1]-taubZ[0])/(matrixOfResults[1][0]-matrixOfResults[0][0])<< "   " << taubZ[1]-taubZ[0]<<"   "<<matrixOfResults[1][0]-matrixOfResults[0][0] << endl;
for(int i=1; i<Nz-1; i++) {
momentumBalance << matrixOfResults[i][1]*sin(alpha) << "   " << -(tauXZ[i+1]-tauXZ[i-1])/(matrixOfResults[i+1][0]-matrixOfResults[i-1][0])<<"   "<<tauXZ[i+1]-tauXZ[i-1]<<"   "<<matrixOfResults[i+1][0]-matrixOfResults[i-1][0] << endl;
energyBalance << matrixOfResults[i][2]*sin(alpha) << "   " << (taubZ[i+1]-taubZ[i-1])/(matrixOfResults[i+1][0]-matrixOfResults[i-1][0]) <<"   "<< taubZ[i+1]-taubZ[i-1]<<"   "<<matrixOfResults[i+1][0]-matrixOfResults[i-1][0] << endl;
}
momentumBalance << matrixOfResults[Nz-1][1]*sin(alpha) << "   " << -(tauXZ[Nz-1]-tauXZ[Nz-2])/(matrixOfResults[Nz-1][0]-matrixOfResults[Nz-2][0])<<"  "<<tauXZ[Nz-1]-tauXZ[Nz-2]<<"  " <<matrixOfResults[Nz-1][0]-matrixOfResults[Nz-2][0] << endl;
energyBalance << matrixOfResults[Nz-1][2]*sin(alpha) << "   " << (taubZ[Nz-1]-taubZ[Nz-2])/(matrixOfResults[Nz-1][0]-matrixOfResults[Nz-2][0])<<"  "<<taubZ[Nz-1]-taubZ[Nz-2]<<"  " <<matrixOfResults[Nz-1][0]-matrixOfResults[Nz-2][0] << endl;

momentumBalance.close();
energyBalance.close();
}//2''


//******* Write all the output files needed ***********//
ofstream ALLmatrix("../output/secondOrder/ALLmatrix.dat");
ofstream brms("../output/secondOrder/brms.dat");
ofstream Urms("../output/secondOrder/Urms.dat");
ofstream UxUzcovar("../output/secondOrder/UxUzcovar.dat");
ofstream Ubcovar("../output/secondOrder/Ubcovar.dat");
ofstream TKE("../output/secondOrder/TKE.dat");
ofstream TKEbalance("../output/secondOrder/TKEbalance.dat");

cout << "Writing all the output files required.\n";

for(int i=0; i<Nz; i++) {
if(i==0) ALLmatrix << "z\tmeanb\tmeanUx\tmeanUy\tmeanUz\trmsb\trmsUx\trmsUy\trmsUz\tcovarUxUz\tcovarUxb\tcovarUyb\tcovarUzb\tTKE\tnut/nu\tdUx/dz\tDIS\tdb/dz\n";
for(int j=0; j<18; j++) ALLmatrix << matrixOfResults[i][j] << "\t";
ALLmatrix << endl;
}
ALLmatrix.close();

for(int i=0; i<Nz; i++) {
brms       << matrixOfResults[i][0] << "  " << matrixOfResults[i][5]  << endl;
Urms       << matrixOfResults[i][0] << "  " << matrixOfResults[i][6]  << "  " << matrixOfResults[i][7]  << "  " << matrixOfResults[i][8]  << endl;
UxUzcovar  << matrixOfResults[i][0] << "  " << matrixOfResults[i][9]  << endl;
Ubcovar    << matrixOfResults[i][0] << "  " << matrixOfResults[i][10] << "  " << matrixOfResults[i][11] << "  " << matrixOfResults[i][12] << endl;
TKE        << matrixOfResults[i][0] << "  " << matrixOfResults[i][13]  << endl;
TKEbalance << matrixOfResults[i][0] << "  " ;

if(i==0){
TKEbalance << -matrixOfResults[i][9]*(matrixOfResults[1][2]-matrixOfResults[0][2])/(matrixOfResults[1][0]-matrixOfResults[0][0]);                     ///    d<U>/dz
} else if(i==Nz-1) {
TKEbalance << -matrixOfResults[i][9]*(matrixOfResults[Nz-1][2]-matrixOfResults[Nz-2][2])/(matrixOfResults[Nz-1][0]-matrixOfResults[Nz-2][0]);
} else {
TKEbalance << -matrixOfResults[i][9]*(matrixOfResults[i+1][2]-matrixOfResults[i-1][2])/(matrixOfResults[i+1][0]-matrixOfResults[i-1][0]);
}
TKEbalance << "  " << sin(alpha)*matrixOfResults[i][10] << "  " << cos(alpha)*matrixOfResults[i][12] << "  " << matrixOfResults[i][16] << endl;

}

}//1

//###########################################################################################################################################################
//Compute autocorrelation

if(autocorrelation_bool) {

vector<double> U={0.0, 0.0, 0.0};
vector<vector<double>> innerMostVectorY;//just for initialization.
vector<vector<vector<double>>> innerVectorY;
vector<vector<vector<double>>> varUY;
for(int i=0; i<Nz; i++) {
innerMostVectorY.push_back({0.0, 0.0, 0.0});
}
for(int i=0; i<Ny; i++) {
innerVectorY.push_back(innerMostVectorY);
varUY.push_back(innerMostVectorY);
}

vector<vector<vector<double>>> autocorrelationY; //first index=y, second index=z, third index=U component
for(int i=0; i<Ny; i++) autocorrelationY.push_back(innerMostVectorY);

vector<vector<double>> innerMostVectorX;//just for initialization.
vector<vector<vector<double>>> innerVectorX;
vector<vector<vector<double>>> varUX;
for(int i=0; i<Nz; i++) {
innerMostVectorX.push_back({0.0, 0.0, 0.0});
}
for(int i=0; i<Nx; i++) {
innerVectorX.push_back(innerMostVectorX);
varUX.push_back(innerMostVectorX);
}

vector<vector<vector<double>>> autocorrelationX; //first index=x, second index=z, third index=U component
for(int i=0; i<Nx; i++) autocorrelationX.push_back(innerMostVectorX);


if(!firstOrder_bool){
ifstream meanUfile("../output/firstOrder/Umean.dat");
if (! meanUfile) {
cout <<"No Umean.dat found in the output directory.\n";
return 1;
}
double useless=0.0;
for(int i=0; i<Nz; i++) {
meanUfile >> useless >> U[0] >> U[1] >> U[2];
matrixOfResults[i][2]=U[0];
matrixOfResults[i][3]=U[1];
matrixOfResults[i][4]=U[2];
}

meanUfile.close();
}

cout << "Beginning of the time loop to calculate the autocorrelation functions.\n";

for(int index=0; index<nsteps; index++){

if(index%20==0) cout << "   current step=" << index << endl;

string filenameU="../stats/U"+to_string(index);
ifstream Ufile(filenameU);

vector<vector<vector<vector<double>>>> matrixOfInput; //first index=x, second index=y, third index=z, fourth index=U component
for(int i=0; i<Nx; i++) matrixOfInput.push_back(innerVectorY);

for(int j=0; j<22; j++) Ufile.ignore(500, '\n');

for(int k=0; k<Nz; k++){
for(int j=0; j<Ny; j++){
for(int i=0; i<Nx; i++){
U={0.0, 0.0, 0.0};
Ufile.ignore(1);
Ufile >> U[0] >> U[1] >> U[2];
Ufile.ignore(500, '\n');
matrixOfInput[i][j][k][0]=U[0]/U_norm;
matrixOfInput[i][j][k][1]=U[1]/U_norm;
matrixOfInput[i][j][k][2]=U[2]/U_norm;
}
}
}

for(int k=0; k<Nz; k++){
for(int j=0; j<Ny; j++){
for(int i=0; i<Nx; i++){

varUY[j][k][0]+=(matrixOfInput[i][j][k][0]-matrixOfResults[k][2])*(matrixOfInput[i][j][k][0]-matrixOfResults[k][2])/(double(Nx*nsteps));
varUY[j][k][1]+=(matrixOfInput[i][j][k][1]-matrixOfResults[k][3])*(matrixOfInput[i][j][k][1]-matrixOfResults[k][3])/(double(Nx*nsteps));
varUY[j][k][2]+=(matrixOfInput[i][j][k][2]-matrixOfResults[k][4])*(matrixOfInput[i][j][k][2]-matrixOfResults[k][4])/(double(Nx*nsteps));

autocorrelationY[j][k][0]+=(matrixOfInput[i][0][k][0]-matrixOfResults[k][2])*(matrixOfInput[i][j][k][0]-matrixOfResults[k][2])/(double(Nx*nsteps));
autocorrelationY[j][k][1]+=(matrixOfInput[i][0][k][1]-matrixOfResults[k][3])*(matrixOfInput[i][j][k][1]-matrixOfResults[k][3])/(double(Nx*nsteps));
autocorrelationY[j][k][2]+=(matrixOfInput[i][0][k][2]-matrixOfResults[k][4])*(matrixOfInput[i][j][k][2]-matrixOfResults[k][4])/(double(Nx*nsteps));
}
}
}

for(int k=0; k<Nz; k++){
for(int j=0; j<Nx; j++){
for(int i=0; i<Ny; i++){
varUX[j][k][0]+=(matrixOfInput[j][i][k][0]-matrixOfResults[k][2])*(matrixOfInput[j][i][k][0]-matrixOfResults[k][2])/(double(Ny*nsteps));
varUX[j][k][1]+=(matrixOfInput[j][i][k][1]-matrixOfResults[k][3])*(matrixOfInput[j][i][k][1]-matrixOfResults[k][3])/(double(Ny*nsteps));
varUX[j][k][2]+=(matrixOfInput[j][i][k][2]-matrixOfResults[k][4])*(matrixOfInput[j][i][k][2]-matrixOfResults[k][4])/(double(Ny*nsteps));

autocorrelationX[j][k][0]+=(matrixOfInput[0][i][k][0]-matrixOfResults[k][2])*(matrixOfInput[j][i][k][0]-matrixOfResults[k][2])/(double(Ny*nsteps));
autocorrelationX[j][k][1]+=(matrixOfInput[0][i][k][1]-matrixOfResults[k][3])*(matrixOfInput[j][i][k][1]-matrixOfResults[k][3])/(double(Ny*nsteps));
autocorrelationX[j][k][2]+=(matrixOfInput[0][i][k][2]-matrixOfResults[k][4])*(matrixOfInput[j][i][k][2]-matrixOfResults[k][4])/(double(Ny*nsteps));
}
}
}


Ufile.close();
}

cout << "Computing the autocorrelation function rho(y, z).\n";

for(int k=0; k<Nz; k++) {
for(int j=0; j<Ny; j++) {
autocorrelationY[j][k][0]/=sqrt(varUY[0][k][0]*varUY[j][k][0]);
autocorrelationY[j][k][1]/=sqrt(varUY[0][k][1]*varUY[j][k][1]);
autocorrelationY[j][k][2]/=sqrt(varUY[0][k][2]*varUY[j][k][2]);
}
}

cout << "Computing the autocorrelation function rho(x, z).\n";

for(int k=0; k<Nz; k++) {
for(int j=0; j<Nx; j++) {
autocorrelationX[j][k][0]/=sqrt(varUX[0][k][0]*varUX[j][k][0]);
autocorrelationX[j][k][1]/=sqrt(varUX[0][k][1]*varUX[j][k][1]);
autocorrelationX[j][k][2]/=sqrt(varUX[0][k][2]*varUX[j][k][2]);
}
}

cout << "Writing all the output required.\n";

ofstream selection1Y("../output/autocorrelation/alongY/selection1Y.dat");
ofstream selection2Y("../output/autocorrelation/alongY/selection2Y.dat");
ofstream selection3Y("../output/autocorrelation/alongY/selection3Y.dat");
ofstream selection4Y("../output/autocorrelation/alongY/selection4Y.dat");
ofstream selection5Y("../output/autocorrelation/alongY/selection5Y.dat");


for(int j=0; j<Ny; j++) {
selection1Y << (j*Ly+Ly/2.0)/(Ny*Ly) << "  " << autocorrelationY[j][k1][0] << "  " << autocorrelationY[j][k1][1] << "  " << autocorrelationY[j][k1][2] << endl;
selection2Y << (j*Ly+Ly/2.0)/(Ny*Ly) << "  " << autocorrelationY[j][k2][0] << "  " << autocorrelationY[j][k2][1] << "  " << autocorrelationY[j][k2][2] << endl;
selection3Y << (j*Ly+Ly/2.0)/(Ny*Ly) << "  " << autocorrelationY[j][k3][0] << "  " << autocorrelationY[j][k3][1] << "  " << autocorrelationY[j][k3][2] << endl;
selection4Y << (j*Ly+Ly/2.0)/(Ny*Ly) << "  " << autocorrelationY[j][k4][0] << "  " << autocorrelationY[j][k4][1] << "  " << autocorrelationY[j][k4][2] << endl;
selection5Y << (j*Ly+Ly/2.0)/(Ny*Ly) << "  " << autocorrelationY[j][k5][0] << "  " << autocorrelationY[j][k5][1] << "  " << autocorrelationY[j][k5][2] << endl;
}

selection1Y.close();
selection2Y.close();
selection3Y.close();
selection4Y.close();
selection5Y.close();

ofstream selection1X("../output/autocorrelation/alongX/selection1X.dat");
ofstream selection2X("../output/autocorrelation/alongX/selection2X.dat");
ofstream selection3X("../output/autocorrelation/alongX/selection3X.dat");
ofstream selection4X("../output/autocorrelation/alongX/selection4X.dat");
ofstream selection5X("../output/autocorrelation/alongX/selection5X.dat");

for(int j=0; j<Nx; j++) {
selection1X << (j*Lx+Lx/2.0)/(Nx*Lx) << "  " << autocorrelationX[j][k1][0] << "  " << autocorrelationX[j][k1][1] << "  " << autocorrelationX[j][k1][2] << endl;
selection2X << (j*Lx+Lx/2.0)/(Nx*Lx) << "  " << autocorrelationX[j][k2][0] << "  " << autocorrelationX[j][k2][1] << "  " << autocorrelationX[j][k2][2] << endl;
selection3X << (j*Lx+Lx/2.0)/(Nx*Lx) << "  " << autocorrelationX[j][k3][0] << "  " << autocorrelationX[j][k3][1] << "  " << autocorrelationX[j][k3][2] << endl;
selection4X << (j*Lx+Lx/2.0)/(Nx*Lx) << "  " << autocorrelationX[j][k4][0] << "  " << autocorrelationX[j][k4][1] << "  " << autocorrelationX[j][k4][2] << endl;
selection5X << (j*Lx+Lx/2.0)/(Nx*Lx) << "  " << autocorrelationX[j][k5][0] << "  " << autocorrelationX[j][k5][1] << "  " << autocorrelationX[j][k5][2] << endl;
}

}
//###########################################################################################################################################################

cout << "END of this program. Bye!\n";
return 0;
}
catch(...) {
cout << "An exception has been caught.\n";
return 1;
}
}//MAIN
