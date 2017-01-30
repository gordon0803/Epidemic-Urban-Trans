#include <ostream>
#include "vnode.h"

using namespace std;
using namespace vnodelp;

template<typename var_type>
void SIER(int n,var_type*yp,const var_typr*y,void*param) //
{
//define parameters
  double sigma; 
  double gamma;
  double betainn; // within zone contact rate
  double betatra[3]; // contactrate in travel mode m, m=1,2,3
  double c12[3];  // ratio of people move from i to j who choose travel mode m, i,j=1,2,3
  double c13[3];
  double c23[3];
  double c21[3];
  double c31[3];
  double c32[3];
  int Np[3]; // population of patch i
  int N[3][3]; // the amount of people currently at patch j who are the residents of patch i
  double g[3]; // total departure rate of patch i
  double m[3][3]; // the rate of movement from patch i to patch j
  double r[3][3]; // return rate of people who have traveled from patch i to patch j

  //ODE

  yp[0]=g[0]*y[0]-r[0][1]*y[1]-r[0][2]*y[2]+beteinn*(y[0]*(y[3]+y[15]+y[27]))/(Np[0]); //S_11 dot
  
  yp[1]=r[0][1]*y[1]-g[0]*m[0][1]*y[0]+betainn*y[1]*(y[4]+y[16]+y[28])/(Np[1])
	  +c12[0]*betatra[0]*y[1]*(c12[0]*y[4]+c13[0]*y[5])/(c12[0]*N[0][1]+c13[0]*N[0][2])
	  +c12[1]*betatra[1]*y[1]*(c12[1]*y[4]+c13[1]*y[5])/(c12[1]*N[0][1]+c13[1]*N[0][2])
	  +c12[2]*betatra[2]*y[1]*(c12[2]*y[4]+c13[2]*y[5])/(c12[2]*N[0][1]+c13[2]*N[0][2]); //S_12 dot
  
  yp[2]=r[0][2]*y[2]-g[0]*m[0][2]*y[0]+betainn*y[2]*(y[5]+y[17]+y[29])/(Np[2])
	  +c13[0]*betatra[0]*y[2]*(c12[0]*y[4]+c13[0]*y[5]+c23[0]*y[17])/(c12[0]*N[0][1]+c13[0]*N[0][2]z+c23[0]*N[1][2])
      +c13[1]*betatra[1]*y[2]*(c12[1]*y[4]+c13[1]*y[5]+c23[1]*y[17])/(c12[1]*N[0][1]+c13[1]*N[0][2]z+c23[1]*N[1][2])
	  +c13[2]*betatra[2]*y[2]*(c12[2]*y[4]+c13[2]*y[5]+c23[2]*y[17])/(c12[2]*N[0][1]+c13[2]*N[0][2]z+c23[2]*N[1][2]); //S_13 dot
  
  yp[3]=g[0]*y[3]-r[0][1]*y[4]-r[0][2]*y[5]-sigma*y[6]+gamma*y[3]; //I_11 dot
  
  yp[4]=r[0][1]*y[4]-g[0]*m[0][1]*y[3]-sigma*y[7]+gamma*y[4]; //I_12 dot
  
  yp[5]=r[0][2]*y[5]-g[0]*m[0][2]*y[4]-sigma*y[8]+gamma*y[5]; //I_13 dot
  
  yp[6]=g[0]*y[6]-r[0][1]*y[7]-r[0][2]*y[8]-beteinn*(y[0]*(y[3]+y[15]+y[27]))/(Np[0])+sigma*y[6]; //E_11 dot
  
  yp[7]=r[0][1]*y[7]-g[0]m[0][1]*y[6]-g[0]*m[0][1]-betainn*y[1]*(y[4]+y[16]+y[28])/(Np[1])
	  -c12[0]*betatra[0]*y[1]*(c12[0]*y[4]+c13[0]*y[5])/(c12[0]*N[0][1]+c13[0]*N[0][2])
	  -c12[1]*betatra[1]*y[1]*(c12[1]*y[4]+c13[1]*y[5])/(c12[1]*N[0][1]+c13[1]*N[0][2])
	  -c12[2]*betatra[2]*y[1]*(c12[2]*y[4]+c13[2]*y[5])/(c12[2]*N[0][1]+c13[2]*N[0][2])+sigma*y[7];//E_12 dot
  
  yp[8]=r[0][2]*y[8]-g[0]m[0][2]*y[6]-betainn*y[2]*(y[5]+y[17]+y[29])/(Np[2])
	  -c13[0]*betatra[0]*y[2]*(c12[0]*y[4]+c13[0]*y[5]+c23[0]*y[17])/(c12[0]*N[0][1]+c13[0]*N[0][2]z+c23[0]*N[1][2])
      -c13[1]*betatra[1]*y[2]*(c12[1]*y[4]+c13[1]*y[5]+c23[1]*y[17])/(c12[1]*N[0][1]+c13[1]*N[0][2]z+c23[1]*N[1][2])
	  -c13[2]*betatra[2]*y[2]*(c12[2]*y[4]+c13[2]*y[5]+c23[2]*y[17])/(c12[2]*N[0][1]+c13[2]*N[0][2]z+c23[2]*N[1][2]); //E_13 dot

  yp[9]=g[0]*y[9]-r[0][1]*y[10]-r[0][2]*y[11]-gamma*y[3]; //R_11 dot
 
  yp[10]=r[0][1]*y[10]-g[0]*m[0][1]*y[9]-gamma*y[4]; //R_12 dot

  yp[11]=r[0][2]*y[11]-g[0]*m[0][1]*y[9]-gamma*y[5]; //R_13 dot
  
  yp[12]=r[1][0]*y[12]-g[1]*m[1][0]*y[13]+betainn*y[12]*(y[3]+y[15]+y[27])/(Np[0])
	  +c21[0]*betatra[0]*y[12]*(c21[0]*y[15]+c31[0]*y[27])/(c21[0]*N[1][0]+c31[0]*N[2][0])
	  +c21[1]*betatra[1]*y[12]*(c21[1]*y[15]+c31[1]*y[27])/(c21[1]*N[1][0]+c31[1]*N[2][0])
	  +c21[2]*betatra[2]*y[12]*(c21[2]*y[15]+c31[2]*y[27])/(c21[2]*N[1][0]+c31[2]*N[2][0]); //S_21 dot

  yp[13]=g[1]*y[13]-r[1][0]*y[12]-r[1][2]*y[14]+betainn*(y[13]*(y[4]+y[16]+y[28]))/(Np[1]); //S_22 dot
  	  
  yp[14]=r[1][2]*y[14]-g[1]*m[1][2]*y[13]+betainn*y[14]*(y[5]+y[17]+y[29])/(Np[2])
	  +c23[0]*betatra[0]*y[12]*(c13[0]*y[5]+c23[0]*y[17])/(c13[0]*N[0][2]+c23[0]*N[1][2])
	  +c23[1]*betatra[1]*y[12]*(c13[1]*y[5]+c23[1]*y[17])/(c13[1]*N[0][2]+c23[1]*N[1][2])
	  +c23[2]*betatra[2]*y[12]*(c13[2]*y[5]+c23[2]*y[17])/(c13[2]*N[0][2]+c23[2]*N[1][2]); //S_23 dot

  yp[15]=r[1][0]*y[15]-g[1]*m[1][0]*y[16]-sigma*y[18]+gamma*y[15]; //I_21 dot
  
  yp[16]=g[1]*y[16]-r[1][0]*y[15]-r[1][2]*y[17]-sigma*y[19]+gamma*y[16]; //I_22 dot

  yp[17]=r[1][2]*y[17]-g[1]*m[1][2]*y[16]-sigma*y[20]+gamma*y[17]; //I_23 dot

  yp[18]=r[1][0]*y[18]-g[1]*m[1][0]*y[19]-betainn*y[12]*(y[3]+y[15]+y[27])/(Np[0])
	  -c21[0]*betatra[0]*y[12]*(c21[0]*y[15]+c31[0]*y[27])/(c21[0]*N[1][0]+c31[0]*N[2][0])
	  -c21[1]*betatra[1]*y[12]*(c21[1]*y[15]+c31[1]*y[27])/(c21[1]*N[1][0]+c31[1]*N[2][0])
	  -c21[2]*betatra[2]*y[12]*(c21[2]*y[15]+c31[2]*y[27])/(c21[2]*N[1][0]+c31[2]*N[2][0])+sigma*y[18]; //E_21 dot

  yp[19]=g[1]*y[19]-r[1][0]*y[18]-r[1][2]*y[20]-betainn*(y[13]*(y[4]+y[16]+y[28]))/(Np[1])+sigma*y[19]; //E_22 dot

  yp[20]=r[1][2]*y[20]-g[1]*m[1][2]*y[19]-betainn*y[14]*(y[5]+y[17]+y[29])/(Np[2]) 
	  -c23[0]*betatra[0]*y[12]*(c13[0]*y[5]+c23[0]*y[17])/(c13[0]*N[0][2]+c23[0]*N[1][2])
	  -c23[1]*betatra[1]*y[12]*(c13[1]*y[5]+c23[1]*y[17])/(c13[1]*N[0][2]+c23[1]*N[1][2])
	  -c23[2]*betatra[2]*y[12]*(c13[2]*y[5]+c23[2]*y[17])/(c13[2]*N[0][2]+c23[2]*N[1][2])+sigma*y[20]; //E_23 dot

  yp[21]=r[1][0]*y[21]-g[1]*m[1][0]*y[22]-gamma*y[15]; //R_21 dot

  yp[22]=g[1]*y[22]-r[1][0]*y[21]-r[1][2]*y[23]-gamma*y[16]; //R_22 dot

  yp[23]=r[1][2]*y[23]-g[1]*m[1][2]*y[22]-gamma*y[17]; //R_23 dot

  yp[24]=r[2][0]*y[24]-g[2]*m[2][0]*y[26]+betainn*y[24]*(y[3]+y[15]+y[27])/(Np[0])
	  +c31[0]*betatra[0]*y[24]*(c31[0]*y[27]+c32[0]*y[28]+c21[0]*y[15])/(c31[0]*N[2][0]+c32[0]*N[2][1]+c21[0]*N[1][0])
	  +c31[1]*betatra[1]*y[24]*(c31[1]*y[27]+c32[1]*y[28]+c21[1]*y[15])/(c31[1]*N[2][0]+c32[1]*N[2][1]+c21[1]*N[1][0])
	  +c31[2]*betatra[2]*y[24]*(c31[2]*y[27]+c32[2]*y[28]+c21[2]*y[15])/(c31[2]*N[2][0]+c32[2]*N[2][1]+c21[2]*N[1][0]); //S_31 dot

  yp[25]=r[2][1]*y[25]-g[2]*m[2][1]*y[26]+betainn*y[25]*(y[4]+y[16]+y[28])/(Np[1])
	  +c32[0]*betatra[0]*y[25]*(c32[0]*y[28]+c31[0]*y[27])/(c32[0]*N[2][1]+c31[0]*N[2][0])
	  +c32[1]*betatra[1]*y[25]*(c32[1]*y[28]+c31[1]*y[27])/(c32[1]*N[2][1]+c31[1]*N[2][0])
	  +c32[2]*betatra[2]*y[25]*(c32[2]*y[28]+c31[2]*y[27])/(c32[2]*N[2][1]+c31[2]*N[2][0]); //S_32 dot

  yp[26]=g[2]*y[26]-r[2][0]*y[24]-r[2][1]*y[25]+betainn*y[26]*(y[5]+y[17]+y[29])/(Np[2]);//S_33 dot

  yp[27]=r[2][0]*y[27]-g[2]*m[2][0]*y[29]-sigma*y[30]+gamma*y[27]; //I_31 dot

  yp[28]=r[2][1]*y[28]-g[2]*m[2][1]*y[29]-sigma*y[31]+gamma*y[28]; //I_32 dot

  yp[29]=g[2]*y[29]-r[2][0]*y[27]-r[2][1]*y[28]-sigma*y[32]+gamma*y[29]; //I_33 dot

  yp[30]=r[2][0]*y[30]+g[2]*m[2][0]*y[32]-betainn*y[24]*(y[3]+y[15]+y[27])/(Np[0])
	  -c31[0]*betatra[0]*y[24]*(c31[0]*y[27]+c32[0]*y[28]+c21[0]*y[15])/(c31[0]*N[2][0]+c32[0]*N[2][1]+c21[0]*N[1][0])
	  -c31[1]*betatra[1]*y[24]*(c31[1]*y[27]+c32[1]*y[28]+c21[1]*y[15])/(c31[1]*N[2][0]+c32[1]*N[2][1]+c21[1]*N[1][0])
	  -c31[2]*betatra[2]*y[24]*(c31[2]*y[27]+c32[2]*y[28]+c21[2]*y[15])/(c31[2]*N[2][0]+c32[2]*N[2][1]+c21[2]*N[1][0])+sigma*y[30]; //E_31 dot

  yp[31]=r[2][1]*y[31]+g[2]*m[2][1]*y[32]-betainn*y[25]*(y[4]+y[16]+y[28])/(Np[1])
	  -c32[0]*betatra[0]*y[25]*(c32[0]*y[28]+c31[0]*y[27])/(c32[0]*N[2][1]+c31[0]*N[2][0])
	  -c32[1]*betatra[1]*y[25]*(c32[1]*y[28]+c31[1]*y[27])/(c32[1]*N[2][1]+c31[1]*N[2][0])
	  -c32[2]*betatra[2]*y[25]*(c32[2]*y[28]+c31[2]*y[27])/(c32[2]*N[2][1]+c31[2]*N[2][0])+sigma*y[31];//E_32 dot

  yp[32]=g[2]*y[32]-r[2][0]*y[30]-r[2][1]*y[31]-betainn*y[26]*(y[5]+y[17]+y[29])/(Np[2])+sigma*y[32]; //E_33 dot

  yp[33]=r[2][0]*y[33]-g[2]*m[2][0]*y[35]-gamma*y[27];//R_31 dot

  yp[34]=r[2][1]*y[34]-g[2]*m[2][1]*y[34]-gamma*y[28];//R_32 dot

  yp[35]=g[2]*y[35]-r[2][0]*y[33]-r[2][1]*y[34]-gamma*y[29];//R_33 dot

  }

int main ( )
{
 const int n=36;
 interval t=0.0, tend=20.0;
 iVector y(n); 
 y[0]=; //S_11
 y[1]=; //S_12
 y[2]=; //S_13
 y[3]=0; //I_11
 y[4]=0; //I_12
 y[5]=0; //I_13
 y[6]=0; //E_11
 y[7]=0; //E_12
 y[8]=0; //E_13
 y[9]=;  //R_11
 y[10]=; //R_12
 y[11]=; //R_13
 
 y[12]=; //S_21
 y[13]=; //S_22
 y[14]=; //S_23
 y[15]=0; //I_21
 y[16]=0; //I_22
 y[17]=0; //I_23
 y[18]=0; //E_21
 y[19]=0; //E_22
 y[20]=0; //E_23
 y[21]=;  //R_21
 y[22]=; //R_22
 y[23]=; //R_23

 y[24]=; //S_21
 y[25]=; //S_22
 y[26]=; //S_23
 y[27]=0; //I_21
 y[28]=0; //I_22
 y[29]=0; //I_23
 y[30]=0; //E_21
 y[31]=0; //E_22
 y[32]=0; //E_23
 y[33]=;  //R_21
 y[34]=; //R_22
 y[35]=; //R_23

AD*ad= new FADBAD_AD(n,SIER,SIER);
VNODE*Solver= new VNODE(ad);

Solver->integrate(t,y,tend);
if ( ! Solver->successful())
 cout<<¡±VNODE-LP could not r each t = ¡±<<tend<<endl;

 cout<<¡±Solution enclosure at t = ¡±<<t<<endl;
 printVector(y);
return 0 ;
}