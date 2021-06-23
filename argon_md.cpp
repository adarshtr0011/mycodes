#include<iostream>
#include<cstdlib>
#include<math.h>
using namespace std;

//constants
#define sigma 3.76e-10                    //sigma of Argon atom
#define eps 1.6556e-21                    //epsilon of Argon in J/atom
#define mAr 6.6335209e-26                 //mass of Argon atom in kg
#define fms 1e-15;                        //femtoseconds

 
int main()
{
    //constant C1 and C2 used in calculating potential energy,    C3 and C4 used in calculating acceleration of atoms,   m= 0.5 times mass of Argon atom
    const double C1=(4*eps*pow(sigma,12)), C2=(4*eps*pow(sigma,6)),  C3=(48*eps*pow(sigma,12))/mAr, C4=(24*eps*pow(sigma,6))/mAr, m= mAr/2;

    cout<<C1<<endl<<C2<<endl<<C3<<endl<<C4<<endl;    
    //dT is in femtoseconds, ns is number of timesteps
    double dT;
    int ns;
    cout<<"enter time step in femto seconds and number of time steps"<<endl;
    cin>>dT>>ns;

    //converting femtoseconds to seconds
    double dt= dT * fms;
    cout<<dt<<endl;

    
   //variable declarations, initialized matrices of size ns*64, 64 atoms with 'ns' number of timesteps
    double vx[ns][64],vy[ns][64], vz[ns][64];
    double x[ns][64], y[ns][64], z[ns][64], ax[ns][64], ay[ns][64], az[ns][64], P[ns][64], K[ns][64], r;
    
    
    int cnt=0;


    //loop to iniitialize our required 2D arrays as 0
    for(int i=0; i<ns; i++)
    {
    for(int j=0; j<64; j++)
     {
      ax[i][j]=0;
      ay[i][j]=0;
      az[i][j]=0;
      P[i][j]=0;
      K[i][j]=0;
     }
    
    }

    cout<<"t1"<<endl;


    
    //loop to generate random velocities ranging from -500 to 500 m/s in each axis using rand() function
    for(int i=0; i<64; i++)
    {
         
        vx[0][i] = -500 + (rand() % 1001);
        vy[0][i] = -500 + (rand() % 1001);
        vz[0][i] = -500 + (rand() % 1001);

    }

cout<<"t2"<<endl;
 
//main loop going through ns amount of timesteps
for (int i=0; i<ns ; i++)
 {
    //initial case 0, that is first cycle
    if(i==0)
    {
      { //loop for initial positions of 64 atoms
         for(int j=0; j<64; j+=2)
         for(int k=0; k<64; k+=2)
          for(int l=0; l<64; l+=2)
          {
           x[i][cnt]=j*sigma;
           y[i][cnt]=k*sigma;
           z[i][cnt]=l*sigma;
           cnt++;
          }
                 
      }

     cout<<"t"<<endl;

     //loop to find initial accelerations of atoms 

     for(int j=0; j<64; j++)
     {
     //calc Kinetic energy of each atom initially
      K[i][j]+= m*(pow(vx[i][j],2) + pow(vx[i][j],2) + pow(vx[i][j],2) );
      cout<<"t3"<<endl;
     

      for(int k=0; k<64; k++)
       
       {
         //finding distance between coordinates
         r= sqrt( (pow((x[i][j]-x[i][k]),2)) + (pow((y[i][j]-y[i][k]),2)) + (pow((z[i][j]-z[i][k]),2)) );

        //calc constant factor in Leonard Jones Potential force
         double temp1= (C3/pow(r,14)) - (C4/pow(r,8));

        //finding acceleration of each atom initially
         ax[i][j]+= temp1 * (x[i][j]-x[i][k]) ;
         ay[i][j]+= temp1 * (y[i][j]-y[i][k]) ;
         az[i][j]+= temp1 * (z[i][j]-z[i][k]) ;

        //calc potential energy of each atom initially
         P[i][j]+=  (C1/pow(r,12)) - (C2/pow(r,6));
       }
     }
    }
     
   
  //for case other than 0, that is from 2nd cycle onwards 
  else

  {
    //calc the next position using verlet algorithm
    for(int j=0; j<64; j++)
    {
        x[i][j]= x[i-1][j] + vx[i-1][j]*dt + 0.5* ax[i-1][j]*dt*dt; 
        y[i][j]= y[i-1][j] + vy[i-1][j]*dt + 0.5* ay[i-1][j]*dt*dt; 
        z[i][j]= z[i-1][j] + vz[i-1][j]*dt + 0.5* az[i-1][j]*dt*dt; 
    }
    cout<<"t5"<<endl;
    //calc the next accelerations and potential energy using the new positions obtained from above 
    for(int j=0; j<64; j++)
    {
        for(int k=0; k<64; k++)
        {
           r= sqrt( (pow((x[i][j]-x[i][k]),2)) + (pow((y[i][j]-y[i][k]),2)) + (pow((z[i][j]-z[i][k]),2)) );

           double temp1= (C3/pow(r,14)) - (C4/pow(r,8));

         ax[i][j]+= temp1 * (x[i][j]-x[i][k]) ;
         ay[i][j]+= temp1 * (y[i][j]-y[i][k]) ;
         az[i][j]+= temp1 * (z[i][j]-z[i][k]) ;

         P[i][j]+=  (C1/pow(r,12)) - (C2/pow(r,6));  

        }
    }

    cout<<"t6"<<endl;

    //calc the new velocities for this cycle and calculating Kinetic energy of each atom
    for(int j=0; j<64; j++)
    {
        //new velocities using Verlet algorithm
        vx[i][j] = vx[i-1][j] + 0.5* ( ax[i][j] + ax[i-1][j]) * dt;
        vx[i][j] = vx[i-1][j] + 0.5* ( ax[i][j] + ax[i-1][j]) * dt;
        vx[i][j] = vx[i-1][j] + 0.5* ( ax[i][j] + ax[i-1][j]) * dt;

        //new Kinetic energy calculation
        K[i][j]+= m*( pow(vx[i][j],2) + pow(vx[i][j],2) + pow(vx[i][j],2) );
    }
    cout<<"t8"<<endl;


 }

 }
cout<<"hey";
 return 0;
}
  




   
 
 
    

 


   
     
        
     
    
