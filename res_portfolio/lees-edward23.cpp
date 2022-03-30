// 0114時間積算のleesedwardに書き換え+sllod法的に運動方程式を変える



#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <cmath>
#include <omp.h>
#include <time.h>
#include <iostream>
#include<fstream>
#include <numeric>
#include "MT.h"
#include <vector>
#include <tuple>
// #include <direct.h>
#include <sys/stat.h>
#include <thread>

using std::this_thread::sleep_for;
using namespace std;

// namespace fs = std::filesystem;
// using std::cout; using std::endl;
// using std::system; using std::string;
// namespace fs = std::filesystem;
#define PI 3.14159265359

// virial.cppにleesedwardを適用させた

// 関数プロトタイプ宣言


void makevorticity_map(
    vector<vector<double>> &Vorticity_MAP,vector<vector<double>> &Vx_MAP,vector<vector<double>> &Vy_MAP,
vector<double> Rx,vector<double> Ry,vector<double> Vx,vector<double> Vy,double vor_l_g,int vor_ngx,int vor_ngy,
int NP);


void initialize_ra(double ave_flow,double temp,double m,double dx,double dy,int Nx,int Ny ,double lx,double ly,
vector<double> &rx,vector<double> &ry,vector<double> &ax,vector<double> &ay,
vector<double> &vx,vector<double> &vy);



double  Maxwell_velocity(double temp,double m);



void force3(double &Pot,int NP,
vector<double> &Ax,vector<double> &Ay,vector<double> Rx,vector<double> Ry,
double  Lx,double Ly,double zeta_M,double K_M,double ave_flow,
	    vector<double> Vx,vector<double> Vy,double &w,double gh,double h,double shear_rate,int time);


////////////////////////////////////////////////////////////

int main(int argc,char *argv[]){
    // const double rho=0.80;
    const double rho=stod(argv[1]);
    const double dx=sqrt(1.0/rho);
    const double dy=dx;
    const int Nx=40;
    const int Ny=40;
    // const double shear_rate=0.4;
    const double shear_rate=stod(argv[2]);
    const string ver=to_string(rho).substr(0,4)+"_"+to_string(shear_rate).substr(0,6)+"_"+to_string(Nx);
 

    cout<<ver<<endl;
    const double lx=Nx*dx;
    const double ly=Ny*dy;
    const double Volume=lx*ly;
    // const int Nx=int(lx/dx);
    // const int Ny=int(ly/dy);
    int N;
    const double temp=1.0;
    const double m=1.0;
    const double ave_flow=0.0;
    // const double epsilon=10.0;


    
    double pot=0.0;
    double kin=0.0;
    double pot_ij=0.0;
    double total_kin=0.0;
    double total_e=0.0;

    const double h=5e-4;
    const double nsteps=3e5;
    const double nout_e=nsteps/10;
    const double t_end=nsteps*h;
    const double h_rev=1/h ;
    const double h2=0.5*h*h;


    // 変数初期値
    vector<double> rx;
    vector<double> ry;
    vector<double> vx;    
    vector<double> vy;
    vector<double> ax;
    vector<double> ay;
        // ランジュバンパラメータ
    const double aa=0.5;
    const double eta=1.0;
    const double kb=1.0;
    const double zeta=6.0*PI*aa*eta;
    const double zeta_M=zeta/m;
    const double K=2.0*kb*temp/zeta;
    const double K_M=sqrt(K)/m;


    // ビリアル
    double w;
    double nktv;
    double p;
    double mvxvy;


    // lees-edward
    // const double shear_rate=0.4;
    const double gh=shear_rate*ly;
    // const int le_ngx=ceil(gh/l_gx_d);

    const double vor_l_g=2.0;
    const int vor_ngx=ceil(lx/vor_l_g);
    const int vor_ngy=ceil(ly/vor_l_g);
    vector<vector<double>> vorticity_MAP(vor_ngy,vector<double>(vor_ngx));
    vector<vector<double>> vx_MAP(vor_ngy,vector<double>(vor_ngx));
    vector<vector<double>> vy_MAP(vor_ngy,vector<double>(vor_ngx));


    // const char *rx_file="./lee2-rx12.dat";
    // const char *ry_file="./lee2-ry12.dat";
    // const char *vx_file="./lee2-vx12.dat";
    // const char *vy_file="./lee2-vy12.dat";
    // const char *vor_file="./lee2-vor12.dat";
    // const char *data_file="./lee2-data12.dat";
    // const char *p_file="./lee2-p12.dat";

    // std::filesystem::create_directory(ver);
    const char *dir_name=ver.c_str();
    mkdir(dir_name,0755);
    // mkdir(dir_name);
    const string rx_file="./"+ver+"/lee-rx.dat";
    const string ry_file="./"+ver+"/lee-ry.dat";
    const string vx_file="./"+ver+"/lee-vx.dat";
    const string vy_file="./"+ver+"/lee-vy.dat";
    const string vor_file="./"+ver+"/lee-vor.dat";
    const string data_file="./"+ver+"/lee-data.dat";
    const string p_file="./"+ver+"/lee-p.dat";
    const string K_file="./"+ver+"/lee-K.dat";
    const string Pot_file="./"+ver+"/lee-Pot.dat";


    ofstream rxofs(rx_file);
    ofstream ryofs(ry_file);
    ofstream vxofs(vx_file);
    ofstream vyofs(vy_file);
    ofstream vorofs(vor_file);
    ofstream dataofs(data_file);
    ofstream pofs(p_file);
    ofstream kofs(K_file);
    ofstream potofs(Pot_file);

    time_t start_time,end_time;

    start_time=time(NULL);

    int cy;
    double cx;
    double KE;
    double ratio;
    double nowT;



    initialize_ra(ave_flow,temp,m,dx,dy,Nx,Ny ,lx,ly,
rx,ry,ax,ay,vx,vy);

    N=rx.size();
    
    dataofs<<N<<" "<<vor_ngx<<" "<<vor_ngy<<endl;
    dataofs<<lx<<" "<<ly<<" "<<shear_rate<<" "<<endl;

    force3(pot,N,ax,ay,rx,ry,lx,ly,zeta_M,K_M,ave_flow,vx,vy,w,gh,h,shear_rate,0);
    nktv=N*kb*temp/Volume;
    // cout<<gh<<endl;
   // メインループ
    for(int t=1;t<nsteps;t++){
        cout<<t<<endl;
        for( int i=0;i<N;i++){
            // sllod法
            rx[i]=rx[i]+vx[i]*h+(ax[i]*h*h)/2+shear_rate*ry[i]*h;
        // normal
            // rx[i]=rx[i]+vx[i]*h+(ax[i]*h*h)/2;
            
            ry[i]=ry[i]+vy[i]*h+(ay[i]*h*h)/2;
            vx[i]=vx[i]+ax[i]*h/2;
            vy[i]=vy[i]+ay[i]*h/2;

        //  lees-edward境界
            cy=round(ry[i]/ly);
            cx=rx[i]-cy*gh*h*t;
            rx[i]=cx-round(cx/lx)*lx;
	    //            vx[i]=vx[i]-cy*gh;
            ry[i]=ry[i]-round(ry[i]/ly)*ly;

        }
        force3(pot,N,ax,ay,rx,ry,lx,ly,zeta_M,K_M,ave_flow,vx,vy,w,gh,h,shear_rate,t);
        // mvxvy=0.0;
        for (int i=0;i<N;i++){
            vx[i]=vx[i]+ax[i]*h/2;
            vy[i]=vy[i]+ay[i]*h/2;
            // mvxvy+=vx[i]*vy[i];
        }


// // 速度スケーリング/////////////////////
         nowT=0.0;
         for (int i=0;i<N;i++){
             nowT+=0.5*(vx[i]*vx[i]+vy[i]*vy[i]);
         }
         nowT/=N;
         ratio=sqrt(temp/nowT);
        
    
         for (int i=0;i<N;i++){
             vx[i]=vx[i]*ratio;
             vy[i]=vy[i]*ratio;
         }
//         /////////////////////////////////////////


        force3(pot,N,ax,ay,rx,ry,lx,ly,zeta_M,K_M,ave_flow,vx,vy,w,gh,h,shear_rate,t);


        mvxvy=0.0;
        for (int i=0;i<N;i++){
            mvxvy+=(vx[i]+shear_rate*ry[i])*vy[i];
        }


        if((t>nsteps*0.5)&&(t%1==0)){
            pofs<<(mvxvy+w)/Volume<<" ";
        }


        if((t>nsteps*0.0)&&(t%10==0)){
            KE=0;

            for(int i=0 ;i<vx.size();i++){
                KE+=vx[i]*vx[i]+vy[i]*vy[i];
            }

            kofs<<0.5*KE/vx.size()<<" ";
            potofs<<pot<<" ";
        }


        if((t>nsteps*0.8)&&(t%100==0)){
            for(size_t i=0;i<rx.size();i++){
                rxofs<<rx[i]<<" ";
                ryofs<<ry[i]<<" ";
                vxofs<<vx[i]<<" ";
                vyofs<<vy[i]<<" ";
            }
            rxofs<<endl;
            ryofs<<endl;
            vxofs<<endl;
            vyofs<<endl;
        }
    }
    rxofs.close();
    ryofs.close();
    vxofs.close();
    vyofs.close();
    // vorofs.close();
    dataofs.close();
    pofs.close();

    end_time=time(NULL);
    cout<<end_time-start_time<<endl;
    return 0;
}




void makevorticity_map(
    vector<vector<double>> &Vorticity_MAP,vector<vector<double>> &Vx_MAP,vector<vector<double>> &Vy_MAP,
vector<double> Rx,vector<double> Ry,vector<double> Vx,vector<double> Vy,double vor_l_g,int vor_ngx,int vor_ngy,
int NP
){
    int vor_gx_map,vor_gy_map;
    // 初期化
    Vx_MAP.assign(vor_ngy,vector<double> (vor_ngx,0));
    Vy_MAP.assign(vor_ngy,vector<double> (vor_ngx,0));
    Vorticity_MAP.assign(vor_ngy,vector<double> (vor_ngx,0));

    for (int i=0;i<NP;i++){
        vor_gx_map=int(Rx[i]/vor_l_g);
        vor_gy_map=int(Ry[i]/vor_l_g);
        Vx_MAP[vor_gy_map][vor_gx_map]+=Vx[i];
        Vy_MAP[vor_gy_map][vor_gx_map]+=Vy[i];
    }
    for (int i=1;i<vor_ngy-1;i++){
        for (int j=1;j<vor_ngx-1;j++){
            Vorticity_MAP[i][j]=(Vy_MAP[i+1][j]-Vy_MAP[i-1][j]-Vx_MAP[i][j+1]+Vx_MAP[i][j-1])/2;
        }
    }

}



void force3(double &Pot,int NP,
vector<double> &Ax,vector<double> &Ay,vector<double> Rx,vector<double> Ry,
double  Lx,double Ly,double zeta_M,double K_M,double ave_flow,
	    vector<double> Vx,vector<double> Vy,double &w,double gh,double h,double shear_rate,int time){

    double Pot_ij;
    double rxij;
    double ryij;
    double temp_rxij;
    double temp_ryij;
    double r2;
    double ir2,ir6,ir12;
    const double sigma=1.0;
    double sigma6=1.0;
    double sigma12=1.0;
    const double rc=pow(2.0,1.0/6.0)*sigma;
    double rc2=rc*rc;
    double epsilon=1.0;
    double fx,fy;
    double wij;
    double rxfy;
    int cy;
    double cdelx;
    double fdotp;
    double gpxpy;
    double sump2;
    double alpha;

    Ax.assign(NP,0.0);
    Ay.assign(NP,0.0);
    
    Pot=0.0;
    Pot_ij=0.0;
    w=0.0;
    
    // cout<<rc2<<endl;

    for(int i=0;i<NP-1;i++){
        // cout<<i<<endl;
        for(int j=i+1;j<NP;j++){
 // minimum image convention
            rxij=Rx[i]-Rx[j];
            ryij=Ry[i]-Ry[j];

            cy=round(ryij/Ly);
            cdelx=rxij-cy*gh*h*time;
            rxij=cdelx-round(cdelx/Lx)*Lx;
            ryij=ryij-cy*Ly;


            r2=rxij*rxij+ryij*ryij;
            // cout<<r2<<endl;
            ir2=1.0/r2;
            ir6=ir2*ir2*ir2;
            ir12=ir6*ir6;




            if(r2>=rc2){
                fx=0.0;
                fy=0.0;
                Pot_ij=0.0;
		        rxfy=0.0;
            }
            else{
                wij=24.0*epsilon*(2.0*sigma12*ir12-sigma6*ir6);
                fx=wij*ir2*rxij;
                fy=wij*ir2*ryij;
                rxfy=(rxij*fy+ryij*fx)/2;
                // cout<<fx<<endl;
                // fx=24.0*epsilon*(2.0*sigma12*ir12-sigma6*ir6)*ir2*rxij;
                // fy=24.0*epsilon*(2.0*sigma12*ir12-sigma6*ir6)*ir2*ryij;
                Pot_ij=4.0*epsilon*(sigma12*ir12-sigma6*ir6)+epsilon;	    
            }
            Ax[i]+=fx;
            Ay[i]+=fy;
            Ax[j]-=fx;
            Ay[j]-=fy;
            Pot+=Pot_ij;
            // w+=wij;
            w+=rxfy;
        }
    }




    // gpxpy=0.0;
    // sump2=0.0;
    // fdotp=0.0;
    // alpha=0.0;
    //  for (int i=0;i<NP;i++){
    // fdotp+=Ax[i]*Vx[i]+Ay[i]*Vy[i];
    //sump2+=Vx[i]*Vx[i]+Vy[i]*Vy[i];
    // gpxpy+=shear_rate*Vx[i]*Vy[i];
    // }

    // alpha=(fdotp-gpxpy)/sump2;
    // alpha=(fdotp-gpxpy)/sump2;
    // if(abs(alpha)>1.0){
    //     cout<<"alpha"<<alpha<<endl;
    // }
    



// Sllod法 
    for(int i=0;i<NP;i++){
         Ax[i]=Ax[i]-shear_rate*Vy[i];
        // gausu熱浴ver
      // Ax[i]=Ax[i]-shear_rate*Vy[i]-alpha*Vx[i];
	//	    Ay[i]=Ay[i]-alpha*Vy[i];
    }

    

// // Lanuven 
//     for(int i=0;i<NP;i++){
//         Ax[i]=Ax[i]-zeta_M*Vx[i]+K_M*genrand_real1();
//         Ay[i]=Ay[i]-zeta_M*Vy[i]+K_M*genrand_real1();
//     }
   




}

void initialize_ra(double ave_flow,double temp,double m,double dx,double dy,int Nx,int Ny ,double lx,double ly,
vector<double> &rx,vector<double> &ry,vector<double> &ax,vector<double> &ay,
vector<double> &vx,vector<double> &vy){
    double vx_sum=0.0;
    double vy_sum=0.0;
    double kin0=0.0;
    int ave_vx,ave_vy;

    for(int i=0;i<Nx/2;i++){
        for (int j=0;j<Ny/2;j++){


            rx.push_back(i*dx+dx/2);
            ry.push_back(j*dy+dy/2);
            vx.push_back(Maxwell_velocity(temp,m));
            vy.push_back(Maxwell_velocity(temp,m));
            ax.push_back(0.0);
            ay.push_back(0.0);

            rx.push_back(-i*dx-dx/2);
            ry.push_back(j*dy+dy/2);
            vx.push_back(Maxwell_velocity(temp,m));
            vy.push_back(Maxwell_velocity(temp,m));
            ax.push_back(0.0);
            ay.push_back(0.0);

            rx.push_back(-i*dx-dx/2);
            ry.push_back(-j*dy-dy/2);
            vx.push_back(Maxwell_velocity(temp,m));
            vy.push_back(Maxwell_velocity(temp,m));
            ax.push_back(0.0);
            ay.push_back(0.0);

            rx.push_back(i*dx+dx/2);
            ry.push_back(-j*dy-dy/2);
            vx.push_back(Maxwell_velocity(temp,m));
            vy.push_back(Maxwell_velocity(temp,m));
            ax.push_back(0.0);
            ay.push_back(0.0);


        }
    }

    for (int i=0;i<vx.size();i++){
        vx_sum=vx_sum+vx[i];
        vy_sum=vy_sum+vy[i];
    }

    ave_vx=vx_sum/(vx.size());
    ave_vy=vy_sum/(vy.size());

    for (int i=0;i<vx.size();i++){
        vx[i]=vx[i]-ave_vx;
        // vx[i]=vx[i]+ave_flow;
        vy[i]=vy[i]-ave_vy;
        kin0=kin0+(pow(vx[i],2)+pow(vy[i],2));
    }
    kin0=kin0/2;
}



double  Maxwell_velocity(double temp,double m){
    double rand1=genrand_real1();
    double rand2=genrand_real1();
    double velocity;

    velocity=(sqrt(-2*(temp/m)*log(rand1)))*cos(2*PI*rand2);
    return velocity;
}


