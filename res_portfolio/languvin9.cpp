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
#include <sstream>

using namespace std;
#define PI 3.14159265359


// 関数宣言////////////////////////////////////////////////////////
void makeNeighborList(int aroundGridNum, vector<int> &neighbor_row, vector<int> &neighbor_col);

double  Maxwell_velocity(double temp,double m);


void initialize_ra(double ave_flow,double temp,double m,int remove_dx_num,double dx,double dy,int Nx,int Ny ,double lx,double ly,
vector<double> &rx,vector<double> &ry,vector<double> &ax,vector<double> &ay,
vector<double> &vx,vector<double> &vy);


void make_pairlist(int move_pair_length,int N_GX,int N_GY,int select_gx,int select_gy,int NEIGHBOR_LEN,
vector<int>NEIGHBOR_ROW,
vector<int> NEIGHBOR_COL,
vector<vector<int>> G_MAP,
vector<vector<int>> &move_PAIRLIST,
vector<int> &static_PAIRLIST);


void gmap_create4(double L_GX,double L_GY ,int static_pair_length,int move_pair_length,int NP,
int N_GX,int N_GY,
vector<vector<int>> &G_MAP,
vector<vector<int>> &move_pairlist,
vector<int> &static_pairlist,
vector<double> RX,
vector<double> RY,
vector<int> neighbor_row_static,
vector<int> neighbor_col_static,
vector<int> neighbor_row_move,
vector<int> neighbor_col_move,int neighbor_len_static,int neighbor_len_move);


void force(int move_pair_length,double &Pot,int NP,
vector<double> &Ax,vector<double> &Ay,
vector<vector<int>> move_PAIRLIST,
vector<int>static_PAIRLIST,
vector<double> Rx,
vector<double> Ry,
double  Lx,double Ly,
double static_dia,
double zeta_M,
double K_M,
double ave_flow,
vector<double> Vx,
vector<double> Vy,
double &staticfx,
double &staticfy
);


int BookKeeping(int NP,vector<double> Vx,vector<double> Vy,int aroundGridNum,
double h,double gridWidth,double rc);



void makevorticity_map(
    vector<vector<double>> &Vorticity_MAP,vector<vector<double>> &Vx_MAP,vector<vector<double>> &Vy_MAP,
vector<double> Rx,vector<double> Ry,vector<double> Vx,vector<double> Vy,double vor_l_g,int vor_ngx,int vor_ngy,
int NP);

///////////////////////////////////////////////////////////////

int main(){

    // int ver=3;
    const double static_dia=40.0;
    const double rho=0.20;
    const double dx=sqrt(1.0/rho);
    const double dy=dx;
    const double lx=static_dia*5*2.5;
    const double ly=static_dia*5;
    const int Nx=int(lx/dx);
    const int Ny=int(ly/dy);
    int N;
    const double temp=1.0;
    const double m=1.0;
    const double ave_flow=1.0;
    // const double epsilon=10.0;
    int reusableCount=0;
    const string ver="T=1_"+to_string(rho).substr(0,4)+"_"+to_string(ave_flow).substr(0,3)+"_"+to_string(static_dia).substr(0,2);
    int reusableCount_static;
    int reusableCount_move;


    cout<<ver<<endl;

   // 配列の宣言
    double pot=0.0;
    double kin=0.0;
    double pot_ij=0.0;
    double total_kin=0.0;
    double total_e=0.0;


    const double h=5e-4;
    const double nsteps=1.0e5;
    const double nout_e=nsteps/10;
    const double t_end=nsteps*h;
    const double h_rev=1/h ;
    const double h2=0.5*h*h;


    double before_vx;
    double before_vy;
    double after_vx;
    double after_vy;

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
    const double bo=(lx/0.4);
    const double l_gx_d=lx/bo;
    const double l_gx=l_gx_d;
    const double l_gy=l_gx;
    const int n_gx=ceil(lx/l_gx);
    const int n_gy=ceil(ly/l_gy);
    const int n_g_all=n_gx*n_gy;
    const int remove_dx_num=ceil(static_dia*0.5/dx);

    const double static_sigma =(static_dia*0.5)+0.5;
    const int aroundGridNum_move=4;
    const int aroundGridNum_staic=ceil(4*static_sigma);
    const int neighbor_len_static=2*2*(aroundGridNum_staic*(aroundGridNum_staic+1));
    const int neighbor_len_move=2*(aroundGridNum_move*(aroundGridNum_move+1));
    const int move_pair_length=15;
    const int static_pair_length=4000;

    vector<int> neighbor_row_static(neighbor_len_static/2);
    vector<int> neighbor_col_static(neighbor_len_static/2);
    // vector<int> temp_neighbor_row_static(neighbor_len_static/2);
    // vector<int> temp_neighbor_col_static(neighbor_len_static/2);
    vector<int> reverse_neighbor_row_static(neighbor_len_static/2);
    vector<int> reverse_neighbor_col_static(neighbor_len_static/2);
    vector<int> neighbor_row_move(neighbor_len_move);
    vector<int> neighbor_col_move(neighbor_len_move);

    vector<vector<int>> g_map(n_gy,vector<int>(n_gx));
    vector<vector<int>> move_pairlist;
    vector<int> static_pairlist(static_pair_length);


    double staticfx;
    double staticfy;


    const double vor_l_g=2.0*static_dia/40;
    const int vor_ngx=ceil(lx/vor_l_g);
    const int vor_ngy=ceil(ly/vor_l_g);
    vector<vector<double>> vorticity_MAP(vor_ngy,vector<double>(vor_ngx));
    vector<vector<double>> vx_MAP(vor_ngy,vector<double>(vor_ngx));
    vector<vector<double>> vy_MAP(vor_ngy,vector<double>(vor_ngx));








    // const char *rx_file="./rx10.dat";
    // const char *ry_file="./ry10.dat";
    // const char *vx_file="./vx10.dat";
    // const char *vy_file="./vy10.dat";
    // const char *vor_file="./vor10.dat";
    // const char *data_file="./data10.dat";
    // const char *changevx_file="./changevx.dat";
    // const char *changevy_file="./changevx.dat";
    

    const string big_dir="languvin11";
    const char *big_dir_name=big_dir.c_str();
    //    mkdir(big_dir_name);
    mkdir(big_dir_name,0755);
    const string main_dir="./languvin11/"+ver;
    const char *dir_name=main_dir.c_str();
    //mkdir(dir_name);
    mkdir(dir_name,0755);
    const string rx_file=main_dir+"/lg11-rx2.dat";
    const string ry_file=main_dir+"/lg11-ry2.dat";
    const string vx_file=main_dir+"/lg11-vx2.dat";
    const string vy_file=main_dir+"/lg11-vy2.dat";
    const string vor_file=main_dir+"/lg11-vor2.dat";
    const string data_file=main_dir+"/lg11-data2.dat";
    const string p_file=main_dir+"/lg11-p2.dat";
    const string losvx_file=main_dir+"/lg11-losvx2.dat";
    const string losvy_file=main_dir+"/lg11-losvy2.dat";
    const string staticfx_file=main_dir+"/lg11-staticfx2.dat";
    const string staticfy_file=main_dir+"/lg11-staticfy2.dat";

    const string vxmap_file=main_dir+"/lg11-vxmap2.dat";
    const string vymap_file=main_dir+"/lg11-vymap2.dat";



    ofstream rxofs(rx_file);
    ofstream ryofs(ry_file);
    ofstream vxofs(vx_file);
    ofstream vyofs(vy_file);
    ofstream vorofs(vor_file);
    ofstream dataofs(data_file);
    ofstream losvx_ofs(losvx_file);
    ofstream losvy_ofs(losvy_file);
    ofstream staticfx_ofs(staticfx_file);
    ofstream staticfy_ofs(staticfy_file);



    ofstream vxmapofs(vxmap_file);
    ofstream vymapofs(vymap_file);



    time_t start_time,end_time;

    start_time=time(NULL);




    makeNeighborList(aroundGridNum_move, neighbor_row_move, neighbor_col_move);
    makeNeighborList(aroundGridNum_staic, neighbor_row_static, neighbor_col_static);
    
    for (size_t i=0;i<neighbor_row_static.size();i++){
        reverse_neighbor_row_static[i]=neighbor_row_static[i]*-1;
    }

    for (size_t i=0;i<neighbor_col_static.size();i++){
        reverse_neighbor_col_static[i]=neighbor_col_static[i]*-1;
    }


    // 結合
    neighbor_row_static.insert(neighbor_row_static.end(),reverse_neighbor_row_static.begin(),reverse_neighbor_row_static.end());
    neighbor_col_static.insert(neighbor_col_static.end(),reverse_neighbor_col_static.begin(),reverse_neighbor_col_static.end());
   
    initialize_ra(ave_flow,temp,m,remove_dx_num,dx,dy,Nx,Ny ,lx,ly,
rx,ry,ax,ay,vx,vy);

    // cout<<Nx*Ny<<endl;
    N=rx.size();
    // cout<<N<<endl;
    move_pairlist.resize(N,vector<int>(move_pair_length));
    

    dataofs<<N<<" "<<vor_ngx<<" "<<vor_ngy<<endl;
    gmap_create4(l_gx,l_gy,static_pair_length,move_pair_length,N,n_gx,n_gy,g_map,move_pairlist,
    static_pairlist,rx,ry,neighbor_row_static,
    neighbor_col_static,neighbor_row_move,neighbor_col_move,neighbor_len_static,neighbor_len_move);


    force(move_pair_length,pot,N,ax,ay,move_pairlist,static_pairlist,
    rx,ry,lx,ly,static_dia,zeta_M,K_M,ave_flow,vx,vy,staticfx,staticfy);



    // メインループ
    for(int t=0;t<nsteps;t++){
        cout<<t<<endl;


        if((t>nsteps*0.5)&&(t%10==0)){
            before_vx=0.0;
            before_vy=0.0;
            for (int i =0;i<static_pairlist.size();i++){
                if(static_pairlist[i]==-1){
                    break;
                }
                before_vx+=vx[static_pairlist[i]];
                before_vy+=vy[static_pairlist[i]];
            }
        }



        /////////// 速度ベルれ法１//////////////////////////////
        // 静止粒子には適用しないので1から


        

        for( int i=1;i<N;i++){
            rx[i]=rx[i]+vx[i]*h+(ax[i]*h*h)/2;
            ry[i]=ry[i]+vy[i]*h+(ay[i]*h*h)/2;
            vx[i]=vx[i]+ax[i]*h/2;
            vy[i]=vy[i]+ay[i]*h/2;
            // 周期境界
            if(rx[i]>lx){
                rx[i]-=lx;
            }
            else if(rx[i]<=0.0){
                rx[i]+=lx;
            }

            if(ry[i]>ly){
                ry[i]-=ly;
            }
            else if(ry[i]<=0.0){
                ry[i]+=ly;
            }
        }



        // ランジュバン熱浴のために二回ブックキーピングをするのは、ちょっと計算効率悪すぎるのでコメントアウト
        // // BookKeeping./////////////////////////////////////////
        // if(reusableCount==0){
        //     gmap_create4(l_gx,l_gy,static_pair_length,move_pair_length,N,n_gx,n_gy,g_map,move_pairlist,
        //     static_pairlist,rx,ry,neighbor_row_static,
        //     neighbor_col_static,neighbor_row_move,neighbor_col_move,neighbor_len_static,neighbor_len_move);
        //     // static
        //     // reusableCount_static=BookKeeping(N,vx,vy,aroundGridNum_staic,h,l_gx_d,pow(2.0,1.0/6.0)*static_sigma);
        //     // move
        //     reusableCount_move=BookKeeping(N,vx,vy,aroundGridNum_move,h,l_gx_d,pow(2.0,1.0/6.0)*1.0);
        //     reusableCount=reusableCount_move;

        //     // if(reusableCount_move>reusableCount_static){
        //     //     cout<<"static"<<endl;
        //     //     reusableCount=reusableCount_static;
        //     // }
        //     // else{
        //     //     reusableCount=reusableCount_move;
        //     //     cout<<"move"<<endl;
        //     // }
        // }
        // else{
        //     reusableCount-=1;
        // }
        // cout<<reusableCount<<endl;


// end_BookKeeping./////////////////////////////////////////

        force(move_pair_length,pot,N,ax,ay,move_pairlist,static_pairlist,
            rx,ry,lx,ly,static_dia,zeta_M,K_M,ave_flow,vx,vy,staticfx,staticfy);

        for (int i=1;i<N;i++){
            vx[i]=vx[i]+ax[i]*h/2;
            vy[i]=vy[i]+ay[i]*h/2;
        }


        if((t>nsteps*0.5)&&(t%10==0)){
            after_vx=0.0;
            after_vy=0.0;
            for (int i =0;i<static_pairlist.size();i++){
                if(static_pairlist[i]==-1){
                    break;
                }
                after_vx+=vx[static_pairlist[i]];
                after_vy+=vy[static_pairlist[i]];
            }
            losvx_ofs<<(after_vx-before_vx)/h<<" ";
            losvy_ofs<<(after_vy-before_vy)/h<<" ";
        }

        if((t>nsteps*0.5)&&(t%10==0)){
            staticfx_ofs<<staticfx<<" ";
            staticfy_ofs<<staticfy<<" ";
        }







        // ランジュバン熱浴のために二回ブックキーピングをするのは、ちょっと計算効率悪すぎるのでコメントアウト
        // BookKeeping./////////////////////////////////////////
        if(reusableCount==0){
            gmap_create4(l_gx,l_gy,static_pair_length,move_pair_length,N,n_gx,n_gy,g_map,move_pairlist,
            static_pairlist,rx,ry,neighbor_row_static,
            neighbor_col_static,neighbor_row_move,neighbor_col_move,neighbor_len_static,neighbor_len_move);
            // static
            // reusableCount_static=BookKeeping(N,vx,vy,aroundGridNum_staic,h,l_gx_d,pow(2.0,1.0/6.0)*static_sigma);
            // move
            reusableCount_move=BookKeeping(N,vx,vy,aroundGridNum_move,h,l_gx_d,pow(2.0,1.0/6.0)*1.0);
            reusableCount=reusableCount_move;
            // if(reusableCount_move>reusableCount_static){
            //     reusableCount=reusableCount_static;
            // }
            // else{
            //     reusableCount=reusableCount_move;
            // }
        }
        else{
            reusableCount-=1;
        }


// end_BookKeeping./////////////////////////////////////////

        force(move_pair_length,pot,N,ax,ay,move_pairlist,static_pairlist,
            rx,ry,lx,ly,static_dia,zeta_M,K_M,ave_flow,vx,vy,staticfx,staticfy);







        if((t>nsteps*0.8)&&(t%10==0)){

            makevorticity_map(vorticity_MAP,vx_MAP,vy_MAP,rx,ry,vx,vy,vor_l_g,vor_ngx,vor_ngy,N);
            for (int i=0;i<vorticity_MAP.size();i++){
                for (int j=0;j<vorticity_MAP.at(0).size();j++){
                    vorofs<<vorticity_MAP[i][j]<<" ";
                    // cout<<vorticity_MAP[i][j]<<endl;
                }
            }
            vorofs<<endl;


            for (int i=0;i<vx_MAP.size();i++){
	      for (int j=0;j<vx_MAP.at(0).size();j++){
		vxmapofs<<vx_MAP[i][j]<<" ";
		// cout<<vorticity_MAP[i][j]<<endl;
	      }
            }
            vxmapofs<<endl;

            for (int i=0;i<vy_MAP.size();i++){
	      for (int j=0;j<vy_MAP.at(0).size();j++){
		vymapofs<<vy_MAP[i][j]<<" ";
		// cout<<vorticity_MAP[i][j]<<endl;
	      }
            }
            vymapofs<<endl;
        }



	

	if((t>nsteps*0.8)&&(t%200==0)){
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
    vorofs.close();
    dataofs.close();

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
    vector<vector<int>> pcount_MAP(vor_ngy,vector<int>(vor_ngx));
    int p_count;




    // 初期化
    pcount_MAP.assign(vor_ngy,vector<int> (vor_ngx,0));
    Vx_MAP.assign(vor_ngy,vector<double> (vor_ngx,0));
    Vy_MAP.assign(vor_ngy,vector<double> (vor_ngx,0));
    Vorticity_MAP.assign(vor_ngy,vector<double> (vor_ngx,0));
    
    for (int i=0;i<NP;i++){
        vor_gx_map=int(Rx[i]/vor_l_g);
        vor_gy_map=int(Ry[i]/vor_l_g);
        pcount_MAP[vor_gy_map][vor_gx_map]+=1;
        Vx_MAP[vor_gy_map][vor_gx_map]+=Vx[i];
        Vy_MAP[vor_gy_map][vor_gx_map]+=Vy[i];
    }

    for (int i=0;i<vor_ngy;i++){
      for (int j=0;j<vor_ngx;j++){
            
	p_count=pcount_MAP[i][j];
	if(p_count==0){
	  p_count=1;
	}
	Vx_MAP[i][j]/=p_count;
	Vy_MAP[i][j]/=p_count; 
      }
    }

    for (int i=1;i<vor_ngy-1;i++){
        for (int j=1;j<vor_ngx-1;j++){
            Vorticity_MAP[i][j]=(Vy_MAP[i+1][j]-Vy_MAP[i-1][j]-Vx_MAP[i][j+1]+Vx_MAP[i][j-1])/2;
        }
    }

}


int BookKeeping(int NP,vector<double> Vx,vector<double> Vy,int aroundGridNum,
double h,double gridWidth,double rc){
    double maxCandidateV;
    double maxV=0.0;
    double tLim;
    int reusableNum;
    for (int i=0;i<NP;i++){
        maxCandidateV=sqrt(Vx[i]*Vx[i]+Vy[i]*Vy[i]);
        if(maxCandidateV>maxV){
            maxV=maxCandidateV;
        }
        tLim=(aroundGridNum*gridWidth-rc)/(2*maxV);
        reusableNum=int(tLim/(1.5*h));
    }
    return reusableNum;
}

void force(int move_pair_length,double &Pot,int NP,
vector<double> &Ax,vector<double> &Ay,
vector<vector<int>> move_PAIRLIST,
vector<int>static_PAIRLIST,
vector<double> Rx,
vector<double> Ry,
double  Lx,double Ly,
double static_dia,
double zeta_M,
double K_M,
double ave_flow,
vector<double> Vx,
vector<double> Vy,
double &staticfx,
double &staticfy
){
    double Pot_ij;
    int roop_num;
    int pair_index;
    double rxij;
    double ryij;
    double r2;
    double ir2,ir6,ir12;
    double sigma;
    double sigma6,sigma12;
    double rc,rc2;
    double epsilon;
    double fx,fy;
    
    
    Ax.assign(NP,0.0);
    Ay.assign(NP,0.0);
    staticfx=0.0;
    staticfy=0.0;
    
    // for (int i=0;i<NP;i++){
    //     Ax[i]=0.0;
    //     Ay[i]=0.0;
    // }
    Pot=0.0;
    Pot_ij=0.0;
    for (int i=0;i<NP;i++){

        if (i==0)
        {
            sigma=(static_dia/2)+0.5;
            sigma6=pow(sigma,6);
            sigma12=sigma6*sigma6;
            epsilon=10.0;
        }
        else{
            sigma=1.0;
            sigma6=1.0;
            sigma12=1.0;
            epsilon=1.0;
        }
        rc=pow(2.0,1.0/6.0)*sigma;
        rc2=rc*rc;
        
        roop_num=move_PAIRLIST[i][move_pair_length-1];
        for (int j=0;j<roop_num;j++){
            if(i==0){
                pair_index=static_PAIRLIST[j];
            }
            else{
                pair_index=move_PAIRLIST[i][j];
            }
            rxij=Rx[i]-Rx[pair_index];


            // minimum image convention
            if (rxij>=Lx/2){
                rxij=rxij-Lx;
            }
            else if (rxij<-Lx/2){
                rxij=rxij+Lx;
            }
            else{
                rxij=rxij;
            }
            ryij=Ry[i]-Ry[pair_index];
            if (ryij>=Ly/2){
                ryij=ryij-Ly;
            }
            else if (ryij<-Ly/2){
                ryij=ryij+Ly;
            }
            else{
                ryij=ryij;
            }
            r2=rxij*rxij+ryij*ryij;
            ir2=1.0/r2;
            ir6=ir2*ir2*ir2;
            ir12=ir6*ir6;
            
            if(r2>=rc2){
                fx=0.0;
                fy=0.0;
                Pot_ij=0.0;
            }
            else{
                fx=24.0*epsilon*(2.0*sigma12*ir12-sigma6*ir6)*ir2*rxij;
                fy=24.0*epsilon*(2.0*sigma12*ir12-sigma6*ir6)*ir2*ryij;
                Pot_ij=4.0*epsilon*(sigma12*ir12-sigma6*ir6)+epsilon;
            }
            Ax[i]+=fx;
            Ay[i]+=fy;
            Ax[pair_index]-=fx;
            Ay[pair_index]-=fy;
            Pot+=Pot_ij;
            if(i==0){
                staticfx+=fx;
                staticfy+=fy;
            }
        }
    }
    for(int i=0;i<NP;i++){
        if(Rx[i]>Lx*0.9){
            Ax[i]=Ax[i]-zeta_M*(Vx[i]-ave_flow)+K_M*genrand_real1();
            Ay[i]=Ay[i]-zeta_M*Vy[i]+K_M*genrand_real1();
        }  
    }





}






void gmap_create4(double L_GX,double L_GY ,int static_pair_length,int move_pair_length,int NP,
int N_GX,int N_GY,
vector<vector<int>> &G_MAP,
vector<vector<int>> &move_pairlist,
vector<int> &static_pairlist,
vector<double> RX,
vector<double> RY,
vector<int> neighbor_row_static,
vector<int> neighbor_col_static,
vector<int> neighbor_row_move,
vector<int> neighbor_col_move,
int neighbor_len_static,int neighbor_len_move){

	int gx_map;
	int gy_map;
    int select_index;



////グリッドとペアリストの初期化//////////
    // for (int i = 0;i < N_GY;i++) {
    //     for (int j = 0;j < N_GX;j++) {
    //         G_MAP[i][j] = -1;
    //     }
    // }
    G_MAP.assign(N_GY,vector<int> (N_GX,-1));

    // for (int i = 0; i < NP; i++) {
    //     for (int j = 0; j < move_pair_length; j++) {
    //         move_pairlist[i][j] = -1;
    //     }
    // }

    move_pairlist.assign(NP,vector<int>(move_pair_length,-1));
    // for (int i=0;i<static_pair_length;i++){
    //     static_pairlist[i]=-1;
    // }
    static_pairlist.assign(static_pair_length,-1);

/////////////////////////////////////////////////
    for (int i=0;i<NP;i++){
        gx_map=int(RX[i]/L_GX);
        gy_map=int(RY[i]/L_GY);

        G_MAP.at(gy_map).at(gx_map)=i;

    }

    for (int i=0;i<N_GY;i++){
        for(int j =0;j<N_GX;j++){
            select_index=G_MAP.at(i).at(j);
            if(select_index!=-1){
                if(select_index==0){
                    make_pairlist(move_pair_length,N_GX,N_GY,j,i,neighbor_len_static,neighbor_row_static,
                                    neighbor_col_static,G_MAP,move_pairlist,static_pairlist);
                }
                else{

                    make_pairlist(move_pair_length,N_GX,N_GY,j,i,neighbor_len_move,neighbor_row_move,
                                    neighbor_col_move,G_MAP,move_pairlist,static_pairlist);
                }
            }
        }
    }

}

void make_pairlist(int move_pair_length,int N_GX,int N_GY,int select_gx,int select_gy,int NEIGHBOR_LEN,
vector<int>NEIGHBOR_ROW,
vector<int> NEIGHBOR_COL,
vector<vector<int>> G_MAP,
vector<vector<int>> &move_PAIRLIST,
vector<int> &static_PAIRLIST){

    int partcle_counter=0;
    int search_gx;
    int search_gy;
    int search_index;
    int select_index=G_MAP.at(select_gy).at(select_gx);





    for(int k=0;k<NEIGHBOR_LEN;k++){
        search_gx=select_gx+NEIGHBOR_COL[k];
        search_gy=select_gy+NEIGHBOR_ROW[k];
        //  グリッドの周期委境界
        if(search_gx>=N_GX){
            search_gx-=N_GX;
        }
        else if (search_gx<0){
            search_gx+=N_GX;
        }
        if (search_gy>=N_GY){
            search_gy-=N_GY;
        }
        else if(search_gy<0){
            search_gy+=N_GY;
        }

        search_index=G_MAP.at(search_gy).at(search_gx);
        if(search_index==0){
            continue;
        }
        if(search_index!=-1){
            if(select_index==0){
                static_PAIRLIST[partcle_counter]=G_MAP[search_gy][search_gx];
            }
            else{
                move_PAIRLIST[select_index][partcle_counter]=search_index;
            }
            partcle_counter+=1;
        }
    }
    move_PAIRLIST[select_index][move_pair_length-1]=partcle_counter;
}



void initialize_ra(double ave_flow,double temp,double m,int remove_dx_num,double dx,double dy,int Nx,int Ny ,double lx,double ly,
vector<double> &rx,vector<double> &ry,vector<double> &ax,vector<double> &ay,
vector<double> &vx,vector<double> &vy){
    double vx_sum=0.0;
    double vy_sum=0.0;
    double kin0=0.0;
    int ave_vx,ave_vy;

    rx.push_back(lx/4.0);
    ry.push_back(ly/2.0);
    vx.push_back(0.0);
    vy.push_back(0.0);
    ax.push_back(0.0);
    ay.push_back(0.0);

    for(int i=0;i<Nx;i++){
        for (int j=0;j<Ny;j++){
            if((j>=int(Ny/2)-remove_dx_num)&&(j<=int(Ny/2)+remove_dx_num)&&(i>=int(Nx/4)-remove_dx_num)&&(i<=int(Nx/4)+remove_dx_num)){
                continue;
            }
            rx.push_back(i*dx+dx/2);
            ry.push_back(j*dy+dy/2);
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

    ave_vx=vx_sum/(vx.size()-1);
    ave_vy=vy_sum/(vy.size()-1);


// 静止粒子の時は、計算しないからi=1から
    for (int i=1;i<vx.size();i++){
        vx[i]=vx[i]-ave_vx;
        vx[i]=vx[i]+ave_flow;
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


void makeNeighborList(int aroundGridNum, vector<int> &neighbor_row, vector<int> &neighbor_col) {
    //xの近接リスト作成
    int neighborRowValue = -aroundGridNum;
    int countIndexRow = 0;
    for (int i = 0; i < aroundGridNum + 1; i++) {
        for (int j = 0; j < aroundGridNum; j ++) {
            neighbor_row[countIndexRow] = neighborRowValue;
            countIndexRow++;
        }
        neighborRowValue++;
    }
    for (int i = 0; i < aroundGridNum; i++) {
        for (int j = 0; j < aroundGridNum + 1; j ++) {
            neighbor_row[countIndexRow] = neighborRowValue;
            countIndexRow++;
        }
        neighborRowValue++;
    }
    //yの近接リスト作成
    int neighborColValue = -aroundGridNum;
    int countIndexCol = 0;
    for (int i = 0; i < aroundGridNum + 1; i++) {
        for (int j = 0; j < aroundGridNum; j ++) {
            neighbor_col[countIndexCol] = neighborColValue;
            countIndexCol++;
            neighborColValue++;
        }
        neighborColValue = -aroundGridNum;
    }
    for (int i = 0; i < aroundGridNum; i++) {
        for (int j = 0; j < aroundGridNum + 1; j ++) {
            neighbor_col[countIndexCol] = neighborColValue;
            countIndexCol++;
            neighborColValue++;
        }
        neighborColValue = -aroundGridNum;
    }
}


