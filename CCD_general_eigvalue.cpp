#include <iostream>
#include "stdio.h"
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iomanip>
#include <time.h>
#include <fstream>
#include "CCD_all.h"
using namespace std;


// solve the generalized non-Hermitian eigenvalue problem(subspace-projected CCSD Hamiltonian is non-Hermitian)
//class general_eigvalue
//{
//    public:
//           general_eigvalue();
//           ~general_eigvalue();
//           void input(double rho, LECs Minnesota_LECs, int subspace_dimension_in);
//           int hhpp_dimension, pp_no, hh_no, subspace_dimension;//subspace_dimension is the number of the new basis 
//    private:
//           void setup_CCD_configuration_space(double rho, LECs Minnesota_LECs);
//           void allocate_all_matrix(int subspace_dimension);
//           void get_N_matrix();
//           void get_H_matrix();
//           void read_all_matrix();
//           
//           tbwf *hh_config, *pp_config;
//           CHANNEL *hhpp_channel_L, *hhpp_channel_R;
//           double **temp_pointer1;//_V, **temp_pointer1_t, **temp_pointer1_H;
//           double *temp_pointer2;//_V, *temp_pointer2_t, *temp_pointer2_H;
//           MATRIX_CHANNEL **V_ijab, **t_ijab, **H_bar_ijab, **l_abij, **x_ijab;
//           double *Ec;
//           double **N_matrix, **H_matrix;
//
//};

general_eigvalue :: general_eigvalue()
{
           hhpp_dimension = 0; 
           pp_no = 0; 
           hh_no = 0;
           Ek    = 0;
}


general_eigvalue :: ~general_eigvalue()
{
}


double general_eigvalue :: inoutput(double rho, LECs Minnesota_LECs, int subspace_dimension_in)
{
//    Minnesota_LECs.Vr = 200;
//    Minnesota_LECs.Vs = -91.85;
//    rho = 0.2;


    subspace_dimension = subspace_dimension_in;
    setup_CCD_configuration_space(rho);
    allocate_all_matrix(subspace_dimension);
    read_all_matrix();


    calculate_H_bar_ij_ab(Minnesota_LECs);
    vacuum_expectation_value_H_bar(Minnesota_LECs);

    get_N_matrix();
    get_H_matrix();

    cout <<"\n N_matrix"<<endl;
    for (int bar = 0; bar < subspace_dimension ; bar++)
    {
        for (int ket = 0; ket < subspace_dimension ; ket++)
        {
            cout << N_matrix[bar][ket]<< " ";
        }
        cout << endl;
    }

    cout <<"\n K_matrix"<<endl;
    for (int bar = 0; bar < subspace_dimension ; bar++)
    {
        for (int ket = 0; ket < subspace_dimension ; ket++)
        {
            cout << K_matrix[bar][ket]<< " ";
        }
        cout << endl;
    }

    cout <<"\n H_matrix"<<endl;
    for (int bar = 0; bar < subspace_dimension ; bar++)
    {
        for (int ket = 0; ket < subspace_dimension ; ket++)
        {
            cout << H_matrix[bar][ket]<< " ";
        }
        cout << endl;
    }

    write_to_file();
    return(0);
}

void general_eigvalue :: write_to_file()
{
    ofstream outfile;
    string file_path;
    file_path = "N_matrix.txt";
    outfile.open(file_path,ios::out);
    if (outfile.is_open())
    {
        outfile << subspace_dimension<<endl;
        for (int bar = 0; bar < subspace_dimension ; bar++)
        {
            for (int ket = 0; ket < subspace_dimension ; ket++)
            {
                outfile << setprecision(10) << N_matrix[bar][ket]<< " ";
            }
            outfile << endl;
        }
    }
    outfile.close();

    file_path = "K_matrix.txt";
    outfile.open(file_path,ios::out);
    if (outfile.is_open())
    {
        outfile << subspace_dimension<<endl;
        for (int bar = 0; bar < subspace_dimension ; bar++)
        {
            for (int ket = 0; ket < subspace_dimension ; ket++)
            {
                outfile << setprecision(10) << K_matrix[bar][ket]<< " ";
            }
            outfile << endl;
        }
    }
    outfile.close();

    file_path = "H_matrix_LEC_1.txt";
    outfile.open(file_path,ios::out);
    if (outfile.is_open())
    {
        outfile << subspace_dimension<<endl;
        for (int bar = 0; bar < subspace_dimension ; bar++)
        {
            for (int ket = 0; ket < subspace_dimension ; ket++)
            {
                outfile << setprecision(10) << H_matrix[bar][ket]<< " ";
            }
            outfile << endl;
        }
    }
    outfile.close();
}


// the overlap matrix elements <psi'|psi> (basis construct with eigenvectors gained by different LECs)
void general_eigvalue :: get_N_matrix()
{
    double N0, N3;
    N_matrix = (double **)malloc(subspace_dimension*sizeof(double*));
    for (int i = 0; i< subspace_dimension ; i++)
    {
        N_matrix[i] = (double *)malloc(subspace_dimension*sizeof(double*));
    }
    for (int bar = 0; bar < subspace_dimension ; bar++)
    {
        for (int ket = 0; ket < subspace_dimension ; ket++)
        {
            N0 = 1;
            N3 = 0;
            for (int k = 0; k < hhpp_dimension; k ++)
            {
                for (int ij = 0; ij < hhpp_channel_L[k].times; ij++)
                {
                    for (int ab = 0; ab < hhpp_channel_R[k].times; ab++)
                    {
                        N3 += 0.25 * l_abij[bar][k].matrix[ij][ab] * (-t_ijab[bar][k].matrix[ij][ab] + t_ijab[ket][k].matrix[ij][ab]);  // x_ijab[ij][ab];
                    }
                }
            }
            N_matrix[bar][ket] = N0 + N3 ;
        }
    }
}


// the H matrix in new basis <psi'|H|psi>(basis construct with eigenvectors gained by different LECs)
void general_eigvalue :: get_H_matrix()
{
    double H0, H3, K0, K3;
    H_matrix = (double **)malloc(subspace_dimension*sizeof(double*));
    K_matrix = (double **)malloc(subspace_dimension*sizeof(double*));
    for (int i = 0; i< subspace_dimension ; i++)
    {
        H_matrix[i] = (double *)malloc(subspace_dimension*sizeof(double*));
        K_matrix[i] = (double *)malloc(subspace_dimension*sizeof(double*));
    }

    for (int bar = 0; bar < subspace_dimension ; bar++)
    {
        for (int ket = 0; ket < subspace_dimension ; ket++)
        {
            //H0 = 0;
            H0 = vacuum_H_bar[ket] * N_matrix[bar][ket];
            K0 = Ek * N_matrix[bar][ket]; 
            H3 = 0;
            K3 = 0;
            for (int k = 0; k < hhpp_dimension; k ++)
            {
                for (int ij = 0; ij < hhpp_channel_L[k].times; ij++)
                {
                    for (int ab = 0; ab < hhpp_channel_R[k].times; ab++)
                    {
                        H3 += 0.25 * l_abij[bar][k].matrix[ij][ab] * H_bar_ijab[ket][k].matrix[ij][ab]; 
                        K3 += 0.25 * l_abij[bar][k].matrix[ij][ab] * kinetic_bar[ket][k].matrix[ij][ab]; 
                    }
                }
            }
            cout <<" H0="<<H0<<"  H3="<<H3<<endl;
            cout <<" K0="<<K0<<"  K3="<<K3<<endl;
            H_matrix[bar][ket] = H0 + H3 ;
            K_matrix[bar][ket] = K0 + K3 ;
        }
    }

}

void general_eigvalue :: vacuum_expectation_value_H_bar(LECs Minnesota_LECs)
{
    double H0, H3;
    H0 = 0;

    vacuum_H_bar = (double *)malloc(subspace_dimension*sizeof(double));
 
    /////////////////////////////////////////////////////////////////////////
    //*********************************************************************//
    //***************************** BENCHMARK *****************************//
    //************** For Kinetic energy and HATREE-FOCK energy ************//
    //********************************************************************///
    /////////////////////////////////////////////////////////////////////////
    double kk = 0.;
    //      cout<<"L = "<<L<<endl;
    for(int i = 0; i<spstates_no; i++) 
    {
            if(config[i].occupied == 1)
            {
                    kk = (double)(pow(config[i].kx,2) + pow(config[i].ky,2) + pow(config[i].kz,2)) * pow((2*pi)/L,2);
                    Ek += pow(h_bar,2) * kk / (2. * mass);
            }
    }
    double hf_energy=0;
    //hf_energy = Ek;
    for(int ij = 0; ij < hh_no; ij++)
    {
            hf_energy += 0.5 * V_ks_AS(L, hh_config[ij].q_kx, hh_config[ij].q_ky, hh_config[ij].q_kz, hh_config[ij].sz1, hh_config[ij].sz2, \
            hh_config[ij].q_kx, hh_config[ij].q_ky, hh_config[ij].q_kz, hh_config[ij].sz1, hh_config[ij].sz2, Minnesota_LECs);
    }
    cout <<"\nEhf = "<<setprecision(7)<<(hf_energy+Ek);    

    for(int ket = 0; ket < subspace_dimension; ket ++)
    {
        H0 = hf_energy;
        H3 = 0;
        for(int i = 0; i < hhpp_dimension; i++)
        {
            for (int ij = 0; ij < hhpp_channel_L[i].times; ij++)
            {
                for (int ab = 0; ab < hhpp_channel_R[i].times; ab++)
                {
                   H3 += 0.25 * V_ijab[i].matrix[ij][ab]*t_ijab[ket][i].matrix[ij][ab];
                }
            }
        }
        vacuum_H_bar[ket] = H0 + H3;
        //vacuum_H_bar[ket] = H3;
    }
}


void general_eigvalue :: setup_CCD_configuration_space(double rho)
{
        L = pow((double)magic_no/rho, 1./3.);
        cout <<"L="<<L<<endl;
        spstates_no = 2 * pow((2*Nmax+1),3);
        config = (spstate *)malloc(spstates_no * sizeof(spstate));


        ////////////////////////////////////////////////////////////////////////
        //********************************************************************//
        //********************************************************************//
        //********************************************************************//
        //*******************  GENERATE ALL THE S.P. STATE  ******************//
        //********************************************************************//
        //********************************************************************//
        //********************************************************************//
        ////////////////////////////////////////////////////////////////////////


        /////////////////////////////////////////////////////////////////////////////////
        //*****************************************************************************//
        //*** Actually it is a MAP more than fullfill the s.p. configuration space. ***//
        //*****************************************************************************//
        /////////////////////////////////////////////////////////////////////////////////

        /*      for(int i = 0; i<spstates_no; i++)
        {
                config[i].sz = (i%2) - 0.5 ;
                config[i].kz = (i/2) % (2*(2*Nmax+1)) - Nmax;
                config[i].ky = ((i/2) / (2*(2*Nmax+1))) % (2*(2*Nmax+1)) - Nmax;
                config[i].kx = (((i/2) / (2*(2*Nmax+1))) / (2*(2*Nmax+1))) % (2*(2*Nmax+1)) - Nmax;
        }*/


        ////////////////////////////////////////////////////////////////////////
        //********************************************************************//
        //*** Alternative choice to generate single-particle configuration ***//
        //*******  More straightforward but might be less efficient.  ********//
        //********************************************************************//
        ////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////
        //********************************************************************//
        //********* Actually we could make the code more efficient. **********//
        //*********** by store the hole and particle seperately. *************//
        //********************************************************************//
        ////////////////////////////////////////////////////////////////////////

        int loop1 = 0;
        double mod;
        do
        {
                for(int kx = -Nmax; kx <= Nmax; kx++)
                {
                        for(int ky = -Nmax; ky <= Nmax; ky++)
                        {
                                for (int kz = -Nmax; kz <= Nmax; kz++)
                                {
                                        for(int sz = -1 ; sz <=1; sz+=2)
                                        {
                                                config[loop1].kx = kx;
                                                config[loop1].ky = ky;
                                                config[loop1].kz = kz;
                                                config[loop1].sz = sz;
                                                mod = pow(pow(config[loop1].kx,2) + pow(config[loop1].ky,2)+pow(config[loop1].kz,2),0.5);
                                                if(mod == 0 || mod == 1) config[loop1].occupied = 1;
                                                else config[loop1].occupied = 0;
                                                loop1++;
                                        }
                                }
                        }
                }
        }while(loop1<spstates_no);

        /////////////////////////////////////////////////////////////////////////
        //*********************************************************************//
        //*** GENERATE THE CONFIGURATIONS IN RELATIVE MOMENTUM PRESENTATION ***//
        //********************************************************************///
        /////////////////////////////////////////////////////////////////////////

//      double pp_no = (spstates_no - magic_no) * (spstates_no - magic_no - 1);//combinator((spstates_no - magic_no) , 2); // number of pp-configurations in relative presentation
        double hh_no_old = 2 * combinator(magic_no, 2);
//      double ph_no = 2 * (spstates_no - magic_no) * magic_no;

//      tbwf *ph_config = (tbwf *)malloc(ph_no * sizeof(tbwf));



        int loop2 = 0;
        int loop3 = 0;
        int loop4 = 0;
        ////////  pp-configuration-number /////
        //double pp_no = 0;
        for(int k1 = 0; k1<spstates_no; k1++)
        {
                for(int k2 = 0; k2<spstates_no; k2++)
                {
                        if((k1 != k2) && (config[k1].occupied == 0) && (config[k2].occupied == 0) &&(config[k1].sz*config[k2].sz == -1))
                        {
                                pp_no ++;
                        }
                }
        }
//      cout <<  "pp_no_new = "<<pp_no<<endl;
        ////////  pp-configuration-setup /////
        pp_config = (tbwf *)malloc(pp_no * sizeof(tbwf)); //configurations in relative presentation

                for(int k1 = 0; k1<spstates_no; k1++)
                {
                        for(int k2 = 0; k2<spstates_no; k2++)
                        {
                                if((k1 != k2) && (config[k1].occupied == 0) && (config[k2].occupied == 0) && (config[k1].sz*config[k2].sz == -1)) // S12 =0 for this V
                                {
                                //      cout<<"!!"<<endl;
                                        pp_config[loop2].Q_kx = config[k1].kx + config[k2].kx;
                                        pp_config[loop2].Q_ky = config[k1].ky + config[k2].ky;
                                        pp_config[loop2].Q_kz = config[k1].kz + config[k2].kz;

                                        pp_config[loop2].q_kx = (double)(config[k1].kx - config[k2].kx) / 2.;
                                        pp_config[loop2].q_ky = (double)(config[k1].ky - config[k2].ky) / 2.;
                                        pp_config[loop2].q_kz = (double)(config[k1].kz - config[k2].kz) / 2.;

                                        pp_config[loop2].sz1 = config[k1].sz;
                                        pp_config[loop2].sz2 = config[k2].sz;

                                        pp_config[loop2].flag = 1;

                                        loop2++;
                                }
                        }
                }
//              cout <<  "loop2 = "<<loop2<<endl;



        ////////  hh-configuration-number /////
        //double hh_no = 0;
        for(int k1 = 0; k1<spstates_no; k1++)
        {
                for(int k2 = 0; k2<spstates_no; k2++)
                {
                        if((k1 != k2) && (config[k1].occupied == 1) && (config[k2].occupied == 1) &&  (config[k1].sz*config[k2].sz == -1))
                        {
                                hh_no ++;
                        }
                }
        }

        ////////  hh-configuration-setup /////
        hh_config = (tbwf *)malloc(hh_no * sizeof(tbwf));
                for(int k1 = 0; k1<spstates_no; k1++)
                {
                        for(int k2 = 0; k2<spstates_no; k2++)
                        {
                                if((k1 != k2) && (config[k1].occupied == 1) && (config[k2].occupied == 1) &&  (config[k1].sz*config[k2].sz == -1))
                                {
                                        hh_config[loop3].Q_kx = config[k1].kx + config[k2].kx;
                                        hh_config[loop3].Q_ky = config[k1].ky + config[k2].ky;
                                        hh_config[loop3].Q_kz = config[k1].kz + config[k2].kz;

                                        hh_config[loop3].q_kx = (double)(config[k1].kx - config[k2].kx) / 2.;
                                        hh_config[loop3].q_ky = (double)(config[k1].ky - config[k2].ky) / 2.;
                                        hh_config[loop3].q_kz = (double)(config[k1].kz - config[k2].kz) / 2.;

                                        hh_config[loop3].sz1 = config[k1].sz;
                                        hh_config[loop3].sz2 = config[k2].sz;
                                        hh_config[loop3].flag = 1;
                                        loop3++;
                                }
                        }
                }
        ////////////////////////////////////////////////////////////////////////
        //********************************************************************//
        //********************************************************************//
        //********************************************************************//
        //***************  GENERATE THE CHANNEL CONFIGURATIONS  **************//
        //********************************************************************//
        //********************************************************************//
        //********************************************************************//
        ////////////////////////////////////////////////////////////////////////



        ////////////////////////////////////////////////////////////////////////
        //********************************************************************//
        //********************************************************************//
        //********************************************************************//
        //*************************** HH-HH - CHANNEL  ***********************//
        //********************************************************************//
        //********************************************************************//
        //********************************************************************//
        ////////////////////////////////////////////////////////////////////////

        int k,loop5;
        int hh_temp_x, hh_temp_y, hh_temp_z;
        hh_dimension = 0;

        for(int i = 0; i < hh_no; i++)
        {
                if(hh_config[i].flag == 1)
                {
                        hh_temp_x = hh_config[i].Q_kx;
                        hh_temp_y = hh_config[i].Q_ky;
                        hh_temp_z = hh_config[i].Q_kz;
                        for(int j = i; j<hh_no; j++)
                        {
                                if(hh_config[j].Q_kx == hh_temp_x && hh_config[j].Q_ky == hh_temp_y && hh_config[j].Q_kz == hh_temp_z)
                                {
                                        hh_config[j].flag = 0;
                                }
                        }
                        hh_dimension++; // 这记录了hh里有多少个不同的总动量
                }
        }

        for(int i = 0; i<hh_no; i++)
        {
                hh_config[i].flag = 1;
        }

        hh_channel = (CHANNEL *)malloc(hh_dimension * sizeof(CHANNEL)); //先声明！！！！！！！
        int hh_temp_dimension = 0;
        for(int i = 0; i< hh_no; i++)
        {
                if(hh_config[i].flag == 1)
                {
                        hh_channel[hh_temp_dimension].Qx = hh_config[i].Q_kx;
                        hh_channel[hh_temp_dimension].Qy = hh_config[i].Q_ky;
                        hh_channel[hh_temp_dimension].Qz = hh_config[i].Q_kz;
                        hh_channel[hh_temp_dimension].times = 0;

                        for(int j = i; j<hh_no; j++)
                        {
                                if(hh_config[j].Q_kx == hh_channel[hh_temp_dimension].Qx && hh_config[j].Q_ky == hh_channel[hh_temp_dimension].Qy && hh_config[j].Q_kz == hh_channel[hh_temp_dimension].Qz)
                                {
                                        hh_config[j].flag = 0;
                                        hh_channel[hh_temp_dimension].times++;
                                }
                        }
                        //为每一个hh_channel[hh_temp_dimension]动态分配
                        tbwf *temp_star_hh;
                        temp_star_hh = (tbwf *)malloc(hh_channel[hh_temp_dimension].times * sizeof(tbwf));
                        hh_channel[hh_temp_dimension].wf = temp_star_hh;
                        hh_temp_dimension++;
                }
        }

        for(int i = 0; i<hh_no; i++) //将hh_config[i]的flag全都赋值回1
        {
                hh_config[i].flag = 1;
        }

        loop5 = 0; 
        k = 0; 
        for(int i = 0; i< hh_no; i++) 
        {
                if(hh_config[i].flag == 1)
                {
                        k = 0; 
                        for(int j = i; j<hh_no; j++) 
                        {
                                if(hh_config[i].Q_kx == hh_config[j].Q_kx && hh_config[i].Q_ky == hh_config[j].Q_ky && hh_config[i].Q_kz == hh_config[j].Q_kz)
                                {
                                        hh_channel[loop5].wf[k].q_kx = hh_config[j].q_kx;
                                        hh_channel[loop5].wf[k].q_ky = hh_config[j].q_ky;
                                        hh_channel[loop5].wf[k].q_kz = hh_config[j].q_kz;
                                        hh_channel[loop5].wf[k].sz1 = hh_config[j].sz1;
                                        hh_channel[loop5].wf[k].sz2 = hh_config[j].sz2;
                                        hh_config[j].flag = 0;     
                                        k++;
                                }
                        }
                        loop5++;
                }
        }



        ////////////////////////////////////////////////////////////////////////
        //********************************************************************//
        //********************************************************************//
        //********************************************************************//
        //*************************  HH-PP - CHANNEL  ************************//
        //********************************************************************//
        //********************************************************************//
        //********************************************************************//
        ////////////////////////////////////////////////////////////////////////
        int hhpp_temp_x, hhpp_temp_y, hhpp_temp_z;
        hhpp_dimension = 0;
        for(int i = 0; i<hh_no; i++) //将hh_config[i]的flag全都赋值回1
        {
                hh_config[i].flag = 1;
        }
        for(int i = 0; i < hh_no; i++)
        {
                if(hh_config[i].flag == 1)
                {
                        hhpp_temp_x = hh_config[i].Q_kx;
                        hhpp_temp_y = hh_config[i].Q_ky;
                        hhpp_temp_z = hh_config[i].Q_kz;
                        for(int j = 0; j<hh_no; j++)
                        {
                                if(hh_config[j].Q_kx == hhpp_temp_x && hh_config[j].Q_ky == hhpp_temp_y && hh_config[j].Q_kz == hhpp_temp_z)
                                {
                                        hh_config[j].flag = 0;
                                }
                        }
                        for(int j = 0; j<pp_no; j++)
                        {
                                if(pp_config[j].Q_kx == hhpp_temp_x && pp_config[j].Q_ky == hhpp_temp_y && pp_config[j].Q_kz == hhpp_temp_z)
                                {
                                        //hh_config[j].flag = 0;
                                        hhpp_dimension++;
                                        break;
                                }
                        }
                        //hhpp_dimension++; // 这记录了hh里有多少个不同的总动量
                }
        }
        //CHANNEL hhpp_channel_L[hhpp_dimension]; //先声明！！！！！！！
        //CHANNEL hhpp_channel_R[hhpp_dimension]; //先声明！！！！！！！
        hhpp_channel_L = (CHANNEL *)malloc(sizeof(CHANNEL) * hhpp_dimension);
        hhpp_channel_R = (CHANNEL *)malloc(sizeof(CHANNEL) * hhpp_dimension);

        for(int i = 0; i<hh_no; i++) //将hh_config[i]的flag全都赋值回1
        {
                hh_config[i].flag = 1;
        }
        for(int i = 0; i<pp_no; i++) //将pp_config[i]的flag全都赋值回1
        {
                pp_config[i].flag = 1;
        }

        int hhpp_temp_dimension = 0;
        int Qx,Qy,Qz,times_L,times_R;
        for(int i = 0; i< hh_no; i++)
        {
                if(hh_config[i].flag == 1)
                {
                        Qx = hh_config[i].Q_kx; //这里先不记录hhpp_channel_L[hhpp_temp_dimension].Q，因为这里事实上并未确定在pp里一定能找到匹配的hh
                        Qy = hh_config[i].Q_ky;
                        Qz = hh_config[i].Q_kz;
                        times_L = 0;
                        times_R = 0;
                        for(int j = i; j<hh_no; j++)
                        {
                                if(hh_config[j].Q_kx == Qx && hh_config[j].Q_ky == Qy && hh_config[j].Q_kz == Qz)
                                {
                                        hh_config[j].flag = 0;
                                        times_L ++;
                                        //hhpp_channel_L[hhpp_temp_dimension].times++;
                                }
                        }
                    for(int j = 0; j<pp_no; j++)
                        {
                                if(pp_config[j].Q_kx == Qx && pp_config[j].Q_ky == Qy && pp_config[j].Q_kz == Qz)
                                {
                                        pp_config[j].flag = 0;
                                        times_R ++;
                                        //hhpp_channel_R[hhpp_temp_dimension].times++;
                                }
                        }
            if(times_R != 0)
            {
                hhpp_channel_L[hhpp_temp_dimension].Qx = Qx;
                hhpp_channel_L[hhpp_temp_dimension].Qy = Qy;
                hhpp_channel_L[hhpp_temp_dimension].Qz = Qz;
                hhpp_channel_L[hhpp_temp_dimension].times = times_L;

                hhpp_channel_R[hhpp_temp_dimension].Qx = Qx;
                hhpp_channel_R[hhpp_temp_dimension].Qy = Qy;
                hhpp_channel_R[hhpp_temp_dimension].Qz = Qz;
                hhpp_channel_R[hhpp_temp_dimension].times = times_R;
            }
            else
            {
                continue; //如果times_R =0,即找不到匹配的pp和hh，那么不进行后面的分配，直接进行下一轮循环
            }
               //      cout<<"here = "<<hhpp_channel_L[hhpp_temp_dimension].times<<endl;

                        //为每一个hhpp_channel[hh_temp_dimension]动态分配
                        tbwf *temp_star_hhpp_L;
                        temp_star_hhpp_L = (tbwf *)malloc(hhpp_channel_L[hhpp_temp_dimension].times * sizeof(tbwf));
                        hhpp_channel_L[hhpp_temp_dimension].wf = temp_star_hhpp_L;

                        tbwf *temp_star_hhpp_R;
                        temp_star_hhpp_R = (tbwf *)malloc(hhpp_channel_R[hhpp_temp_dimension].times * sizeof(tbwf));
                        hhpp_channel_R[hhpp_temp_dimension].wf = temp_star_hhpp_R;
                        hhpp_temp_dimension++;
                }
        }


        for(int i = 0; i<hh_no; i++) //将hh_config[i]的flag全都赋值回1
        {
                hh_config[i].flag = 1;
        }
        for(int i = 0; i<hh_no; i++) //将hh_config[i]的flag全都赋值回1
        {
                pp_config[i].flag = 1;
        }

        //////////////////////////////////////////
        //**************************************//
        ///////////重要判定！！！！！/////////////
        int xx;
        if(hh_dimension != hhpp_dimension)
        {
                cout<<"hh_dimension = "<<hh_dimension;
                cout<<"hhpp_dimension = "<<hhpp_dimension;
                cout<<"\nERROR!"<<endl;
                cin>>xx;
        }
    for(int i = 0; i<hhpp_dimension; i++) //循环所有hhpp_channel 并且对hhpp_channel[].L_wf R_wf 赋值
    {
        k = 0;
        //cout<<"channel = "<<i<<endl;
        //cout <<" qx="<<hhpp_channel[i].Qx<<" qy="<<hhpp_channel[i].Qy<<" qz="<<hhpp_channel[i].Qz<<endl;

        for(int j = 0; j< hh_no; j++)
        {

                //cout <<" hqx="<<hh_config[j].Q_kx<<" hqy="<<hh_config[j].Q_ky<<" hqz="<<hh_config[j].Q_kz<<endl;

                if(hh_config[j].Q_kx == hhpp_channel_L[i].Qx && hh_config[j].Q_ky == hhpp_channel_L[i].Qy && hh_config[j].Q_kz == hhpp_channel_L[i].Qz)
                        {
                                //cout <<" hqx="<<hh_config[j].Q_kx<<" hqy="<<hh_config[j].Q_ky<<" hqz="<<hh_config[j].Q_kz<<endl;
                   // cout <<" hqkx="<<hh_config[j].q_kx<<" hqky="<<hh_config[j].q_ky<<" hqkz="<<hh_config[j].q_kz<<endl;
                                hhpp_channel_L[i].wf[k].q_kx = hh_config[j].q_kx;
                                hhpp_channel_L[i].wf[k].q_ky = hh_config[j].q_ky;
                                hhpp_channel_L[i].wf[k].q_kz = hh_config[j].q_kz;
                                hhpp_channel_L[i].wf[k].sz1  = hh_config[j].sz1;
                                hhpp_channel_L[i].wf[k].sz2  = hh_config[j].sz2;
                k ++ ;
                        }


        }
         if (k != hhpp_channel_L[i].times)
                        {
                                cout<<"wtf erro1#"<<endl;
                        }

        k = 0;
        for(int j = 0; j< pp_no; j++)
        {
                if(pp_config[j].Q_kx == hhpp_channel_R[i].Qx && pp_config[j].Q_ky == hhpp_channel_R[i].Qy && pp_config[j].Q_kz == hhpp_channel_R[i].Qz)
                        {
                                hhpp_channel_R[i].wf[k].q_kx = pp_config[j].q_kx;
                                hhpp_channel_R[i].wf[k].q_ky = pp_config[j].q_ky;
                                hhpp_channel_R[i].wf[k].q_kz = pp_config[j].q_kz;
                                hhpp_channel_R[i].wf[k].sz1  = pp_config[j].sz1;
                                hhpp_channel_R[i].wf[k].sz2  = pp_config[j].sz2;
                                k ++ ;
                        }
        }
        if (k != hhpp_channel_R[i].times)
                {
                        cout<< " k = "<<k<<" hhpp_times="<<hhpp_channel_R[i].times<<endl;
                        cout<<"wtf erro2#"<<endl;
                }
    }



        ////////////////////////////////////////////////////////////////////////
        //********************************************************************//
        //********************************************************************//
        //********************************************************************//
        //*************************  PP-PP- CHANNEL  *************************//
        //********************************************************************//
        //********************************************************************//
        //********************************************************************//
        ////////////////////////////////////////////////////////////////////////
        pp_channel = (CHANNEL *)malloc( hhpp_dimension * sizeof(CHANNEL)); //先声明！！！！！！！

        for(int i = 0; i<pp_no; i++)
        {
                pp_config[i].flag = 1;
        }

        for(int i = 0; i< hhpp_dimension; i++)
        {
                pp_channel[i].Qx = hhpp_channel_L[i].Qx;
                pp_channel[i].Qy = hhpp_channel_L[i].Qy;
                pp_channel[i].Qz = hhpp_channel_L[i].Qz;
                pp_channel[i].times = hhpp_channel_R[i].times;

                //为每一个pp_channel[pp_temp_dimension]动态分配
                tbwf *temp_star;
                temp_star = (tbwf *)malloc(pp_channel[i].times * sizeof(tbwf));
                pp_channel[i].wf = temp_star;
        }

        for(int i = 0; i< hhpp_dimension; i++)
        {
                for(int j = 0; j < pp_channel[i].times; j++)
                {
                        pp_channel[i].wf[j].q_kx = hhpp_channel_R[i].wf[j].q_kx;
                        pp_channel[i].wf[j].q_ky = hhpp_channel_R[i].wf[j].q_ky;
                        pp_channel[i].wf[j].q_kz = hhpp_channel_R[i].wf[j].q_kz;
                        pp_channel[i].wf[j].sz1  = hhpp_channel_R[i].wf[j].sz1;
                        pp_channel[i].wf[j].sz2  = hhpp_channel_R[i].wf[j].sz2;
                }
        }






}


void general_eigvalue :: allocate_all_matrix(int subspace_dimension)
{

        ///////////////////////////////////////////////////////////
        //*******************************************************//
        //******** allocate all the matrix that are needed ******//
        //*******************************************************//
        ///////////////////////////////////////////////////////////

        t_ijab     = (MATRIX_CHANNEL **)malloc(subspace_dimension * sizeof(MATRIX_CHANNEL *));
        x_ijab     = (MATRIX_CHANNEL **)malloc(subspace_dimension * sizeof(MATRIX_CHANNEL *));
        l_abij     = (MATRIX_CHANNEL **)malloc(subspace_dimension * sizeof(MATRIX_CHANNEL *));
        H_bar_ijab = (MATRIX_CHANNEL **)malloc(subspace_dimension * sizeof(MATRIX_CHANNEL *));
        kinetic_bar= (MATRIX_CHANNEL **)malloc(subspace_dimension * sizeof(MATRIX_CHANNEL *));
 
        for(int i = 0; i< subspace_dimension; i++)
        {
            t_ijab[i]     = (MATRIX_CHANNEL *)malloc(hhpp_dimension * sizeof (MATRIX_CHANNEL));
            x_ijab[i]     = (MATRIX_CHANNEL *)malloc(hhpp_dimension * sizeof (MATRIX_CHANNEL));
            l_abij[i]     = (MATRIX_CHANNEL *)malloc(hhpp_dimension * sizeof (MATRIX_CHANNEL));
            H_bar_ijab[i] = (MATRIX_CHANNEL *)malloc(hhpp_dimension * sizeof (MATRIX_CHANNEL));
            kinetic_bar[i]= (MATRIX_CHANNEL *)malloc(hhpp_dimension * sizeof (MATRIX_CHANNEL));
        }

        cout<< "subspace_dimension = "<<subspace_dimension<<endl;
        cout<< "hhpp_dimension = "<<hhpp_dimension<<endl;


        for(int i = 0; i< subspace_dimension; i++)
        {
            for(int j = 0; j< hhpp_dimension; j++)
            {
                temp_pointer1 = (double **)malloc(hhpp_channel_L[j].times * sizeof(double*));
                t_ijab[i][j].matrix = temp_pointer1;

        //      temp_pointer1 = (double **)malloc(hhpp_channel_L[j].times * sizeof(double*));
        //      x_ijab[i][j].matrix = temp_pointer1;

                temp_pointer1 = (double **)malloc(hhpp_channel_L[j].times * sizeof(double*));
                l_abij[i][j].matrix = temp_pointer1;

                temp_pointer1 = (double **)malloc(hhpp_channel_L[j].times * sizeof(double*));
                H_bar_ijab[i][j].matrix = temp_pointer1;

                temp_pointer1 = (double **)malloc(hhpp_channel_L[j].times * sizeof(double*));
                kinetic_bar[i][j].matrix = temp_pointer1;

                for(int ij = 0; ij < hhpp_channel_L[j].times; ij++)
                {
                    temp_pointer2 = (double *)malloc(hhpp_channel_R[j].times * sizeof(double));
                    t_ijab[i][j].matrix[ij] = temp_pointer2;

          //          temp_pointer2 = (double *)malloc(hhpp_channel_R[j].times * sizeof(double));
          //          x_ijab[i][j].matrix[ij] = temp_pointer2;

                    temp_pointer2 = (double *)malloc(hhpp_channel_R[j].times * sizeof(double));
                    l_abij[i][j].matrix[ij] = temp_pointer2;

                    temp_pointer2 = (double *)malloc(hhpp_channel_R[j].times * sizeof(double));
                    H_bar_ijab[i][j].matrix[ij] = temp_pointer2;

                    temp_pointer2 = (double *)malloc(hhpp_channel_R[j].times * sizeof(double));
                    kinetic_bar[i][j].matrix[ij] = temp_pointer2;
               }
            }
        }
       // t_ijab[0][1].matrix[0][0] = 0;

        Ec = (double *)malloc(subspace_dimension *sizeof(double));

        int store_space = 0;
        for (int i = 0; i < subspace_dimension; i++)
        {
            for (int j = 0; j < hhpp_dimension; j++)
            {
                store_space += hhpp_channel_L[j].times * hhpp_channel_R[j].times ;//* sizeof(double);
            }
        }



        //// Assignment for the V_ijkl ////
        //为每一个pp_channel[pp_temp_dimension]动态分配
        V_ijkl = (MATRIX_CHANNEL *)malloc(hh_dimension * sizeof(MATRIX_CHANNEL));
        for(int i = 0; i< hh_dimension; i++)
        {
                temp_pointer1 = (double **)malloc(hh_channel[i].times * sizeof(double*));
                V_ijkl[i].matrix = temp_pointer1;

                for(int ij = 0; ij < hh_channel[i].times; ij++)
                {
                        temp_pointer2 = (double *)malloc(hh_channel[i].times * sizeof(double));
                        V_ijkl[i].matrix[ij] = temp_pointer2;
                }
        }

        //// Assignment for the V_cdab ////
        //为每一个pp_channel[pp_temp_dimension]动态分配
        V_cdab = (MATRIX_CHANNEL *)malloc(hhpp_dimension * sizeof(MATRIX_CHANNEL));
        for(int i = 0; i< hhpp_dimension; i++)
        {
                temp_pointer1 = (double **)malloc(pp_channel[i].times * sizeof(double*));
                V_cdab[i].matrix = temp_pointer1;

                for(int cd = 0; cd < pp_channel[i].times; cd++)
                {
                        temp_pointer2 = (double *)malloc(pp_channel[i].times * sizeof(double));
                        V_cdab[i].matrix[cd] = temp_pointer2;
                }
        }

        //// Assignment for the V_ijab ////
        V_ijab = (MATRIX_CHANNEL *)malloc(hhpp_dimension * sizeof(MATRIX_CHANNEL));
        for(int i = 0; i< hhpp_dimension; i++)
        {
                temp_pointer1 = (double **)malloc(hhpp_channel_L[i].times * sizeof(double*));
                V_ijab[i].matrix = temp_pointer1;

                for(int ij = 0; ij < hhpp_channel_L[i].times; ij++)
                {
                        temp_pointer2 = (double *)malloc(hhpp_channel_R[i].times * sizeof(double));
                        V_ijab[i].matrix[ij] = temp_pointer2;
                }
        }

//        cout  << "size_of_double :"<< sizeof(double) << endl;
        cout  << "subspace_matrix store : "<< store_space/subspace_dimension << endl;
}



//read in all the t, H_bar, matrix for further calculation
void general_eigvalue :: read_all_matrix()
{
    ifstream infile_1;
    string file_path, ss_temp;
    int count;
    cout << "subspace_dimension= "<<subspace_dimension<<endl;
    for (int i= 0 ; i < subspace_dimension; i++)
    {
            file_path = to_string(i) +".txt";
            infile_1.open(file_path,ios::in);
            if(!infile_1)
            {
                cerr<<"open error! #1"<<endl;
                cerr<<file_path<<endl;
                exit(1);
            }

            getline(infile_1,ss_temp);
            infile_1 >> Ec[i];
            infile_1 >> ss_temp;

            for(int j = 0; j< hhpp_dimension; j++)
            {
                for(int ij = 0; ij < hhpp_channel_L[j].times; ij++)
                {
                    for(int ab = 0; ab < hhpp_channel_R[j].times; ab++)
                    {
                        infile_1 >> t_ijab[i][j].matrix[ij][ab];
                        l_abij[i][j].matrix[ij][ab]= t_ijab[i][j].matrix[ij][ab];
                    }
                }
            }
            infile_1.close();
            cout << Ec[i]<<endl;
    }
}


void general_eigvalue :: calculate_H_bar_ij_ab(LECs Minnesota_LECs)  //introduce the similarity-transformed target Hamiltonian (with target LECs) with respect to right |psi>.
{
    int iter = 0;
    int sp, kp_x, kp_y, kp_z;
    double t_no_this, t_last;
    double alpha = 0.9; //*** Here we could try different alpha!!! ****//
    double temp1, temp2, delta;
    double f_ii, f_jj, f_aa, f_bb;  
    double f_ii_kinetic, f_jj_kinetic, f_aa_kinetic, f_bb_kinetic;  

    double h0, h1, h2, h3, h4, k34;
    


//V_ij_ab and all the one body f term are calclulate under new LECs (c_i = 1, c_other = 0)
//This mean we project each part of the hamiltonian in the given right |psi> (with the correspoding t matrix)
    /////////////////////////////////////////////////////////////////////////
    //*********************************************************************//
    //******************** FULFILL THE V_ij_ab MATRIX *********************//
    //********************************************************************///
    /////////////////////////////////////////////////////////////////////////
    for(int i = 0; i< hhpp_dimension; i++)
    {
        // cout<<"Q = "<<pp_channel[i].Qx<<" "<<pp_channel[i].Qy<<" "<<pp_channel[i].Qz;
        // cout<<" times = "<<pp_channel[i].times<<endl;
        V_ijab[i].Qx = hhpp_channel_L[i].Qx;
        V_ijab[i].Qy = hhpp_channel_L[i].Qy;
        V_ijab[i].Qz = hhpp_channel_L[i].Qz;
        for(int ij = 0; ij < hhpp_channel_L[i].times; ij++)
        {
            for(int ab = 0; ab < hhpp_channel_R[i].times; ab++)
            {
                V_ijab[i].matrix[ij][ab]= V_ks_AS(L, hhpp_channel_L[i].wf[ij].q_kx, hhpp_channel_L[i].wf[ij].q_ky, hhpp_channel_L[i].wf[ij].q_kz, hhpp_channel_L[i].wf[ij].sz1, hhpp_channel_L[i].wf[ij].sz2,\
                hhpp_channel_R[i].wf[ab].q_kx, hhpp_channel_R[i].wf[ab].q_ky, hhpp_channel_R[i].wf[ab].q_kz, hhpp_channel_R[i].wf[ab].sz1, hhpp_channel_R[i].wf[ab].sz2, Minnesota_LECs);
    //         cout<<V_ijab[i].matrix[ij][ab]<<" ";
            }
        }
    }
   
        ///////////////////////////////////////////////////////////
        //*******************************************************//
        //************** FULFILL THE V_ij_kl MATRIX *************//
        //*******************************************************//
        //*******************************************************//
        for(int i = 0; i< hh_dimension; i++)
        {
//              cout<<"Q = "<<pp_channel[i].Qx<<" "<<pp_channel[i].Qy<<" "<<pp_channel[i].Qz;
//              cout<<" times = "<<pp_channel[i].times<<endl;
                V_ijkl[i].Qx = hh_channel[i].Qx;
                V_ijkl[i].Qy = hh_channel[i].Qy;
                V_ijkl[i].Qz = hh_channel[i].Qz;

                for(int ab = 0; ab < hh_channel[i].times; ab++)
                {
                    for(int cd = 0; cd < hh_channel[i].times; cd++)
                    {
                        V_ijkl[i].matrix[cd][ab]= V_ks_AS(L, hh_channel[i].wf[cd].q_kx, hh_channel[i].wf[cd].q_ky, hh_channel[i].wf[cd].q_kz, hh_channel[i].wf[cd].sz1, hh_channel[i].wf[cd].sz2,\
                        hh_channel[i].wf[ab].q_kx,hh_channel[i].wf[ab].q_ky, hh_channel[i].wf[ab].q_kz, hh_channel[i].wf[ab].sz1, hh_channel[i].wf[ab].sz2, Minnesota_LECs);
                       //if (i ==1) {cout<<"ab="<< ab<<" cd="<<cd<<" V_ijkl= "<<V_ijkl[i].matrix[cd][ab]<<" "<<"L= "<<L<<" q_kx="<<hh_channel[i].wf[cd].q_kx<<" q_ky="<<hh_channel[i].wf[cd].q_ky<<" q_kz="<< hh_channel[i].wf[cd].q_kz <<endl;}
                    }
                }
        }


        ///////////////////////////////////////////////////////////
        //*******************************************************//
        //************** FULFILL THE V_cd_ab MATRIX *************//
        //*******************************************************//
        //*******************************************************//
        for(int i = 0; i< hhpp_dimension; i++)
        {
//              cout<<"Q = "<<pp_channel[i].Qx<<" "<<pp_channel[i].Qy<<" "<<pp_channel[i].Qz;
//              cout<<" times = "<<pp_channel[i].times<<endl;
                V_cdab[i].Qx = pp_channel[i].Qx;
                V_cdab[i].Qy = pp_channel[i].Qy;
                V_cdab[i].Qz = pp_channel[i].Qz;

                for(int ab = 0; ab < pp_channel[i].times; ab++)
                {
                    for(int cd = 0; cd < pp_channel[i].times; cd++)
                    {
                        V_cdab[i].matrix[cd][ab]= V_ks_AS(L, pp_channel[i].wf[cd].q_kx, pp_channel[i].wf[cd].q_ky, pp_channel[i].wf[cd].q_kz, pp_channel[i].wf[cd].sz1, pp_channel[i].wf[cd].sz2,\
                        pp_channel[i].wf[ab].q_kx,pp_channel[i].wf[ab].q_ky, pp_channel[i].wf[ab].q_kz, pp_channel[i].wf[ab].sz1, pp_channel[i].wf[ab].sz2, Minnesota_LECs);
                        //cout<<V_cdab[i].matrix[cd][ab]<<" ";
                    }
                }
        }



    /////////////////////////////////////////////////////////////////////////
    //*********************************************************************//
    //************* Calculate the H_bar_ij_ab matrix  *********************//
    //********************************************************************///
    /////////////////////////////////////////////////////////////////////////
    for (int ket = 0; ket < subspace_dimension; ket ++)
    {
        for(int i = 0; i< hhpp_dimension; i++)
        {
            H_bar_ijab[ket][i].Qx = hhpp_channel_L[i].Qx;
            H_bar_ijab[ket][i].Qy = hhpp_channel_L[i].Qy;
            H_bar_ijab[ket][i].Qz = hhpp_channel_L[i].Qz;

            for(int ij = 0; ij < hhpp_channel_L[i].times; ij++)
            {
                kp_x = (int)(((double)hhpp_channel_L[i].Qx + 2 * hhpp_channel_L[i].wf[ij].q_kx) / 2.);
                kp_y = (int)(((double)hhpp_channel_L[i].Qy + 2 * hhpp_channel_L[i].wf[ij].q_ky) / 2.);
                kp_z = (int)(((double)hhpp_channel_L[i].Qz + 2 * hhpp_channel_L[i].wf[ij].q_kz) / 2.);
                sp = hhpp_channel_L[i].wf[ij].sz1;
                f_ii = f_ks_2(config, spstates_no, L, kp_x, kp_y, kp_z, sp, kp_x, kp_y, kp_z, sp, Minnesota_LECs);
                f_ii_kinetic = f_ks_1(config, spstates_no, L, kp_x, kp_y, kp_z, sp, kp_x, kp_y, kp_z, sp, Minnesota_LECs);

                kp_x = (int)(((double)hhpp_channel_L[i].Qx - 2 * hhpp_channel_L[i].wf[ij].q_kx) / 2.);
                kp_y = (int)(((double)hhpp_channel_L[i].Qy - 2 * hhpp_channel_L[i].wf[ij].q_ky) / 2.);
                kp_z = (int)(((double)hhpp_channel_L[i].Qz - 2 * hhpp_channel_L[i].wf[ij].q_kz) / 2.);
                sp = hhpp_channel_L[i].wf[ij].sz2;
                f_jj = f_ks_2(config, spstates_no, L, kp_x, kp_y, kp_z, sp, kp_x, kp_y, kp_z, sp, Minnesota_LECs);
                f_jj_kinetic = f_ks_1(config, spstates_no, L, kp_x, kp_y, kp_z, sp, kp_x, kp_y, kp_z, sp, Minnesota_LECs);
    
                for(int ab = 0; ab < hhpp_channel_R[i].times; ab++)
                {
                    H_bar_ijab[ket][i].matrix[ij][ab] = 0.;
                    h0 = V_ijab[i].matrix[ij][ab];
                    h1 = 0.;
                    h2 = 0.;
                    for(int kl = 0; kl< hh_channel[i].times; kl++)
                    {
                        h1 += 0.5 * t_ijab[ket][i].matrix[kl][ab] * V_ijkl[i].matrix[ij][kl];
                    }
    
                    for(int cd = 0; cd < pp_channel[i].times; cd++)
                    {
                        h2 += 0.5 * V_cdab[i].matrix[cd][ab] * t_ijab[ket][i].matrix[ij][cd];
                    }
    
                    kp_x = (int)(((double)hhpp_channel_R[i].Qx + 2 * hhpp_channel_R[i].wf[ab].q_kx) / 2.);
                    kp_y = (int)(((double)hhpp_channel_R[i].Qy + 2 * hhpp_channel_R[i].wf[ab].q_ky) / 2.);
                    kp_z = (int)(((double)hhpp_channel_R[i].Qz + 2 * hhpp_channel_R[i].wf[ab].q_kz) / 2.);
                    sp = hhpp_channel_R[i].wf[ab].sz1;
                    f_aa = f_ks_2(config, spstates_no, L, kp_x, kp_y, kp_z, sp, kp_x, kp_y, kp_z, sp, Minnesota_LECs);
                    f_aa_kinetic = f_ks_1(config, spstates_no, L, kp_x, kp_y, kp_z, sp, kp_x, kp_y, kp_z, sp, Minnesota_LECs);
                    kp_x = (int)(((double)hhpp_channel_R[i].Qx - 2 * hhpp_channel_R[i].wf[ab].q_kx) / 2.);
                    kp_y = (int)(((double)hhpp_channel_R[i].Qy - 2 * hhpp_channel_R[i].wf[ab].q_ky) / 2.);
                    kp_z = (int)(((double)hhpp_channel_R[i].Qz - 2 * hhpp_channel_R[i].wf[ab].q_kz) / 2.);
                    sp = hhpp_channel_R[i].wf[ab].sz2;
                    f_bb = f_ks_2(config, spstates_no, L, kp_x, kp_y, kp_z, sp, kp_x, kp_y, kp_z, sp, Minnesota_LECs);
                    f_bb_kinetic = f_ks_1(config, spstates_no, L, kp_x, kp_y, kp_z, sp, kp_x, kp_y, kp_z, sp, Minnesota_LECs);
                    h3 =  (f_aa + f_bb) * t_ijab[ket][i].matrix[ij][ab];
                    h4 = -(f_ii + f_jj) * t_ijab[ket][i].matrix[ij][ab];
  
                    k34 =  ((f_aa_kinetic + f_bb_kinetic) - (f_ii_kinetic + f_jj_kinetic)) * t_ijab[ket][i].matrix[ij][ab];
 
                    H_bar_ijab[ket][i].matrix[ij][ab] = h0 + h1 + h2 + h3 + h4;
                    kinetic_bar[ket][i].matrix[ij][ab] = k34;
                }
            }
        }
    }

}



