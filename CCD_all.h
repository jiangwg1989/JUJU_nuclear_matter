#define pi 3.1415926

const int magic_no = 14;
const int Nmax = 2;
const double mass = 938.92; // in unit of MeV
const double h_bar = 197.33;

//********************************************************************//
//****************************  STRUCTURE  ***************************//
//********************************************************************//
struct spstate
{
        int kx; //Nmax = 4, so kx, ky ,kz have (2*Nmax+1).
        int ky;
        int kz;
        int sz; //spin_up = +1, spin_down = -1;
        int occupied; //occupied = 1, unoccupied = 0;
};

struct tbwf//two-body wave function in relative momentum presentation
{
        int Q_kx;
        int Q_ky;
        int Q_kz;
        double q_kx;
        double q_ky;
        double q_kz;
        int sz1;
        int sz2;
        int flag;
};

struct CHANNEL
{
        int Qx; //record the same total momentum
        int Qy;
        int Qz;
        int times; //record how many times this Q occurs.
        tbwf *wf;
};

struct MATRIX_CHANNEL
{
        int Qx;
        int Qy;
        int Qz;
        double **matrix;
};

struct LECs
{
        double Vr;
        double Vs;
};

//********************************************************************//
//****************************   CLASSES   ***************************//
//********************************************************************//
class general_eigvalue
{
    public:
           general_eigvalue();
           ~general_eigvalue();
           void inoutput(double rho, LECs Minnesota_LECs, int subspace_dimension_in);
           int hhpp_dimension, hh_dimension, pp_no, hh_no, subspace_dimension;//subspace_dimension is the number of the new basis 
           double L;         
           double Ek;

    private:
           void setup_CCD_configuration_space(double rho);     
           void allocate_all_matrix(int subspace_dimension);
           void get_N_matrix();
           void get_H_matrix();
           void read_all_matrix();
           void calculate_H_bar_ij_ab(LECs Minnesota_LECs);
           void vacuum_expectation_value_H_bar(LECs Minnesota_LECs);
           void write_to_file();

           spstate *config;   
           int spstates_no;       
           tbwf *hh_config, *pp_config;
           CHANNEL *pp_channel, *hh_channel, *hhpp_channel_L, *hhpp_channel_R;
           double **temp_pointer1;//_V, **temp_pointer1_t, **temp_pointer1_H;
           double *temp_pointer2;//_V, *temp_pointer2_t, *temp_pointer2_H;
           MATRIX_CHANNEL **t_ijab, **H_bar_ijab, **l_abij, **x_ijab, **kinetic_bar;
           MATRIX_CHANNEL *V_ijkl, *V_cdab, *V_ijab;
           double *Ec;
           double **N_matrix, **H_matrix, **K_matrix;
           double *vacuum_H_bar;
};


//********************************************************************//
//****************************  FUNCTIONS  ***************************//
//********************************************************************//

double factorial(double n); // calculate the factorial
double combinator(double n, double m); //calculate the combinator
double kronecker(int s1, int s2);

double V_ks(double L, double q_rs_x, double q_rs_y, double q_rs_z, int sz_r, int sz_s, double q_pq_x, double q_pq_y, double q_pq_z, int sz_p, int sz_q, LECs Minnesota_LECs);

double V_ks_AS(double L, double q_rs_x, double q_rs_y, double q_rs_z, int sz_r, int sz_s, double q_pq_x, double q_pq_y, double q_pq_z, int sz_p, int sz_q, LECs Minnesota_LECs);

double f_ks(struct spstate *config, int spstates_no, double L, int kp_x, int kp_y, int kp_z, int sp, int kq_x, int kq_y, int kq_z, int sq, LECs Minnesota_LECs);

double f_ks_1(struct spstate *config, int spstates_no, double L, int kp_x, int kp_y, int kp_z, int sp, int kq_x, int kq_y, int kq_z, int sq, LECs Minnesota_LECs);

double f_ks_2(struct spstate *config, int spstates_no, double L, int kp_x, int kp_y, int kp_z, int sp, int kq_x, int kq_y, int kq_z, int sq, LECs Minnesota_LECs);
