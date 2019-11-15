#include <iostream>
#include "stdio.h"
#include <stdlib.h>
#include <math.h>

using namespace std;


/*改int类型为long型*/
int count1Bits(int N) //Calculate the "1" counts of an integer written in form of binary form. 求一个整数写成二进制后共有多少个1，即全部核子数
{
	int n=0;
	if(N==0) {n=0;}
	while(N!=0)
	{
		N &= N-1;
		n++;
	}

	return n;
}

double factorial(double n) // calculate the factorial
{
	double result = n;
	do
	{
		result *= (n-1);
		n--;
	} while(n != 1);

	return result;
}

double combinator(double n, double m) //calculate the combinator
{
	double result;
	double numerator, denorminator;
	double temp_n;
	temp_n = n;
	numerator = n; //分子
	denorminator = factorial(m);//分母
	do
	{
		temp_n--;
		numerator *= temp_n;
	}while(temp_n >= (n-m+2));
	result = numerator / denorminator;
//	result = factorial(n) / (factorial(m) * factorial(n-m)); //too slow
	return result;
}

void print_bin(int N, int n_max) //the bin-representation, type an integer in form of binary form. 以二进制形式打印一个整数，即显示成bit-representation形式的组态
{
	double f = 0;
	int digits = 0;
	f = frexp(n_max, &digits);

    for(int i = digits-1; i >= 0; i--)
    {
        if( ( N & (1<<i) ) != 0) cout<<"1";
        else cout<<"0";
    }
    cout<<endl;
}

int judge_digit_no(int N)//判断一个数的最低位是第几位，这个数只有一个“1”或者两个“1”
{
	int i = 0;
	if( (N & 1) !=0 ) i=0;
	else
	{
		do
		{
			i++;
		}while( ((N & (1<<i)) == 0) );
	}
	return i;
}

//前后波函数只差一位的情况，给出前后两个波函数相差的那一位的位数
int delta_1_digit(int N1, int N2, int &i, int &j, int particle_no)//前后波函数只差一位的情况，给出前后两个波函数相差的那一位的位数
{
	int N;
	N = (N1 & N2);
//	cout<<"1bits have "<<count1Bits(0)<<endl;
	if(( particle_no - count1Bits(N) )  == 1)
	{
		i = judge_digit_no(N ^ N1);
		j = judge_digit_no(N ^ N2);
		return 1;
	}
	else {return 0;}
}

int reversebit(int N, int i) //N的i位取反
{
	return (N ^= (1<<i));
}

//前后波函数相差两位的情况，给出前后两个波函数相差的那两位的位数
int delta_2_digit(int N1, int N2, int &p1, int &q1, int &p2, int &q2, int particle_no)
{
	int N, NN1, NN2;
	N = (N1 & N2);
	if( (particle_no - count1Bits(N) )  == 2)
	{
		p1 = judge_digit_no(N ^ N1);
		NN1 = reversebit(N ^ N1,p1);
		q1 = judge_digit_no(NN1);//p1,q1 are the index in N1, p2 and q2 are the index in N2

		p2 = judge_digit_no(N ^ N2);
		NN2 = reversebit(N ^ N2,p2);
		q2 = judge_digit_no(NN2);
		return 1;
	}
	else {return 0;}
}

int count_between(int N, int i, int j) //count how many 1 between i bit and j bit, define i<j
{
	int count = 0;
	int ii, jj;
	if(i <= j) 
	{
		ii = i;
		jj = j;
	}
	if(i > j)
	{
		ii = j;
		jj = i;
	}
	for(int k = ii+1; k < jj; k++)
	{
		if( (N & (1<<k)) != 0) count++;
	}
	return count;
}

double two_body_diagonal(int N, int particle_no)
{
	double sum = 0.;
	double two_body_diag = 0.2; // actually <pq|v|qp> is not a constant.
	for(int i =0; i<particle_no-1; i++)
	{
		if( (N & (1<<i)) != 0 )
		{
			for(int j = i+1; j<particle_no; j++)
			{
				if( (N & (1<<j)) != 0 )
				{
					sum += two_body_diag;
				}
			}
		}
	}
	return sum;

}

/* Creat the configuration space. The reference state correspond to the N_max. All the other states are beyond */
/* this state, and at the same time the N values are smaller than N_max. So we use N-- here.*/

void creat_config(int particle_no, int degen_lev)//creat the full configurations.生成全部组态空间
{
	double n = combinator(2*degen_lev, particle_no);

//////////////////////////////////////////////////
//**********************************************//
//*  there is something wrong with this n!!!   *//
//**********************************************//	
//////////////////////////////////////////////////
	int *config = (int *) malloc (n*sizeof(int));
	int N_max;//N_max is corresponding value of the reference state (phi_0)
	N_max = pow(2, 2*degen_lev) - pow(2, 2*degen_lev - particle_no);//equals to 2^(2*P)-2^(2P-A). according to the geometric series sum formula.
	int N = N_max;
	int i=0;
	int number_1;

	double epsilon = 0.5; //epsilon is just a constant;
	double g_pq = -0.3; //is the coefficient of <q|p>.
	double g_ip_iq = -0.4; // v_iq1_iq2 is <iq1|v|iq2>. we assume it is a constant at first as well.
	double g_p1q1_p2q2 = -0.6; // v_p1q1_p2q2 is <p1q1|v|p2q2>. we assume it is a constant.
 	
 	/* Define the constant matrix of <q|p>, <ip|V|qi>, <p1,p2|q1,q2> */
 	double **e_pq = (double **)malloc(n*sizeof(double*));
 	double **v_ip_iq = (double**)malloc(n*sizeof(double*));
 	double **v_p1q1_p2q2 = (double**)malloc(n*sizeof(double*));
 	for(int i=0; i<n; i++)
	{
		e_pq[i] = (double*)malloc(n*sizeof(double));
		v_ip_iq[i] = (double*)malloc(n*sizeof(double));
		v_p1q1_p2q2[i] = (double*)malloc(n*sizeof(double));
	}

 	for(int i = 0; i < n; i++)
 	{
 		for(int j = 0; j < n; j++)
 		{
 			
 			if(j==i) e_pq[i][j] = ((i/2)+1) * epsilon; //s-p energies
 			else e_pq[i][j] = g_pq;
 		}
 	}


	for(int i = 0; i<n; i++)
	{
		for(int j = 0; j<n; j++)
		{
			if(j!=i) v_ip_iq[i][j] = g_ip_iq;
			else v_ip_iq[i][j] = 0.;
		}
	}


	for(int i = 0; i<n; i++)
	{
		for(int j = 0; j<n; j++)
		{
			if(j!=i) v_p1q1_p2q2[i][j] = g_p1q1_p2q2;
			else v_p1q1_p2q2[i][j] = 0.;
		}
	}


	///////// Generate the full configuration space. ////////
	do
	{
		number_1 = count1Bits(N);
		if(number_1==particle_no) // find the configurations have the same "1" occupied numbers who equal A.
		{
			config[i] = N;
			N--;
			i++;
		}
		else
		{
			N--;
		}
	}while(N != 0);



	for(int i = 0; i < n ; i++)
	{
		print_bin(config[i], config[0]);//print the configuration in form of binary form.
	}
	

	/****** Hamiltonian Matrix ******/

	/*declare a dynamic 2-D matrix, n × n dimension,声明2维哈密顿量矩阵*/
	double **Matrix_element;
	Matrix_element = (double **)malloc(n*sizeof(double*));
	for(int i=0; i<n; i++)
	{
		Matrix_element[i] = (double*)malloc(n*sizeof(double));
	}

	//initial the Hamiltonian Matrix elements.
	for(int i =0; i<n; i++)
	{
		for(int j =0; j<n; j++)
		{
			Matrix_element[i][j]=0.;
		}
	}

	/* 1-body operator, i=j*/
	for(int i = 0; i<n; i++)
	{
		Matrix_element[i][i] += e_pq[i][i];
	}
	

	/* 1-body operator, i!=j*/
	for(int i =0; i<n; i++)
	{
		int config_temp_i = config[i];
		for(int j =i+1; j<n; j++)
		{
			int config_temp_j = config[j];
			int p, q;//the different index p and q in the i and j states, respectively.
			int flag =0;
			int N_temp, Permut;
			flag = delta_1_digit(config_temp_i,config_temp_j,p,q,particle_no);
			if(flag == 1) 
			{
				N_temp = (config_temp_i | (1<<q));
				Permut = count_between(N_temp, p, q);
				Matrix_element[i][j] = pow(-1, Permut) * e_pq[i][j]; //!!!Here h0 actually should not be a constant. in fact it is the integral value of <p|h0|q>.
			}
		}
	}


	/* 2-body operator, diagonal matrix elements */
	for(int i = 0; i <n; i++)
	{
		Matrix_element[i][i] += two_body_diagonal(config[i], particle_no);
	}


	/* 2-body operator,  |ip> and |iq> */
	for(int i = 0; i<n; i++)
	{
		int config_temp_i1 = config[i];
		for(int j =0; j<n ;j++)
		{
			int config_temp_j1 = config[j];
			int p,q; //different index p,q.
			int flag1 = 0;
			int N_temp1, Permut1;
			flag1 = delta_1_digit(config_temp_i1, config_temp_j1,p,q,particle_no);
			if(flag1 == 1)
			{
				N_temp1 = (config_temp_i1 | (1<<q));
				Permut1 = count_between(N_temp1, p, q);
				Matrix_element[i][j] += (particle_no-2) * pow(-1, Permut1) * v_ip_iq[i][j]; 
			}
		}
	}


	/* 2-body operator, |p1q1> and |p2q2> */

	for(int i = 0; i<n; i++)
	{
		int config_temp_i2 = config[i];
		for(int j =0; j<n ;j++)
		{
			int config_temp_j2 = config[j];
			int p1, q1, p2, q2; //p1 and q1 are in N1, p2 and q2 are in N2.
			int flag2 = 0;
			int N_temp2p, Permut2p, N_temp2q, Permut2q;//2p是p1和p2之间相差的1的个数，2q是q1和q2位之间相差的1的个数。
			flag2 = delta_2_digit(config_temp_i2, config_temp_j2,p1,q1,p2,q2,particle_no);
			if(flag2 == 1)
			{
				N_temp2p = (config_temp_i2 | (1<<p2));
				Permut2p = count_between(N_temp2p, p1, p2);

				N_temp2q = (config_temp_i2 | (1<<q2));
				Permut2q = count_between(N_temp2q, q1, q2);

				Matrix_element[i][j] = pow(-1, Permut2p) * pow(-1, Permut2q) * v_p1q1_p2q2[i][j]; //这里的系数需要再仔细推导！！！
			}
		}
	}
	
	
	for(int i =0; i<n; i++)
	{
		for(int j=0; j<n; j++)
		{
			printf("%.2f ",Matrix_element[i][j]);
		}
		cout<<endl;
	}

	
	////////// FREE ///////////
	for(int i=0; i<n; i++)
	{
		free(Matrix_element[i]);
		free(e_pq[i]);
		free(v_ip_iq[i]);
		free(v_p1q1_p2q2[i]);
	}
	free(Matrix_element);
	free(v_p1q1_p2q2);
	free(e_pq);
	free(v_ip_iq);

	free(config);

}



int main()
{
	creat_config(4,4);
//	cout<<combinator(32,16)<<endl;

	return 0;
}
