/*
 * Georgios Karagiannis 
 * Postdoctoral research associate
 * Department of Mathematics, Purdue University
 * 150 N. University Street
 * West Lafayette, IN 47907-2067, USA
 *
 * Telephone: +1 (765) 496-1007
 *
 * Email: gkaragia@purdue.edu
 *
 * Contact email: georgios.stats@gmail.com
 *
 * Georgios Karagiannis © 2014  
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "nrutil.h"

void setseedrng(unsigned long) ;
double uniformrng(void) ;
int integerrng(int,int) ;
void permutrng(int*,int) ;

void get_data(char[],int*) ;
double cost(int*,int**,int,double*,int*,int) ;

void self_adj_grid_points(double*,int,double,double) ;
void self_adj_desired_freq(double*,int,double) ;
void self_adj_theta_update(double*,int,double*,double*,int,double*,double,double*) ;
void self_adj_theta_norm(double*,double*,double*,int,double) ;

void Mutation_TemporalOrderChange(int*,int**,double*,double*,int,
		int*,int*,double*,double*,int,double,double*,int*,int**,double*) ;
void Mutation_SkeletalChange(int*,int**,double*,double*,int,
		int*,int*,double*,double*,int,double,double*,int*,int**,double*) ;
void Mutation_DoubleSkeletalChange(int*,int**,double*,double*,int,
		int*,int*,double*,double*,int,double,double*,int*,int**,double*) ;

void Crossover_int_Kpoint(int**,int***,double *,double **,int,int ,int*,int**,
							double*,double*,int, double,double*,
							int**,int**,int*,double*,double*) ;

double mean_vec(double *x, int d) {

	int i ;
	double x_mean ;

	x_mean = 0.0 ;
	for (i=1; i<=d; i++) x_mean += x[i] ;

	return x_mean/d ;
}

void flags_usage(void) {
	printf("Usage:\n") ;
	printf(" -ID value \n") ;
	printf(" -Niter value \n") ;
	printf(" -Npop value \n") ;
	printf(" -Data value \n") ;
	printf(" -Gwarm value \n") ;
	printf(" -Ghigh value \n") ;
	printf(" -Gpow value \n") ;
	printf(" -Hlow value \n") ;
	printf(" -Hhigh value \n") ;
	printf(" -Hsize value \n") ;
	printf(" -Hzeta value \n") ;
	printf(" -Hconst value \n") ;
	printf(" -Twarm value \n") ;
	printf(" -Tlow value \n") ;
	printf(" -Thigh value \n") ;
	printf(" -Tlow value \n") ;
	printf(" -Tpow value \n") ;
	printf(" -flags \n") ;
	exit (8) ;
}

void update_best_value(int *z_best, int **zmat_best,
						double *fz_best, double *fz_cond_best,
						int **x, int ***xmat,
						double *fx, double **fx_cond,
						int N_node, int N_population ){
	int n ;
	int i ;
	int j ;
	int n_best ;

	n_best = 0 ;
	for (n=1; n<=N_population; n++)
		if (fx[n] < *fz_best) {
			n_best = n ;
			*fz_best = fx[n] ;
		}

	if (n_best > 0)
		for(i=1; i<=N_node; i++) {
			z_best[i] = x[n_best][i] ;
			for(j=1; j<=N_node; j++)  zmat_best[i][j] = xmat[n_best][i][j] ;
			fz_cond_best[i] = fx_cond[n_best][i] ;
		}
}


#define MAX(x,y) (((x)>(y))?(x):(y))
#define MIN(x,y) (((x)<(y))?(x):(y))

int main(int argc, char *argv[]){

	double sumv ;
	double avev ;
	double minv ;
	double maxv ;
	double un ;

	int N_empty ;

	int lab_mcmc ;
	int i ;
	int k ;
	int j ;
	int pop ;
	int iter ;
	int ID ;
	int N_sample ;
	int N_thinning ;
	int N_population ;
	int N_iteration ;

	int N_node ;

	double z_min ;
	double z_max ;

	double *freq_ref ;
	double *freq_est ;
	double freq_zeta ;
	double *theta ;
	double theta_norm_const ;
	double *grid_points ;
	int grid_size ;
	double grid_high ;
	double grid_low ;

	double gain ;
	double gain_pow ;
	int gain_warmup ;
	double gain_high ;

	double temp_saa ;
	double temp_saa_pow ;
	int temp_saa_warmup ;
	double temp_saa_high ;
	double temp_saa_low ;

	int **x ;
	int ***xmat ;
	double *fx ;
	double **fx_cond ;

	int *changelength ;
	int **changelist ;

	int *z_best ;
	int **zmat_best ;
	double fz_best ;
	double *fz_cond_best ;

	double accpr ;

	double *work0_d1 ;
	int *work1_i1 ;
	double *work2_d1 ;
	int **work3_i2 ;
	double *work4_d1 ;
	int **work5_i2 ;

	double exp_accpr_mh1 ;
	double exp_accpr_mh2 ;
	double exp_accpr_mh3 ;
	double exp_accpr_co1 ;

	int n_mh1 ;
	int n_mh2 ;
	int n_mh3 ;
	int n_co1 ;

	FILE *ins_hist = NULL ;
	FILE *ins_fz_best_trace = NULL ;
	FILE *ins_fz_best = NULL ;

	char file_name[50] ;
	char datapath[50] ;

	/*INITIALIZE THE RNG -------------------------------------------------- */
	printf("\n INITIALIZE THE RNG \n") ;

	setseedrng( (unsigned long) time(NULL) ) ;
	for ( i=1 ; i<=10 ; i++ ) un = uniformrng() ;

	/* SET DEFAULT ALGORITHMIC SETTINGS ----------------------------------- */
	printf("\n SET ALGORITHMIC SETTINGS \n") ;

	snprintf(datapath, sizeof datapath, "./SPECT.dat") ;

	N_population = 5 ;
	N_iteration = 200000000 ;
	ID = 1 ;
	N_sample = 20000 ;
	N_thinning = (int) (N_iteration/N_sample) ;

	gain_warmup = 500000 ;
	gain_pow = 1.0 ;
	gain_high = 1.0 ;

	temp_saa_warmup = 500000 ;
	temp_saa_pow = 0.5 ;
	temp_saa_high = 1.0 ;
	temp_saa_low = 0.1 ;

	freq_zeta = 0.1 ;
	grid_low = 2000.0 ;
	grid_high = 3999.0 ;
	grid_size = 2000 ;

	theta_norm_const = 100.0 ;

	/*PASS EXTERNAL ALGORITHMIC ARGUMENTS -------------------------------- */

	for (i = 1; i < argc; i++) {
		if (strcmp("-ID", argv[i]) == 0)					/*..REPEAT NO*/
			ID = atoi(argv[++i]) ;
		else if (strcmp("-Data", argv[i]) == 0) 			/*..DATA*/
			snprintf(datapath, sizeof datapath, argv[++i]) ;
		else if (strcmp("-Niter", argv[i]) == 0)			/*..COUNTERS*/
		{
			N_iteration = atoi(argv[++i]) ;
			N_thinning = (N_sample==0)?1:(N_iteration/N_sample) ;
		}
		else if (strcmp("-Npop", argv[i]) == 0)
			N_population = atoi(argv[++i]) ;
		else if (strcmp("-Nsam", argv[i]) == 0)
		{
			N_sample = atoi(argv[++i]) ;
			N_thinning = (N_sample==0)?1:(N_iteration/N_sample) ;
		}
		else if (strcmp("-Gwarm", argv[i]) == 0)			/*..GAIN*/
			gain_warmup = atoi(argv[++i]) ;
		else if (strcmp("-Ghigh", argv[i]) == 0)
			gain_high = atof(argv[++i]) ;
		else if (strcmp("-Gpow", argv[i]) == 0)
			gain_pow = atof(argv[++i]) ;
		else if (strcmp("-Hlow", argv[i]) == 0)				/*..HIST*/
			grid_low = atof(argv[++i]) ;
		else if (strcmp("-Hhigh", argv[i]) == 0)
			grid_high = atof(argv[++i]) ;
		else if (strcmp("-Hsize", argv[i]) == 0)
			grid_size = atoi(argv[++i]) ;
		else if (strcmp("-Hzeta", argv[i]) == 0)
			freq_zeta = atof(argv[++i]) ;
		else if (strcmp("-Hconst", argv[i]) == 0)
			theta_norm_const = atof(argv[++i]) ;
		else if (strcmp("-Twarm", argv[i]) == 0) 		/*..TEMPERATURE LADER*/
			temp_saa_warmup = atoi(argv[++i]) ;
		else if (strcmp("-Tlow", argv[i]) == 0)
			temp_saa_low = atof(argv[++i]) ;
		else if (strcmp("-Thigh", argv[i]) == 0)
			temp_saa_high = atof(argv[++i]) ;
		else if (strcmp("-Tpow", argv[i]) == 0)
			temp_saa_pow = atof(argv[++i]) ;
		else if (strcmp("-flags", argv[i]) == 0)
			flags_usage() ;
   }

	/*PRINT INPUTS OF THE ALGORITHM*/
    printf( "\n PRINT INPUTS OF THE ALGORITHM \n") ;

	printf( "Data path: \t\t %s \n", datapath) ;

	printf( "N_iteration: \t\t %d \n", N_iteration) ;
	printf( "N_population: \t\t %d \n", N_population) ;
	printf( "ID: \t\t\t %d \n", ID) ;
	printf( "N_sample: \t\t %d \n", N_sample) ;
	printf( "N_thinning: \t\t %d \n", N_thinning) ;

	printf( "gain_warmup: \t\t %d \n", gain_warmup) ;
	printf( "gain_high: \t\t %f \n", gain_high) ;
	printf( "gain_pow: \t\t %f \n", gain_pow) ;

	printf( "grid_low: \t\t %f \n", grid_low) ;
	printf( "grid_high: \t\t %f \n", grid_high) ;
	printf( "grid_size: \t\t %d \n", grid_size) ;
	printf( "freq_zeta: \t\t %f \n", freq_zeta) ;
	printf( "theta_norm_const: \t %f \n", theta_norm_const) ;

	printf( "temp_saa_warmup: \t %d \n", temp_saa_warmup) ;
	printf( "temp_saa_low: \t\t %f \n", temp_saa_low) ;
	printf( "temp_saa_high: \t\t %f \n", temp_saa_high) ;
	printf( "temp_saa_pow: \t\t %f \n", temp_saa_pow) ;

	printf( "\n") ;

	/*GET DATA -----------------------------------------------------------*/
	printf("\n GET DATA\n") ;
	get_data(datapath,&N_node) ;

	/*OPEN FILES -----------------------------------------------------------*/
	printf("\n OPEN FILES\n") ;

	snprintf(file_name, sizeof file_name, "./output_files/hist-n=%d-r=%d",
			N_population, ID);
    ins_hist = fopen( file_name , "w" ) ;

    if ( N_sample>0 ) {
    snprintf(file_name, sizeof file_name, "./output_files/fz_best_trace-n=%d-r=%d",
    		N_population, ID);
    ins_fz_best_trace = fopen( file_name , "w" ) ;
    }

    snprintf(file_name, sizeof file_name, "./output_files/fz_best-n=%d-r=%d",
    		N_population, ID);
    ins_fz_best = fopen( file_name , "w" ) ;

	/* ALLOCATE SPACE FOR THE ARRAYS -------------------------------------- */
	printf("\n ALLOCATE SPACE FOR THE ARRAYS \n") ;

	grid_points = dvector( 1 , grid_size ) ;
	freq_ref = dvector( 1 , grid_size+1 ) ;
	theta = dvector( 1 , grid_size+1 ) ;
	freq_est = dvector( 1 , grid_size+1 ) ;

	x = imatrix(1, N_population, 1, N_node) ;
	xmat = i3tensor(1, N_population, 1, N_node, 1, N_node) ;
	fx = dvector(1, N_population) ;
	fx_cond = dmatrix(1, N_population, 1, N_node) ;

	changelength = ivector(1, N_population) ;
	changelist = imatrix(1, N_population, 1, N_node) ;

	z_best = ivector( 1, N_node ) ;
	zmat_best = imatrix(1, N_node, 1, N_node) ;
	fz_cond_best = dvector( 1, N_node ) ;

	work0_d1 = dvector(1, grid_size+1) ;
	work1_i1 = ivector(1, N_node) ;
	work2_d1 = dvector(1, N_node) ;
	work3_i2 = imatrix(1, N_node, 1, N_node) ;
	work4_d1 = dvector(1, N_node) ;
	work5_i2 = imatrix(1, N_node, 1, N_node) ;

	/* INITIALIZE THE POPULATION -------------------------------------- */
	printf("\n INITIALIZE THE POPULATION \n") ;

	/*Get a random initialization*/
	for (pop=1 ; pop<=N_population; pop++) {

		/*Initialize values*/

		for (i=1;i<=N_node; i++) x[pop][i] = i ;
		permutrng(x[pop], N_node) ;

		for (i=1;i<=N_node; i++)
			for (j=1;j<=N_node; j++)
				xmat[pop][i][j] = ((i==j) ? 1 : 0) ;

		changelength[pop] = N_node ;
		for (i=1;i<=N_node; i++) changelist[pop][i] = i ;

		fx[pop] = cost(x[pop], xmat[pop], N_node, fx_cond[pop],
				changelist[pop], changelength[pop]) ;

	}

	for (pop=1 ; pop<=N_population; pop++) {

		printf("Popul.  : %d \n",pop) ;

		printf("Cost    : %f \n",fx[pop]) ;

		printf("Cost[i] : ") ;

		for (i=1 ; i<=N_node; i++) printf("%f ", fx_cond[pop][i]) ;
		printf("\n") ;

		printf("x[i]    : ") ;
		for (i=1;i<=N_node; i++) printf("%d ",x[pop][i]) ; printf("\n") ;

		printf("xmat    : \n") ;
		for (i=1; i<=N_node; i++){
			for (j=1; j<=N_node; j++) printf("%d ",xmat[pop][i][j]) ;
			printf("\n") ;
		}
	}
	printf("\n") ;

	/* INITIALIZE BEST VALUES ----------------------------------------- */
	printf("\n INITIALIZE BEST VALUES \n") ;

	k = 1 ; fz_best = fx[1] ;
	for ( pop=2 ; pop <= N_population ; pop++ )
		if ( fx[pop] < fz_best ){
			k = pop ;
			fz_best = fx[k] ;
		}
	for (i=1; i<=N_node; i++) fz_cond_best[i] = fx_cond[k][i] ;
	for (i=1; i<=N_node; i++) z_best[i] = x[k][i] ;
	for (i=1; i<=N_node; i++)
		for (j=1; j<=N_node; j++)
			zmat_best[i][j] = xmat[k][i][j] ;

	/*PRINT*/

	printf("Best value : \n") ;

	printf( "fz_best : %f \n",fz_best) ;

	printf( "fz_best_cond  : \n") ;
	for (i=1; i<=N_node; i++) printf("%f ", fz_cond_best[i]) ;
	printf("\n") ;

	printf( "z_best  : \n") ;
	for ( i=1 ; i<=N_node ; i++ ) printf("%d ", z_best[i]) ;
	printf("\n") ;

	printf( "zmat_best  : \n") ;
	for ( i=1 ; i<=N_node ; i++ ) {
		for ( j=1 ; j<=N_node ; j++ ) printf("%d ", zmat_best[i][j]) ;
		printf("\n") ;
	}

	/* INITIALIZE THE SELF-ADJUSTING MECHANISM ------------------------- */
	printf("\n INITIALIZE THE SELF-ADJUSTING MECHANISM \n") ;

	/* INITIALIZE THE REFERENCE FREQUENCE */
	printf("\n INITIALIZE THE REFERENCE FREQUENCE \n") ;

	self_adj_desired_freq(freq_ref, grid_size, freq_zeta) ;

	/* INITIALIZE EMPIRICAL FREQUENCE */
	printf("\n INITIALIZE EMPIRICAL FREQUENCE \n") ;

	for (i=1 ; i<=grid_size+1 ; i++) freq_est[i] = 0.0 ;

	/* INITIALIZE GRID POINTS */
	printf("\n INITIALIZE GRID POINTSE \n") ;

	self_adj_grid_points(grid_points, grid_size, grid_low, grid_high) ;

	/* INITIALIZE THETA */
	printf("\n INITIALIZE THETA \n") ;
	for (i=1 ; i<=grid_size+1 ; i++) theta[i] = 0.0 ;

	/* 	PRINT INITIAL SETTINGS OF THE SELF-ADJUSTMENT PROSEDURE */
	if( grid_size < 500 ) {
		printf("\n PRINT INITIAL SETTINGS OF "
				"THE SELF-ADJUSTMENT PROSEDURE \n") ;
		for (i=1 ; i<=grid_size ; i++)
			printf("%d \t %f \t %f \t %f \t %f \t %f \n",
					i,
					grid_points[i],
					theta[i], exp(theta[i]),
					freq_ref[i], freq_est[i] ) ;
		i = grid_size+1 ;
		printf("%d \t %f \t %f \t %f \t %f \t %f \n",
				i,
				0.0,
				theta[i], exp(theta[i]),
				freq_ref[i], freq_est[i] ) ;
	}

	/* INITIALIZE THE COUNTERS ---------------------------------------- */
	printf("\n INITIALIZE THE COUNTERS \n") ;

	/* INITIALIZE THE EXPETED ACCEPTANCE PROBABLITY COUNTERS */

	exp_accpr_mh1 = 0.0 ; n_mh1 = 0 ;
	exp_accpr_mh2 = 0.0 ; n_mh2 = 0 ;
	exp_accpr_mh3 = 0.0 ; n_mh3 = 0 ;
	exp_accpr_co1 = 0.0 ; n_co1 = 0 ;

	/* PERFORM PISAA ITERATIONS ------------------------------------------*/
	printf("\n PERFORM PISAA ITERATIONS \n") ;

	for( iter=-100; iter<=N_iteration; iter++ ) {

		/* UPDATE THE TEMPERATURE -------------------------------------- */

		if( iter <= temp_saa_warmup )
			temp_saa = temp_saa_high + temp_saa_low ;
		else{
			temp_saa = temp_saa_warmup/((double) iter) ;
			temp_saa = temp_saa_high*pow(temp_saa,temp_saa_pow)
							+ temp_saa_low ;
		}

		/* SAMPLING UPDATE -------------------------------------------- */

		lab_mcmc = integerrng(1, 3) ;

		switch ( lab_mcmc ) {
		case 0 :
			break ;
		case 1 :
			/*sample*/
			sumv = 0.0 ;
			for (pop=1 ; pop<=N_population ; pop++ ) {
				Mutation_TemporalOrderChange(x[pop], xmat[pop],
										&fx[pop], fx_cond[pop],
										N_node,
										&changelength[pop], changelist[pop],
										theta, grid_points, grid_size,
										temp_saa, &accpr,
										work1_i1, work3_i2, work2_d1) ;
				sumv += accpr ;
			}
			accpr = sumv/pop ;
			exp_accpr_mh1 += accpr ; n_mh1++ ;
			break ;
		case 2 :
			/*sample*/
			sumv = 0.0 ;
			for (pop=1 ; pop<=N_population ; pop++ ) {
				Mutation_SkeletalChange(x[pop], xmat[pop],
										&fx[pop], fx_cond[pop],
										N_node,
										&changelength[pop], changelist[pop],
										theta, grid_points, grid_size,
										temp_saa, &accpr,
										work1_i1, work3_i2, work2_d1) ;
				sumv += accpr ;
			}
			accpr = sumv/pop ;
			exp_accpr_mh2 += accpr ; n_mh2++ ;
			break ;
		case 3 :
			/*sample*/
			sumv = 0.0 ;
			for (pop=1 ; pop<=N_population ; pop++ ) {
				Mutation_DoubleSkeletalChange(x[pop], xmat[pop],
										&fx[pop], fx_cond[pop],
										N_node,
										&changelength[pop], changelist[pop],
										theta, grid_points, grid_size,
										temp_saa, &accpr,
										work1_i1, work3_i2, work2_d1) ;
				sumv += accpr ;
			}
			accpr = sumv/pop ;
			exp_accpr_mh3 += accpr ; n_mh3++ ;
			break ;
		case 4 :
			/*sample*/
			Crossover_int_Kpoint(x, xmat,
										fx, fx_cond,
										N_node, N_population,
										changelength, changelist,
										theta, grid_points, grid_size,
										temp_saa, &accpr,
										work3_i2, work5_i2, work1_i1,
										work2_d1, work4_d1) ;
			exp_accpr_co1 += accpr ; n_co1++ ;
			break ;
		default :
			break ;
		}

		/* UPDATE THE GAIN FUNCTION ------------------------------------ */

		if( iter <= gain_warmup )
			gain = gain_high ;
		else {
			gain = gain_warmup/((double) iter) ;
			gain = gain_high*pow( gain,gain_pow) ;
		}

		/* WEIGHTING UPDATE -------------------------------------------- */

		self_adj_theta_update(fx, N_population,
								theta, grid_points, grid_size,
								freq_ref, gain, work0_d1) ;

		for (i=1;i<=grid_size+1;i++) freq_est[i] += work0_d1[i] ;

		/* BEST VALUE UPDATE ------------------------------------------- */

		update_best_value(z_best, zmat_best, &fz_best, fz_cond_best,
							x, xmat, fx, fx_cond, N_node, N_population ) ;


		/* KEEP RECORDS ------------------------------------------------ */

		if ( iter > 0 && iter%N_thinning == 0  && N_sample>0)
					printf("iter=%d gain=%f temp_saa=%f min=%f \n",
							iter, gain, temp_saa, fz_best ) ;

		if (ins_fz_best_trace != NULL)
		if (iter>0 && iter%N_thinning==0 && N_sample>0) {
			fprintf(ins_fz_best_trace,"%d %f \n", iter, fz_best) ;
			fflush(ins_fz_best_trace) ;
		}

	} /* ... ITERATIONS*/

	/*THETA-NORMALIZATION STEP ----------------------------------------- */
	printf("\n THETA-NORMALIZATION STEP \n") ;

	for ( i=1 ; i<=grid_size+1 ; i++ ) freq_est[i] /= N_iteration ;

	self_adj_theta_norm(theta,freq_ref,freq_est,grid_size,theta_norm_const) ;

	/* RECORD THE SELF-ADJUSTMENT MECHANISM */

	if (ins_hist != NULL){

		for ( i=1 ; i<=grid_size ; i++ )
			fprintf(ins_hist, "%d %f %f %15.15f %f %f \n",
					i,grid_points[i],
					theta[i], exp(theta[i]),
					freq_ref[i], freq_est[i]) ;

		i = grid_size+1 ;
		fprintf(ins_hist, "%d %f %f %f %f %f \n",
				i,0.0,
				theta[i], exp(theta[i]),
				freq_ref[i], freq_est[i]) ;

		fflush(ins_hist) ;
	}

	/* RECORD BEST VALUES ---------------------------------------------- */
	printf("\n RECORD BEST VALUES \n") ;

	if (ins_fz_best != NULL){

		fprintf(ins_fz_best, "%f \n", fz_best) ;

		for (i=1; i<=N_node; i++)
			fprintf(ins_fz_best, "%f ", fz_cond_best[i]) ;
		fprintf(ins_fz_best, "\n") ;

		for (i=1; i<=N_node; i++)
			fprintf(ins_fz_best, "%d ", z_best[i]) ;
		fprintf(ins_fz_best, "\n") ;

		for (i=1; i<=N_node; i++) {
			for (j=1; j<=N_node; j++)
				fprintf(ins_fz_best, "%d ", zmat_best[i][j]) ;
			fprintf(ins_fz_best, "\n") ;
		}

		fflush(ins_fz_best) ;
	}

	/* SUMMARY --------------------------------------------------------- */

	printf("\n SUMMARY \n") ;

	/*Result*/
	printf( "\n ...Result \n");

	printf( "fz_best : %f \n",fz_best) ;

	printf( "fz_best_cond  : \n") ;
	for (i=1; i<=N_node; i++) printf("%f ", fz_cond_best[i]) ;
	printf("\n") ;

	printf( "z_best  : \n") ;
	for ( i=1 ; i<=N_node ; i++ ) printf("%d ", z_best[i]) ;
	printf("\n") ;

	printf( "zmat_best  : \n") ;
	for ( i=1 ; i<=N_node ; i++ ) {
		for ( j=1 ; j<=N_node ; j++ ) printf("%d ", zmat_best[i][j]) ;
		printf("\n") ;
	}

	/*Acceptance ratios*/
	printf( "\n ...Acceptance ratios \n");

	exp_accpr_mh1 /= n_mh1 ;
	exp_accpr_mh2 /= n_mh2 ;
	exp_accpr_mh3 /= n_mh3 ;
	exp_accpr_co1 /= n_co1 ;

	printf( "mh1 rate=%f \n",exp_accpr_mh1);
	printf( "mh2 rate=%f \n",exp_accpr_mh2);
	printf( "mh3 rate=%f \n",exp_accpr_mh3);
	printf( "co1 rate=%f \n",exp_accpr_co1);

	/* CLOSE FILES -------------------------------------------------------- */

	printf("\n CLOSE FILES \n") ;

	if (ins_hist != NULL) fclose( ins_hist ) ;
	if (ins_fz_best_trace != NULL) fclose( ins_fz_best_trace ) ;
	if (ins_fz_best != NULL) fclose( ins_fz_best ) ;

	/* FREE SPACE --------------------------------------------------------- */

	printf("\n FREE SPACE \n") ;

	free_dvector(grid_points, 1, grid_size) ;
	free_dvector(freq_ref, 1, grid_size+1) ;
	free_dvector(freq_est, 1, grid_size+1) ;
	free_dvector(theta, 1, grid_size+1) ;

	free_dvector(fx, 1, N_population) ;
	free_imatrix(x, 1, N_population, 1, N_node) ;
	free_i3tensor(xmat, 1, N_population, 1, N_node, 1, N_node) ;
	free_dmatrix(fx_cond, 1, N_population, 1, N_node) ;

	free_ivector(changelength, 1, N_population) ;
	free_imatrix(changelist, 1, N_population, 1, N_node) ;

	free_ivector(z_best, 1, N_node) ;
	free_imatrix(zmat_best, 1, N_node, 1, N_node) ;
	free_dvector(fz_cond_best, 1, N_node ) ;

	free_dvector(work0_d1, 1, grid_size+1) ;
	free_ivector(work1_i1, 1, N_node) ;
	free_dvector(work2_d1, 1, N_node) ;
	free_imatrix(work3_i2, 1, N_node, 1, N_node) ;
	free_dvector(work4_d1, 1, N_node) ;
	free_imatrix(work5_i2, 1, N_node, 1, N_node) ;

	printf("\n\nDONE\n") ;

	return 0 ;
}






