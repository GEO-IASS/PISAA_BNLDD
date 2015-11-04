/*
 * Copyright 2014 Georgios Karagiannis
 *
 * Georgios Karagiannis
 * Postdoctoral research associate
 * Department of Mathematics, Purdue University
 * 150 N. University Street
 * West Lafayette, IN 47907-2067, USA
 *
 * Telephone: +1 765 494-3405
 *
 * Email: gkaragia@purdue.edu
 *
 * Contact email: georgios.stats@gmail.com
*/

#include <math.h>

#define MAX(x,y) (((x)>(y))?(x):(y))
#define MIN(x,y) (((x)<(y))?(x):(y))

double uniformrng(void) ;

int integerrng(int,int) ;

double cost(int*,int**,int,double*,int*,int) ;

void self_adj_index_search(int*,double,double*,int) ;

/* THIS IS THE TEMPORAL ORDER CHANGE METROPOLIS HASTINGS UPDATE */

void Mutation_TemporalOrderChange(int *z, int **zmat,
							double *fz, double *fz_cond,
							int N_node,
							int *changelength, int *changelist,
							double *theta, double *grid_points, int grid_size,
							double temp, double *accpr_pop,
							int *z_new, int **zmat_new, double *fz_cond_new) {

	int k_old ;
	int k_new ;
	int j ;
	int i ;
	int node ;

	double fz_new ;
	double rat ;
	double un ;

	/* get the partition label */

	self_adj_index_search(&k_old, *fz, grid_points, grid_size) ;

	/* generate the proposal */

	for (i=1; i<=N_node; i++) {
		z_new[i] = z[i] ;
		fz_cond_new[i] = fz_cond[i] ;
		for (j=1; j<=N_node; j++) zmat_new[i][j] = zmat[i][j] ;
	}

	node = integerrng(1, N_node-1) ;
	z_new[node] = z[node+1] ;
	z_new[node+1] = z[node] ;

    *changelength = 2 ;
    changelist[1] = node ; changelist[2] = node+1 ;
    for(j = node+2; j<=N_node; j++){
        if(zmat_new[node][j]==1 || zmat_new[node+1][j]==1){
           (*changelength)++ ;
           changelist[*changelength] = j ;
          }
        }

    fz_new = cost(z_new, zmat_new, N_node, fz_cond_new,
    		changelist, *changelength) ;

	/* get the partition label */

    self_adj_index_search(&k_new, fz_new, grid_points, grid_size) ;

	/* Accept reject */

    rat = -theta[k_new] -fz_new/temp +theta[k_old] + *fz/temp  ;
	*accpr_pop = exp( MIN(0.0, rat) ) ;

	un = uniformrng() ;
	if ( *accpr_pop>un ){
		*fz = fz_new;
		for (i=1; i<=N_node; i++) {
			z[i] = z_new[i] ;
			fz_cond[i] = fz_cond_new[i] ;
			for (j=1; j<=N_node; j++) zmat[i][j] = zmat_new[i][j] ;
		}
	}

}

/* THIS IS THE SKELETAL CHANGE METROPOLIS HASTINGS UPDATE */

void Mutation_SkeletalChange(int *z, int **zmat,
							double *fz, double *fz_cond,
							int N_node,
							int *changelength, int *changelist,
							double *theta, double *grid_points, int grid_size,
							double temp, double *accpr_pop,
							int *z_new, int **zmat_new, double *fz_cond_new) {

	int k_old ;
	int k_new ;
	int j ;
	int i ;
	int node1 ;
	int node2 ;

	double fz_new ;
	double rat ;
	double un ;

	/* get the partition label */

	self_adj_index_search(&k_old, *fz, grid_points, grid_size) ;

	/* generate the proposal */

	for (i=1; i<=N_node; i++) {
		z_new[i] = z[i] ;
		fz_cond_new[i] = fz_cond[i] ;
		for (j=1; j<=N_node; j++) zmat_new[i][j] = zmat[i][j] ;
	}

	node1 = integerrng(1, N_node) ;
	node2 = integerrng(1, N_node-1) ;
	if (node2 >= node1) node2++ ;
	if (node1 > node2) { i=node1; node1=node2; node2=i; }

	zmat_new[node1][node2] = 1-zmat_new[node1][node2] ;
	*changelength = 1 ;
	changelist[1] = node2 ;

    fz_new = cost(z_new, zmat_new, N_node, fz_cond_new,
    				changelist, *changelength) ;

	/* get the partition label */

    self_adj_index_search(&k_new, fz_new, grid_points, grid_size) ;

	/* Accept reject */

    rat = -theta[k_new] -fz_new/temp +theta[k_old] + *fz/temp  ;
	*accpr_pop = exp( MIN(0.0, rat) ) ;

	un = uniformrng() ;
	if ( *accpr_pop>un ){
		*fz = fz_new;
		for (i=1; i<=N_node; i++) {
			z[i] = z_new[i] ;
			fz_cond[i] = fz_cond_new[i] ;
			for (j=1; j<=N_node; j++) zmat[i][j] = zmat_new[i][j] ;
		}
	}

}

/* THIS IS THE DOUBLE SKELETAL CHANGE METROPOLIS HASTINGS UPDATE */

void Mutation_DoubleSkeletalChange(int *z, int **zmat,
							double *fz, double *fz_cond,
							int N_node,
							int *changelength, int *changelist,
							double *theta, double *grid_points, int grid_size,
							double temp, double *accpr_pop,
							int *z_new, int **zmat_new, double *fz_cond_new) {

	int k_old ;
	int k_new ;
	int j ;
	int i ;
	int nodeA1 ;
	int nodeA2 ;
	int nodeB1 ;
	int nodeB2 ;

	double fz_new ;
	double rat ;
	double un ;

	/* get the partition label */

	self_adj_index_search(&k_old, *fz, grid_points, grid_size) ;

	/* generate the proposal */

	for (i=1; i<=N_node; i++) {
		z_new[i] = z[i] ;
		fz_cond_new[i] = fz_cond[i] ;
		for (j=1; j<=N_node; j++) zmat_new[i][j] = zmat[i][j] ;
	}

	nodeA1 = 0 ;
	nodeA2 = 0 ;
	nodeB1 = 0 ;
	nodeB2 = 0 ;

	while (nodeA1 == nodeB1 && nodeA2 == nodeB2) {

		nodeA1 = integerrng(1, N_node) ;
		nodeA2 = integerrng(1, N_node-1) ;
		if (nodeA2 >= nodeA1) nodeA2++ ;
		if (nodeA1 > nodeA2) {i = nodeA1; nodeA1 = nodeA2; nodeA2 = i; }

		nodeB1 = integerrng(1, N_node) ;
		nodeB2 = integerrng(1, N_node-1) ;
		if (nodeB2 >= nodeB1) nodeB2++ ;
		if (nodeB1 > nodeB2) {i = nodeB1; nodeB1 = nodeB2; nodeB2 = i; }

	}

	zmat_new[nodeA1][nodeA2] = 1 -zmat_new[nodeA1][nodeA2];
	zmat_new[nodeB1][nodeB2] = 1 -zmat_new[nodeB1][nodeB2];

    if(nodeA2 == nodeB2) {
    	*changelength = 1 ;
    	changelist[1] = nodeA2 ;
    } else {
    	*changelength = 2 ;
    	changelist[1] = nodeA2 ;
    	changelist[2] = nodeB2 ;
    }

    fz_new = cost(z_new, zmat_new, N_node, fz_cond_new,
    						changelist, *changelength) ;

	/* get the partition label */

    self_adj_index_search(&k_new, fz_new, grid_points, grid_size) ;

	/* Accept reject */

    rat = -theta[k_new] -fz_new/temp +theta[k_old] + *fz/temp  ;
	*accpr_pop = exp( MIN(0.0, rat) ) ;

	un = uniformrng() ;
	if ( *accpr_pop>un ){
		*fz = fz_new;
		for (i=1; i<=N_node; i++) {
			z[i] = z_new[i] ;
			fz_cond[i] = fz_cond_new[i] ;
			for (j=1; j<=N_node; j++) zmat[i][j] = zmat_new[i][j] ;
		}
	}

}







