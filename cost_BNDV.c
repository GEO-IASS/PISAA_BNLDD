/*
 * Georgios Karagiannis (2014)
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

/*
 * A modified version of a routine borrowed from the source code of the book
 *
 * Liang, Faming, Chuanhai Liu, and Raymond Carroll.
 * Advanced Markov chain Monte Carlo methods: learning from past samples.
 * Vol. 714. John Wiley & Sons, 2011.
 *
 * Chapter 7
 *
 * */

#include <stdlib.h>
#include <math.h>
#include "nrutil.h"
#include <stdio.h>

#ifndef INFINITY
	#include <float.h>
	#define INFINITY DBL_MAX
#endif

double lgamma(double) ;
#define gammln(X) (lgamma((double)(X)))

static int **data = NULL ;
static int N_data ;
static int *state = NULL ;

void get_data(char datapath[], int *N_node){

	int i ;
	int j ;
	FILE *ins ;

	/*GET NAME OF THE FILE*/

	ins = fopen(datapath, "r") ;

	if (ins==NULL) {
		printf("No data set\n") ;
		abort() ;
	}

	/*GET THE DATASET*/

	fscanf(ins,"%d %d ", &N_data, N_node)  ;

	data = imatrix(1, N_data, 1, *N_node) ;
	state = ivector(1, *N_node) ;

	for(i=1; i<=N_data; i++)
	    for(j=1; j<=*N_node; j++)
	    	fscanf(ins,"%d ",&data[i][j]) ;

	/*for(i=1; i<=N_data; i++)
	    for(j=1; j<=*N_node; j++)
	    	data[i][j]++ ;*/

	for(i=1; i<=*N_node; i++) state[i] = 2 ;

	fclose(ins) ;
}

double cost(int *x, int **mat, int node_num, double *fvalue,
		int *changelist, int changelength) {

	int **tab, **tab00;
	int count, count00, num, num00, i, j, k, l, m, s, tep;
	int numparent, parstate, *parlist;
	double sum, alphaijk, alphaik;

	int limparent = 5 ;
	double prior_alpha = 1.0 ;
	double prior_gamma = 0.1111 ;

	parlist = (int *) malloc((size_t) (node_num*sizeof(int))) ;

for(m=1; m<=changelength; m++){

    i=changelist[m];

    parstate=1; s=0;
    for(j=1; j<i; j++){
        if(mat[j][i]==1){
           parstate*=state[x[j]];
           s++;
           parlist[s]=x[j];
         }
       }
    numparent=s;
    fvalue[i]=numparent*log(prior_gamma);

    if(numparent>limparent){
       for(j=1; j<=node_num; j++) fvalue[j]=1.0e+100;
       goto ABC;
      }

    /* parent state table */
    tab=imatrix(1,(parstate+1)*state[x[i]],1,2);
    tab00=imatrix(1,parstate+1,1,2);

    for(j=1; j<=parstate*state[x[i]]; j++)
       for(l=1; l<=2; l++) tab[j][l]=0;
    for(j=1; j<=parstate; j++)
       for(l=1; l<=2; l++) tab00[j][l]=0;


    /* data summary: count N_{ijk} */
    count=count00=0;
    for(k=1; k<=N_data; k++){

        for(num00=0,s=numparent; s>=1; s--){
            tep=1;
            for(j=1; j<s; j++) tep*=10;
            num00+=data[k][parlist[s]]*tep;
           }
        num=num00*10+data[k][x[i]];

        if(count00==0){
            count00++;
            tab00[count00][1]=num00;
            tab00[count00][2]+=1;
           }
          else{
            j=1;
            while(tab00[j][1]<num00 && j<=count00) j++;

            if(tab00[j][1]==num00) tab00[j][2]+=1;
               else{
                for(l=count00; l>=j; l--){
                    tab00[l+1][1]=tab00[l][1];
                    tab00[l+1][2]=tab00[l][2];
                   }
                tab00[j][1]=num00; tab00[j][2]+=1;
                count00++;
               }
           }


        if(count==0){
            count++;
            tab[count][1]=num;
            tab[count][2]+=1;
           }
          else{
            j=1;
            while(tab[j][1]<num && j<=count) j++;

            if(tab[j][1]==num) tab[j][2]+=1;
               else{
                for(l=count; l>=j; l--){
                    tab[l+1][1]=tab[l][1];
                    tab[l+1][2]=tab[l][2];
                   }
                tab[j][1]=num; tab[j][2]=1;
                count++;
               }
           }
       } /* end data summary */

     alphaijk=prior_alpha/parstate/state[x[i]];
     alphaik=prior_alpha/parstate;

     fvalue[i]-=count*gammln(alphaijk);
     for(k=1; k<=count; k++) fvalue[i]+=gammln(tab[k][2]+alphaijk);

     fvalue[i]+=count00*gammln(alphaik);
     for(k=1; k<=count00; k++) fvalue[i]-=gammln(tab00[k][2]+alphaik);
     fvalue[i]*=-1.0;

     free_imatrix(tab,1,(parstate+1)*state[x[i]],1,2);
     free_imatrix(tab00,1,parstate+1,1,2);
   }

ABC:
  for(sum=0.0,m=1; m<=node_num; m++){
      sum+=fvalue[m];
     }

  free( (char*) parlist ) ;

  return sum;
}


