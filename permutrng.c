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

#include <stdlib.h>
#include <math.h>

double genrand_real3(void) ;

void permutrng(int *ivec, int n) {

	int i ;
	int j ;
	double u ;
	int itmp ;

	for (i=1; i<=n-1; i++) {
		u = genrand_real3() ;
		j = i + (int) floor((n-i+1)*u) ;
        itmp = ivec[j] ;
        ivec[j] = ivec[i] ;
        ivec[i] = itmp ;
	}

}

