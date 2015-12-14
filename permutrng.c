/*
 * Copyrigtht 2014 Georgios Karagiannis
 *
 * This file is part of PISAA_BNLDD.
 *
 * PISAA_BNLDD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 2 of the License.
 *
 * PISAA_BNLDD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISAA_BNLDD.  If not, see <http://www.gnu.org/licenses/>.
*/

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

