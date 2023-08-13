/* Smooth estimate of instantaneous phase. */
/*
  Copyright (C) 2011 Jilin University
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <math.h>

#include <rsf.h>

int main (int argc, char* argv[])
{
    int nh, n1,n2,i1,i2, i, n12, niter, dim, n[SF_MAX_DIM], rect[SF_MAX_DIM];
    float *trace, *hilb, *num, *den, *phase, a,b,c, mean;
    char key[6];
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    
    dim = sf_filedims (in,n);
    n1 = n[0];
    n12 = 1;
    for (i=0; i < dim; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (!sf_getint(key,rect+i)) rect[i]=1;
	/*( rect#=(1,1,...) smoothing radius on #-th axis )*/ 
	n12 *= n[i];
    }
    n2 = n12/n1;

    trace = sf_floatalloc(n1);
    hilb = sf_floatalloc(n1);

    num = sf_floatalloc(n12);
    den = sf_floatalloc(n12);
	
    phase = sf_floatalloc(n12);

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getint("order",&nh)) nh=10;
    /* Hilbert transformer order */
    if (!sf_getfloat("ref",&c)) c=1.;
    /* Hilbert transformer reference (0.5 < ref <= 1) */

    sf_hilbert_init(n1, nh, c);

    mean=0.;
    for (i=i2=0; i2 < n2; i2++) {
	sf_floatread(trace,n1,in);
	sf_hilbert(trace,hilb);


	for (i1=0; i1 < nh; i1++, i++) {
	    num[i] = 0.;
	    den[i] = 0.;
	}	
	
	for (i1=nh; i1 < n1-nh; i1++, i++) {
	    a = trace[i1];
	    b = hilb[i1];
	    num[i] = b;
	    den[i] = a;
		
	    mean += den[i]*den[i] ;
	}
	
	for (i1=n1-nh; i1 < n1; i1++, i++) {
	    num[i] = 0.;
	    den[i] = 0.;
	}
	
    } /* i2 */
    
    
    mean = sqrtf(n12/mean);
    
    for (i=0; i < n12; i++) {
	num[i] *= mean;
	den[i] *= mean;
    }
    
    sf_divn_init(dim, n12, n, rect, niter, true);
    sf_divn (num, den, phase);
    
    for (i=0; i < n12; i++) {
	
	phase[i] = atanf(phase[i])*90./SF_PI;
    }
    
    sf_floatwrite(phase,n12,out);
    
    exit(0);
}


/* 	$Id$	 */
