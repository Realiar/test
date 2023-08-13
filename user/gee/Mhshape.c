/* Helical autoregressive shaping */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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
#include <rsf.h>

#include "regrid.h"

int main(int argc, char* argv[])
{
    int dim, na, nx, i, j, n[SF_MAX_DIM], m[SF_MAX_DIM], ir, nr;
    float *data, *tmp1, *tmp2, eps, r;
    char *lagfile;
    sf_filter pef;
    sf_file inp, out, fil, lag;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");
    fil = sf_input("filt");

    if (!sf_getfloat("eps",&eps)) eps=1.0f;
    /* regularization parameter */
    eps *= eps;

    if (!sf_getint("rect",&nr)) nr=1;
    /* shaping radius */

    /* input data, output model */
    dim = sf_filedims(inp,n);
    
    if (!sf_histint(fil,"n1",&na)) sf_error("No n1= in sfilt");

    pef = sf_allocatehelix(na);

    if (NULL == (lagfile = sf_histstring(fil,"lag")) &&
	NULL == (lagfile = sf_getstring("lag"))) 
	sf_error("Need lag=");
    lag = sf_input(lagfile);
    if (!sf_histints(lag,"n",m,dim)) sf_error("No n= in %s",lagfile);
    sf_intread(pef->lag,na,lag);
    regrid(dim,m,n,pef);
    sf_fileclose(lag);

    sf_helicon_init (pef);

    /* input data, output model */
    nx=1;
    for(j=0; j < dim; j++) {
	nx *= n[j];
    }

    data = sf_floatalloc(nx);
    tmp1 = sf_floatalloc(nx);
    tmp2 = sf_floatalloc(nx);

    sf_floatread(pef->flt,na,fil);
    sf_floatread(data,nx,inp);

    for (ir=1; ir < nr; ir++) {
	r = 2*sinf(SF_PI*ir/nr);
	r = -0.5*eps/(r*r);

	for (i=0; i < nx; i++) {
	    tmp1[i] = r*data[i];
	}
	    
	sf_helicon_lop(false,false,nx,nx,tmp1,tmp2);
	sf_helicon_lop(true,true,nx,nx,data,tmp2);
    }
	
    sf_floatwrite(data,nx,out);
    
    exit(0);
}

