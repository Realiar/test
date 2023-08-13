//   Ajoint test of prop1Pa and prop1P, prop2Na and prop2N
//   Copyright (C) 2010 University of Texas at Austin
//  
//   This program is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 2 of the License, or
//   (at your option) any later version.
//  
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//  
//   You should have received a copy of the GNU General Public License
//   along with this program; if not, write to the Free Software
//   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include <rsf.hh>
#include "vecmatop.hh"

extern "C" {
#include "prop1Pa.h"
#include "prop1P.h"
#include "prop2Na.h"
#include "prop2N.h"
#include "cfft2nsps.h"
}

using namespace std;
using std::cerr;

static int nz, nx, nzx, nkz, nkx, nkzx, m2;
static float dz, dx, z0, x0, dkz, dkx, kz0, kx0;

//------------------------------------------------------------

int main(int argc, char** argv)
{   
    sf_init(argc,argv); // Initialize RSF

    iRSF fft, model("model"), data("data"), left("left"), right("right"); // Get input
    fft.get("n1",nkz);
    fft.get("d1",dkz);
    fft.get("o1",kz0);
    fft.get("n2",nkx);
    fft.get("d2",dkx);
    fft.get("o2",kx0);
    left.get("n2",m2);

    model.get("n1",nz);
    model.get("d1",dz);
    model.get("o1",z0);
    model.get("n2",nx);
    model.get("d2",dx);
    model.get("o2",x0);

    nzx = nz*nx; 
    nkzx = nkz*nkx;

    std::valarray<sf_complex> mod(nzx), dat(nzx), lt(nzx*m2), rt(m2*nkzx);
    model >> mod;
    data >> dat;
    left >> lt;
    right >> rt;

    sf_complex *cleft, *cright, *cmod, *cdat, *cpmod, *cpadat, *cnmod, *cnadat;
    cleft = sf_complexalloc(nzx*m2);    
    cright = sf_complexalloc(m2*nkzx);
    cmod = sf_complexalloc(nzx);
    cdat = sf_complexalloc(nzx);   
    cpmod = sf_complexalloc(nzx);   
    cpadat = sf_complexalloc(nzx);  
    cnmod = sf_complexalloc(nzx);   
    cnadat = sf_complexalloc(nzx);    
    
    for (int i=0; i<nzx*m2; i++) cleft[i] = lt[i];
    for (int i=0; i<m2*nkzx; i++) cright[i] = rt[i];
    for (int i=0; i<nzx; i++) cmod[i] = mod[i];
    for (int i=0; i<nzx; i++) cdat[i] = dat[i];


    // Check inner product <cpmod,cdat>?=<mod,cpadat>
    prop1P( cmod, cpmod, cleft, cright, nz, nx, nkzx, m2);  
    prop1Pa( cdat, cpadat, cleft, cright, nz, nx, nkzx, m2);  
      
    double aar1=0., aai1=0., bbr1=0., bbi1=0.;
    for (int i=0; i<nzx; i++) {
        aar1 += crealf(sf_cmul(cpmod[i],conjf(cdat[i])));
        aai1 += cimagf(sf_cmul(cpmod[i],conjf(cdat[i])));
        bbr1 += crealf(sf_cmul(cmod[i],conjf(cpadat[i])));
	bbi1 += cimagf(sf_cmul(cmod[i],conjf(cpadat[i])));
    }

    cerr<<"aa1= "<<aar1<<" "<<aai1<<endl;
    cerr<<"bb1= "<<bbr1<<" "<<bbi1<<endl;

    
    // Check inner product <cnmod,cdat>?=<mod,cnadat>
    prop2N( cmod, cnmod, cleft, cright, nz, nx, nkzx, m2);  
    prop2Na( cdat, cnadat, cleft, cright, nz, nx, nkzx, m2);  
      
    double aar2=0., aai2=0., bbr2=0., bbi2=0.;
    for (int i=0; i<nzx; i++) {
        aar2 += crealf(sf_cmul(cnmod[i],conjf(cdat[i])));
        aai2 += cimagf(sf_cmul(cnmod[i],conjf(cdat[i])));
        bbr2 += crealf(sf_cmul(cmod[i],conjf(cnadat[i])));
	bbi2 += cimagf(sf_cmul(cmod[i],conjf(cnadat[i])));
    }

    cerr<<"aa2= "<<aar2<<" "<<aai2<<endl;
    cerr<<"bb2= "<<bbr2<<" "<<bbi2<<endl;
    

    exit(0);
}
