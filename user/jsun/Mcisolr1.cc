// Complex lowrank decomposition for 2-D isotropic wave propagation. 
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
#include <time.h>
#include <rsf.hh>

#include "vecmatop.hh"
#include "serialize.hh"

using namespace std;

static std::valarray<float> vs;
static std::valarray<float> ks;
static float dt;

int sample(vector<int>& rs, vector<int>& cs, CpxNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
    setvalue(res,cpx(0.0f,0.0f));
    for(int a=0; a<nr; a++) {
	for(int b=0; b<nc; b++) {
	    float phase = vs[rs[a]]*ks[cs[b]]*dt; 
	    res(a,b) = cpx(cos(phase),sin(phase)); 
//	    sf_warning("real=%g, imag=%g", real(res(a,b)),imag(res(a,b)));
	}
    }
    return 0;
}

int main(int argc, char** argv)
{   
    sf_init(argc,argv); // Initialize RSF

    iRSF par(0);
    int seed;

    par.get("seed",seed,time(NULL)); // seed for random number generator
    srand48(seed);

    float eps;
    par.get("eps",eps,1.e-4); // tolerance

    int npk;
    par.get("npk",npk,20); // maximum rank

    par.get("dt",dt); // time step

    iRSF vel;

    int nz;
    vel.get("n1",nz);
    int m = nz;
    std::valarray<float> vels(m);
    vel >> vels;
    vs.resize(m);
    vs = vels;
    
    iRSF fft("fft");

    int nkz;
    fft.get("n1",nkz);

    float dkz;
    fft.get("d1",dkz);
    
    float kz0;
    fft.get("o1",kz0);

    float kz;

    int n = nkz;
    std::valarray<float> k(n);
    for (int iz=0; iz < nkz; iz++) {
	kz = kz0+iz*dkz;
	k[iz] = 2*SF_PI*abs(kz);
    }
    ks.resize(n);
    ks = k;

    vector<int> lidx, ridx;
    CpxNumMat mid;

    iC( lowrank(m,n,sample,eps,npk,lidx,ridx,mid) );

    int n2=mid.n();
    int m2=mid.m();

    vector<int> midx(m), nidx(n);
    for (int k=0; k < m; k++) 
	midx[k] = k;
    for (int k=0; k < n; k++) 
	nidx[k] = k;    

    CpxNumMat lmat(m,m2);
    iC ( sample(midx,lidx,lmat) );

    CpxNumMat lmat2(m,n2);
    iC( zgemm(1.0, lmat, mid, 0.0, lmat2) );

    cpx *ldat = lmat2.data();
    std::valarray<sf_complex> ldata(m*n2);
    for (int k=0; k < m*n2; k++) {
	ldata[k] = sf_cmplx(real(ldat[k]),imag(ldat[k]));
//	sf_warning("real of ldat=%g, imag of ldat=%g", real(ldat[k]),imag(ldat[k]));
    }

    oRSF left("left");
    left.type(SF_COMPLEX);
    left.put("n1",m);
    left.put("n2",n2);
    left << ldata;

    CpxNumMat rmat(n2,n);
    iC ( sample(ridx,nidx,rmat) );

    cpx *rdat = rmat.data();
    std::valarray<sf_complex> rdata(n2*n);    
    for (int k=0; k < n2*n; k++) {
	rdata[k] = sf_cmplx(real(rdat[k]),imag(rdat[k]));
//    	sf_warning("real of rdat=%g, imag of rdat=%g", real(rdat[k]),imag(rdat[k]));
    }
    oRSF right;
    right.type(SF_COMPLEX);
    right.put("n1",n2);
    right.put("n2",n);
    right << rdata;

    exit(0);
}
