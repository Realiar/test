/* 2-D two-components wavefield modeling using pseudo-pure mode P-wave equation in TTI media.

   Copyright (C) 2012 Tongji University, Shanghai, China 
   Authors: Jiubing Cheng, Wei Kang and Tengfei Wang
     
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

/* prepared head files by myself */
#include "_fd.h"
#include "_cjb.h"

/* head files aumatically produced from *.c */
#include "ricker.h"
#include "puthead.h"
#include "zero.h"
#include "fdcoef.h"
#include "clipsmthspec.h"

/*wavefield propagators */
#include "fwpttipseudop.h"

/*****************************************************************************************/
int main(int argc, char* argv[])
{
	int	i,j,im, jm,nx,nz,nxpad,nzpad,it,ii,jj;

        float   A, f0, t, t0, dx, dz, dt, dt2, div;
        int     mm, nvx, nvz, ns;
        int     isx, isz, isxm, iszm; /*source location */
        float   *coeff_1dx, *coeff_1dz, *coeff_2dx, *coeff_2dz; /* finite-difference coefficient */

        float **vp0, **vs0, **epsi, **del, **theta;         /* velocity model */
        float **p1, **p2, **p3, **q1, **q2, **q3;  /* wavefield array */

        float   fx,fz; 

        sf_init(argc,argv);

        sf_file Fo;
       
        /*  wavelet parameter for source definition */
        f0=30.0;                  
        t0=0.04;                  
        A=1.0;                  

        /* time samping paramter */
        if (!sf_getint("ns",&ns)) ns=301;
        if (!sf_getfloat("dt",&dt)) dt=0.001;
        sf_warning("ns=%d dt=%f",ns,dt);
        sf_warning("read velocity model parameters");

        /* setup I/O files */
        sf_file Fvp0, Fvs0, Feps, Fdel, Fthe;

        Fvp0 = sf_input ("in");  /* vp0 using standard input */
        Fvs0 = sf_input ("vs0");  /* vs0 */
        Feps = sf_input ("epsi");  /* epsi */
        Fdel = sf_input ("del");  /* delta */
        Fthe = sf_input ("the");  /* theta */

        /* Read/Write axes */
        sf_axis az, ax;
        az = sf_iaxa(Fvp0,1); nvz = sf_n(az); dz = sf_d(az)*1000.0;
        ax = sf_iaxa(Fvp0,2); nvx = sf_n(ax); dx = sf_d(ax)*1000.0;
        fx=sf_o(ax);
        fz=sf_o(az);

        /* source definition */
        isx=nvx/2;
        isz=nvz/2;
        //isz=nvz*2/5;

        /* wave modeling space */
	nx=nvx;
	nz=nvz;
        nxpad=nx+2*_m;
        nzpad=nz+2*_m;

        sf_warning("dx=%f dz=%f",dx,dz);

        sf_warning("nx=%d nz=%d nxpad=%d nzpad=%d", nx,nz,nxpad,nzpad);

	vp0=sf_floatalloc2(nz,nx);	
	vs0=sf_floatalloc2(nz,nx);	
	epsi=sf_floatalloc2(nz,nx);	
	del=sf_floatalloc2(nz,nx);	
	theta=sf_floatalloc2(nz,nx);	

        int nxz=nx*nz;
        mm=2*_m+1;

        dt2=dt*dt;
        isxm=isx+_m;  /* source's x location */
        iszm=isz+_m;  /* source's z-location */

        /* read velocity model */
        sf_floatread(vp0[0],nxz,Fvp0);
        sf_floatread(vs0[0],nxz,Fvs0);
        sf_floatread(epsi[0],nxz,Feps);
        sf_floatread(del[0],nxz,Fdel);
        sf_floatread(theta[0],nxz,Fthe);

        for(i=0;i<nx;i++)
        for(j=0;j<nz;j++)
           theta[i][j] *= SF_PI/180.0;

        /* setup I/O files */
        Fo = sf_output("out"); /* pseudo-pure P-wave x-component */
        puthead3(Fo, nz, nx, 1, dz/1000.0, dx/1000.0, dt, fz, fx, 0.0);

      /****************begin to calculate wavefield****************/
       coeff_2dx=sf_floatalloc(mm);
       coeff_2dz=sf_floatalloc(mm);
       coeff_1dx=sf_floatalloc(mm);
       coeff_1dz=sf_floatalloc(mm);

      coeff2d(coeff_2dx,dx);
      coeff2d(coeff_2dz,dz);
      coeff1d(coeff_1dx,dx);
      coeff1d(coeff_1dz,dz);

      p1=sf_floatalloc2(nzpad, nxpad);
      p2=sf_floatalloc2(nzpad, nxpad);
      p3=sf_floatalloc2(nzpad, nxpad);

      q1=sf_floatalloc2(nzpad, nxpad);
      q2=sf_floatalloc2(nzpad, nxpad);
      q3=sf_floatalloc2(nzpad, nxpad);

      zero2float(p1, nzpad, nxpad);
      zero2float(p2, nzpad, nxpad);
      zero2float(p3, nzpad, nxpad);
        
      zero2float(q1, nzpad, nxpad);
      zero2float(q2, nzpad, nxpad);
      zero2float(q3, nzpad, nxpad);
        
      sf_warning("==================================================");
      sf_warning("==  Propagation Using Pseudo-Pure P-Wave Eq.    ==");
      sf_warning("==================================================");
      for(it=0;it<ns;it++)
      {
		t=it*dt;

		p2[isxm][iszm]+=Ricker(t, f0, t0, A);
		q2[isxm][iszm]+=Ricker(t, f0, t0, A);

        /* fwpttipseudop: forward-propagating in TTI media with pseudo-pure P-wave equation */
		fwpttipseudop(dt2, p1, p2, p3, q1, q2, q3, coeff_2dx, coeff_2dz,
                      dx, dz, nx, nz, nxpad, nzpad, vp0, vs0, epsi, del, theta);

          /******* output wavefields: component and divergence *******/
	       if(it==ns-1)
		{
		   for(i=0;i<nx;i++)
                   {
                      im=i+_m;
		      for(j=0;j<nz;j++)
		      {
                          jm=j+_m;

		          div=p3[im][jm]+q3[im][jm];

	                  sf_floatwrite(&div,1,Fo);
		       }
                     }/* i loop*/
                        
                }/* (it+1)%ntstep==0 */

                /**************************************/
 	        for(i=0,ii=_m;i<nx;i++,ii++)
	        for(j=0,jj=_m;j<nz;j++,jj++)
		{
				p1[ii][jj]=p2[ii][jj];	
				p2[ii][jj]=p3[ii][jj];	

				q1[ii][jj]=q2[ii][jj];	
				q2[ii][jj]=q3[ii][jj];	
		}

		if(it%50==0)
			sf_warning("Pseudo: it= %d",it);
	}/* it loop */

        free(*p1);
        free(*p2);
        free(*p3);
        free(*q1);
        free(*q2);
        free(*q3);

        free(*vp0);
        free(*vs0);
        free(*epsi);
        free(*del);
        free(*theta);

	exit(0);
}
