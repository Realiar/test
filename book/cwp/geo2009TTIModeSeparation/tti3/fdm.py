from rsf.proj import *
import pplot

# default parameters
def param(par):
    if('lt' not in par):       par['lt']='t'
    if('lz' not in par):       par['lz']='z'
    if('lx' not in par):       par['lx']='x'

    if('ut' not in par):       par['ut']='s'
    if('uz' not in par):       par['uz']='m'
    if('ux' not in par):       par['ux']='m'

    if('nb' not in par):       par['nb']=0
    if('nbz' not in par):      par['nbz']=100
    if('nbx' not in par):      par['nbx']=100
    if('tz' not in par):       par['tz']=0.0035
    if('tx' not in par):       par['tx']=0.0035

    if('nbell' not in par):    par['nbell']=1

    if('snap' not in par):     par['snap']='y'
    if('jsnap' not in par):    par['jsnap']=100
    if('jdata' not in par):    par['jdata']=1

    if('ompchunk' not in par): par['ompchunk']=1
    if('ompnth' not in par):   par['ompnth']=0
    if('free' not in par):     par['free']='n'
   
    if('ot' not in par):       par['ot']=0.
    if('nt' not in par):       par['nt']=1
    if('dt' not in par):       par['dt']=1.

    if('tmin' not in par):     par['tmin']=par['ot']
    if('tmax' not in par):     par['tmax']=par['ot'] + (par['nt']-1) * par['dt']
    if('xmin' not in par):     par['xmin']=par['ox']
    if('xmax' not in par):     par['xmax']=par['ox'] + (par['nx']-1) * par['dx']
    if('ymin' not in par):     par['ymin']=par['oy']
    if('ymax' not in par):     par['ymax']=par['oy'] + (par['ny']-1) * par['dy']
    if('zmin' not in par):     par['zmin']=par['oz']
    if('zmax' not in par):     par['zmax']=par['oz'] + (par['nz']-1) * par['dz']

    if('ratio' not in par):    par['ratio']=1.0*(par['zmax']-par['zmin'])/(par['xmax']-par['xmin'])

    if(par['ratio']>=1):
        par['height']=10
    else:
        par['height']=13*par['ratio']

    dx=par['xmax']-par['xmin'];
    dy=par['ymax']-par['ymin'];
    dz=par['zmax']-par['zmin'];
    yxratio=dx/(dx+dy);
    yzratio=dz/(dz+dy);
    
    par['ratio3d']=(dz+dy)/(dx+dy);
    par['pointz']=yzratio;
    par['pointx']=yxratio;
    if(par['ratio3d']>=1):
        par['height3d']=10
    else:
        par['height3d']=13*par['ratio3d']
            
    if(('scalebar') not in par): par['scalebar']='n'
    if(('labelattr') not in par): par['labelattr']=' labelsz=6 labelfat=3 titlesz=12 titlefat=3 '
    # parallel2=n
    if(('nqz') not in par): par['nqz']=par['nz']
    if(('oqz') not in par): par['oqz']=par['oz']
    if(('dqz') not in par): par['dqz']=par['dz']

    if(('nqx') not in par): par['nqx']=par['nx']
    if(('oqx') not in par): par['oqx']=par['ox']
    if(('dqx') not in par): par['dqx']=par['dx']
    
    par['labelattr']=' '+par['labelattr']+' '

# ------------------------------------------------------------
# plotting functions
def cgrey(custom,par):
    return '''
    grey labelrot=n wantaxis=y title=""
    pclip=100
    min1=%g max1=%g label1=%s unit1=%s
    min2=%g max2=%g label2=%s unit2=%s
    screenratio=%g screenht=%g wantscalebar=%s
    %s
    ''' % (par['zmin'],par['zmax'],par['lz'],par['uz'],
           par['xmin'],par['xmax'],par['lx'],par['ux'],
           par['ratio'],par['height'],par['scalebar'],
           par['labelattr']+' '+custom)

def cgrey3d(custom,par):
    return '''
    window min1=%g max1=%g min2=%g max2=%g min3=%g max3=%g |
    byte gainpanel=a pclip=100 %s |
    grey3 title="" framelabel=n
    label1=%s unit1=%s
    label2=%s unit2=%s
    label3=%s unit3=%s
    frame1=%d frame2=%d frame3=%d
    flat=y screenratio=%g screenht=%g point1=%g point2=%g
    %s
    ''' % (par['zmin'],par['zmax'],
           par['xmin'],par['xmax'],
           par['ymin'],par['ymax'],
           custom,
           par['lz'],par['uz'],
           par['lx'],par['ux'],
           par['ly'],par['uy'],
           par['nz']/2,par['nx']/2,par['ny']/2,
           par['ratio3d'],par['height3d'],par['pointz'],par['pointx'],
           par['labelattr']+' '+custom)



def dgrey3d(custom,par):
    return '''
    window min1=%g max1=%g min2=%g max2=%g min3=%g max3=%g |
    byte gainpanel=a pclip=100 %s |
    grey3 title="" framelabel=n
    label1=%s unit1=%s
    label2=%s unit2=%s
    label3=%s unit3=%s
    frame1=%d frame2=%d frame3=%d
    flat=y screenratio=%g screenht=%g point1=%g point2=%g
    %s
    ''' % (par['tmin'],par['tmax'],
           par['xmin'],par['xmax'],
           par['ymin'],par['ymax'],
           custom,
           par['lt'],par['ut'],
           par['lx'],par['ux'],
           par['ly'],par['uy'],
           par['nt']/2,par['nx']/2,par['ny']/2,
           par['ratio3d'],par['height3d'],par['pointz'],par['pointx'],
           par['labelattr']+' '+custom)

def center3d(x,y,z,par):
    return '''
    frame1=%d frame2=%d frame3=%d
    ''' % ( (float(z)-par['oz'])/par['dz']+1,
            (float(x)-par['ox'])/par['dx']+1,
            (float(y)-par['oy'])/par['dy']+1 )

def wgrey(custom,par):
    return '''
    window min1=%g max1=%g min2=%g max2=%g |
    grey labelrot=n wantaxis=y title=""
    pclip=100 gainpanel=a
    label1=%s unit1=%s
    label2=%s unit2=%s
    screenratio=%g screenht=%g wantscalebar=%s
    %s
    ''' % (par['zmin'],par['zmax'],
           par['xmin'],par['xmax'],
           par['lz'],par['uz'],
           par['lx'],par['ux'],
           par['ratio'],par['height'],par['scalebar'],
           par['labelattr']+' '+custom)

def cgraph(custom,par):
    return '''
    graph labelrot=n wantaxis=n title="" yreverse=y
    min2=%g max2=%g label2=%s unit2=%s
    min1=%g max1=%g label1=%s unit1=%s
    screenratio=%g screenht=%g wantscalebar=%s
    %s
    ''' % (par['zmin'],par['zmax'],par['lz'],par['uz'],
           par['xmin'],par['xmax'],par['lx'],par['ux'],
           par['ratio'],par['height'],par['scalebar'],
           par['labelattr']+' '+custom)

def ccont(custom,par):
    return '''
    contour labelrot=n wantaxis=n title=""
    min1=%g max1=%g label1=%s unit1=%s
    min2=%g max2=%g label2=%s unit2=%s
    screenratio=%g screenht=%g wantscalebar=%s
    plotcol=2 plotfat=3
    %s
    ''' % (par['zmin'],par['zmax'],par['lz'],par['uz'],
           par['xmin'],par['xmax'],par['lx'],par['ux'],
           par['ratio'],par['height'],par['scalebar'],
           par['labelattr']+' '+custom)

def dgrey(custom,par):
    return '''
    grey labelrot=n wantaxis=y title=""
    pclip=100
    min1=%g max1=%g label1=%s unit1=%s
    min2=%g max2=%g label2=%s unit2=%s
    %s
    ''' % (par['tmin'],par['tmax'],par['lt'],par['ut'],
           par['xmin'],par['xmax'],par['lx'],par['ux'],
           par['labelattr']+' '+custom)

def egrey(custom,par):
    return '''
    grey labelrot=n wantaxis=y title=""
    pclip=100
    min2=%g max2=%g label2=%s unit2=%s
    min1=%g max1=%g label1=%s unit1=%s
    %s
    ''' % (par['tmin'],par['tmax'],par['lt'],par['ut'],
           par['zmin'],par['zmax'],par['lz'],par['uz'],
           par['labelattr']+' '+custom)

# ------------------------------------------------------------
# plot wavelet
def waveplot(custom,par):
    return '''
    graph min2=-1 max2=+1 title=""
    plotfat=5 plotcol=2
    label1=%s unit1=%s label2="" unit2=""
    %s
    ''' % (par['lt'],par['ut'],
           par['labelattr']+' '+custom)

# ------------------------------------------------------------
# create wavelet
def wavelet(wav,frequency,par):
    par['frequency'] = frequency
    
    Flow(wav,None,
         '''
         spike nsp=1 mag=1 n1=%(nt)d d1=%(dt)g o1=%(ot)g k1=%(kt)d |
         pad end1=%(nt)d |
         ricker1 frequency=%(frequency)g |
         window n1=%(nt)d |
         scale axis=123 |
         put label1=t 
         ''' % par)    

# ------------------------------------------------------------
def horizontal(cc,coord,par):
    Flow(cc+'_',None,'math n1=%(nx)d d1=%(dx)g o1=%(ox)g output=0' % par)
    Flow(cc+'_z',cc+'_','math output="%g" ' % coord)
    Flow(cc+'_x',cc+'_','math output="x1" ')
    Flow(cc,[cc+'_x',cc+'_z'],
         '''
         cat axis=2 space=n
         ${SOURCES[0]} ${SOURCES[1]} | transp
         ''', stdin=0)

def horizontal3d(cc,coord,par):
    Flow(cc+'_',None,
         'math n1=%(nx)d d1=%(dx)g o1=%(ox)g n2=%(ny)d d2=%(dy)g o2=%(oy)g output=0' % par)
    Flow(cc+'_z',cc+'_','math output="%g" | put n1=%d n2=1' % (coord,par['nx']*par['ny']) )
    Flow(cc+'_x',cc+'_','math output="x1" | put n1=%d n2=1' % (      par['nx']*par['ny']) )
    Flow(cc+'_y',cc+'_','math output="x2" | put n1=%d n2=1' % (      par['nx']*par['ny']) )
    Flow(cc,[cc+'_x',cc+'_y',cc+'_z'],
         '''
         cat axis=2 space=n
         ${SOURCES[0:3]} | transp
         ''', stdin=0)

def vertical(cc,coord,par):
    Flow(cc+'_',None,'math n1=%(nz)d d1=%(dz)g o1=%(oz)g output=0' % par)
    Flow(cc+'_x',cc+'_','math output="%g" '% coord)
    Flow(cc+'_z',cc+'_','math output="x1" ')
    Flow(cc,[cc+'_x',cc+'_z'],
         '''
         cat axis=2 space=n
         ${SOURCES[0]} ${SOURCES[1]} | transp
         ''', stdin=0)

def point(cc,xcoord,zcoord,par):
    Flow(cc+'_',None,'math n1=1 d1=1 o1=0 output=0' % par)
    Flow(cc+'_z',cc+'_','math output="%g"' % zcoord)
    Flow(cc+'_x',cc+'_','math output="%g"' % xcoord)
    Flow(cc,[cc+'_x',cc+'_z'],
         '''
         cat axis=2 space=n
         ${SOURCES[0]} ${SOURCES[1]} | transp
         ''', stdin=0)
def point3d(cc,xcoord,ycoord,zcoord,par):
    Flow(cc+'_',None,'math n1=1 d1=1 o1=0 output=0' % par)
    Flow(cc+'_x',cc+'_','math output="%g"' % xcoord)
    Flow(cc+'_y',cc+'_','math output="%g"' % ycoord)
    Flow(cc+'_z',cc+'_','math output="%g"' % zcoord)

    Flow(cc,[cc+'_x',cc+'_y',cc+'_z'],
         '''
         cat axis=2 space=n
         ${SOURCES[0:3]} | transp
         ''', stdin=0)
    

def point3(cc,xcoord,zcoord,magn,par):
    Flow(cc+'_',None,'math n1=1 d1=1 o1=0 output=0' % par)
    Flow(cc+'_z',cc+'_','math output="%g"' % zcoord)
    Flow(cc+'_x',cc+'_','math output="%g"' % xcoord)
    Flow(cc+'_r',cc+'_','math output="%g"' % magn)
    Flow(cc,[cc+'_x',cc+'_z',cc+'_r'],
         '''
         cat axis=2 space=n
         ${SOURCES[0]} ${SOURCES[1]} ${SOURCES[2]} | transp
         ''', stdin=0)

def circle(cc,xcenter,zcenter,radius,sampling,par):
    Flow(cc+'_x',None,
         'math n1=%d d1=%g o1=%g output="%g+%g*cos(x1/180*3.14)"'
         % (sampling,360./sampling,0.,xcenter,radius) )
    Flow(cc+'_z',None,
         'math n1=%d d1=%g o1=%g output="%g+%g*sin(x1/180*3.14)"'
         % (sampling,360./sampling,0.,zcenter,radius) )
    Flow(cc,[cc+'_x',cc+'_z'],
         '''
         cat axis=2 space=n
         ${SOURCES[0]} ${SOURCES[1]} | transp
         ''', stdin=0)

def boxarray(cc,nz,oz,dz,nx,ox,dx,par):
    Flow(cc+'_',None,
         '''
         math output=1
         n1=%d d1=%g o1=%g
         n2=%d d2=%g o2=%g
         ''' % (nz,dz,oz,nx,dx,ox) )
    Flow(cc+'_z',cc+'_','math output="x1" | put n1=%d n2=1' % (nz*nx))
    Flow(cc+'_x',cc+'_','math output="x2" | put n1=%d n2=1' % (nz*nx))
    Flow(cc,[cc+'_x',cc+'_z'],
         '''
         cat axis=2 space=n
         ${SOURCES[0]} ${SOURCES[1]} | transp
         ''', stdin=0)

def ssplot(custom,par):
    return '''
    window n1=2 |
    dd type=complex |
    ''' + cgraph('symbol=o plotcol=6 plotfat=10 wantaxis=n %s' % custom,par)

def rrplot(custom,par):
    return '''
    window n1=2 |
    dd type=complex |
    ''' + cgraph('symbol=. plotcol=3 plotfat=5 wantaxis=n %s' % custom,par)

def qqplot(custom,par):
    return '''
    window n1=2 |
    dd type=complex |
    ''' + cgraph('symbol=. plotcol=1 plotfat=3 wantaxis=n %s' % custom,par)

def qqwin(par):
    return '''
    nqz=%(nqz)d oqz=%(oqz)g dqz=%(dqz)g
    nqx=%(nqx)d oqx=%(oqx)g dqx=%(dqx)g
    ''' % par

# ------------------------------------------------------------
# rays plot
def rayplot(hwt,j1ray,j2ray,j1wft,j2wft,custom,par):

    Plot(hwt+'ray',hwt,'window squeeze=n j1=%d j2=%d f2=%d | transp |' %(j1ray,j2ray,j2wft)
         + cgraph('plotcol=1 wantaxis=n '+custom,par))
    Plot(hwt+'wft',hwt,'window j1=%d j2=%d f2=%d |'          %(j1wft,j2wft,j2wft)
         + cgraph('plotcol=2 squeeze=n wantaxis=n symbol=. '+custom,par))

    Plot  (hwt,[hwt+'ray',hwt+'wft'],'Overlay')
  
# ------------------------------------------------------------
# execute acoustic finite-differences modeling
def amodel(data,wfld,  wavl,velo,dens,sou,rec,custom,par):
    par['fdcustom'] = custom
    
    Flow( [data,wfld],[wavl,velo,dens,sou,rec],
          '''
          afmod ompchunk=%(ompchunk)d
          verb=y abc=y free=n dens=y
          snap=%(snap)s jsnap=%(jsnap)d
          nbz=%(nbz)d tz=%(tz)g
          nbx=%(nbx)d tx=%(tx)g
          vel=${SOURCES[1]}
          den=${SOURCES[2]}
          sou=${SOURCES[3]}
          rec=${SOURCES[4]}
          wfl=${TARGETS[1]}
          %(fdcustom)s
          ''' % par)

def lmodel(data,wfld,ldata,lwfld,  wavl,velo,refl,sou,rec,custom,par):
    par['fdcustom'] = custom

    Flow([ data,wfld,ldata,lwfld],[wavl,velo,refl,sou,rec],
         '''
         aborn ompchunk=%(ompchunk)d
         verb=y abc=y free=n
         snap=%(snap)s jsnap=%(jsnap)d
         nbz=%(nbz)d tz=%(tz)g
         nbx=%(nbx)d tx=%(tx)g
         vel=${SOURCES[1]}
         ref=${SOURCES[2]}
         sou=${SOURCES[3]}
         rec=${SOURCES[4]}
         wfl=${TARGETS[1]}
         lid=${TARGETS[2]}
         liw=${TARGETS[3]}
         %(fdcustom)s
         ''' % par)

# ------------------------------------------------------------
# isotropic stiffness tensor
def isotropic(cc,vp,vs,ro,par):
    
    # Lame parameters
    # lambda = ro * (vp^2 - 2 vs^2)
    Flow(cc+'lambda',[ro,vp,vs],
         '''
         math output="ro*(vp*vp-2*vs*vs)"
         ro=${SOURCES[0]}
         vp=${SOURCES[1]}
         vs=${SOURCES[2]}     
         ''')
    
    #     mu = ro *           vs^2
    Flow(cc+'mu',[ro,vs],
         '''
         math output="ro*vs*vs"
         ro=${SOURCES[0]}
         vs=${SOURCES[1]}
         ''')
    
    # c11 = lambda + 2 mu
    # c33 = lambda + 2 mu
    # c13 = lambda
    # c44 = mu
    Flow(cc+'11',[cc+'lambda',cc+'mu'],
         'math l=${SOURCES[0]} m=${SOURCES[1]} output="l+2*m"')
    Flow(cc+'13',[cc+'lambda',cc+'mu'],
         'math l=${SOURCES[0]} m=${SOURCES[1]} output="l"')
    Flow(cc+'33',[cc+'lambda',cc+'mu'],
         'math l=${SOURCES[0]} m=${SOURCES[1]} output="l+2*m"')
    Flow(cc+'44',[cc+'lambda',cc+'mu'],
         'math l=${SOURCES[0]} m=${SOURCES[1]} output="m"')
    Flow(cc,[cc+'11',cc+'13',cc+'33',cc+'44'],
         'cat axis=3 space=n ${SOURCES[1:4]}')

# ------------------------------------------------------------
# anisotropic stiffness tensor
def anisotropic(cc,vp,vs,ro,epsilon,delta,par):
    Flow(cc+'33',[vp,ro],
         '''
         math output="ro*vp*vp"
         vp=${SOURCES[0]}
         ro=${SOURCES[1]}
         ''')    
    Flow(cc+'44',[vs,ro],
         '''
         math output="ro*vs*vs"
         vs=${SOURCES[0]}
         ro=${SOURCES[1]}
         ''')
    Flow(cc+'11',[cc+'33',epsilon],
         '''
         math output="2*epsilon*c33+c33"
         c33=${SOURCES[0]}
         epsilon=${SOURCES[1]}
         ''')
    Flow(cc+'13',[cc+'33',cc+'44',delta],
         '''
         math output="sqrt(2*c33*(c33-c44)*delta+(c33-c44)*(c33-c44))-c44"
         c33=${SOURCES[0]}
         c44=${SOURCES[1]}
         delta=${SOURCES[2]}
         ''')
    
    Flow(cc,[cc+'11',cc+'13',cc+'33',cc+'44'],
         'cat axis=3 space=n ${SOURCES[1:4]}')
    
# ------------------------------------------------------------
# acoustic modeling
def awefd(odat,owfl,idat,velo,dens,sou,rec,custom,par):
    par['fdcustom'] = custom
    
    Flow([odat,owfl],[idat,velo,dens,sou,rec],
         '''
         awefd
         ompchunk=%(ompchunk)d ompnth=%(ompnth)d 
         verb=y free=n snap=%(snap)s jsnap=%(jsnap)d
         nb=%(nb)d
         vel=${SOURCES[1]}
         den=${SOURCES[2]}
         sou=${SOURCES[3]}
         rec=${SOURCES[4]}
         wfl=${TARGETS[1]}
         %(fdcustom)s
         ''' % par)

def awefd2d(odat,owfl,idat,velo,dens,sou,rec,custom,par):
    par['fdcustom'] = custom
    
    Flow([odat,owfl],[idat,velo,dens,sou,rec,'./AWEFD2d.x'],
         '''
         ./AWEFD2d.x
         ompchunk=%(ompchunk)d ompnth=%(ompnth)d 
         verb=y free=n snap=%(snap)s jsnap=%(jsnap)d
         nb=%(nb)d dabc=%(dabc)s
         vel=${SOURCES[1]}
         den=${SOURCES[2]}
         sou=${SOURCES[3]}
         rec=${SOURCES[4]}
         wfl=${TARGETS[1]}
         %(fdcustom)s
         ''' % par)
    
def awefd3d(odat,owfl,idat,velo,dens,sou,rec,custom,par):
    par['fdcustom'] = custom
    
    Flow([odat,owfl],[idat,velo,dens,sou,rec,'./AWEFD3d.x'],
         '''
         ./AWEFD3d.x
         ompchunk=%(ompchunk)d ompnth=%(ompnth)d 
         verb=y free=n snap=%(snap)s jsnap=%(jsnap)d
         nb=%(nb)d dabc=%(dabc)s
         vel=${SOURCES[1]}
         den=${SOURCES[2]}
         sou=${SOURCES[3]}
         rec=${SOURCES[4]}
         wfl=${TARGETS[1]}
         %(fdcustom)s
         ''' % par)
    
def awefd1(odat,owfl,idat,velo,dens,sou,rec,custom,par):
    awefd(odat,owfl,idat,velo,dens,sou,rec,custom+' expl=y ',par)
    

def lwefd(bdat,bwfl,sdat,swfl,idat,velo,dens,refl,sou,rec,custom,par):
    par['fdcustom'] = custom
    
    Flow([bdat,bwfl,sdat,swfl],[idat,velo,dens,refl,sou,rec],
         '''
         lwefd
         ompchunk=%(ompchunk)d ompnth=%(ompnth)d 
         verb=y free=n snap=%(snap)s jsnap=%(jsnap)d
         nb=%(nb)d
         vel=${SOURCES[1]}
         den=${SOURCES[2]}
         ref=${SOURCES[3]}
         sou=${SOURCES[4]}
         rec=${SOURCES[5]}
         wfl=${TARGETS[1]}
         lid=${TARGETS[2]}
         liw=${TARGETS[3]}
         %(fdcustom)s
         ''' % par)
def lwefd1(bdat,bwfl,sdat,swfl,idat,velo,dens,refl,sou,rec,custom,par):
    lwefd(bdat,bwfl,sdat,swfl,idat,velo,dens,refl,sou,rec,custom+' expl=y ',par)

# ------------------------------------------------------------
# elastic modeling
def ewefd(odat,owfl,idat,cccc,dens,sou,rec,custom,par):
    par['fdcustom'] = custom
    
    Flow( [odat,owfl],[idat,cccc,dens,sou,rec],
         '''
         ewefd
         ompchunk=%(ompchunk)d  ompnth=%(ompnth)d 
         verb=y free=n snap=%(snap)s jsnap=%(jsnap)d nb=%(nb)d nbell=%(nbell)d
         ccc=${SOURCES[1]}
         den=${SOURCES[2]}
         sou=${SOURCES[3]}
         rec=${SOURCES[4]}
         wfl=${TARGETS[1]}
         %(fdcustom)s
         ''' % par)

def ewefd2d(odat,owfl,idat,cccc,dens,sou,rec,custom,par):
    par['fdcustom'] = custom
    
    Flow( [odat,owfl],[idat,cccc,dens,sou,rec,'./EWEFD2d.x'],
         '''
         ./EWEFD2d.x
         ompchunk=%(ompchunk)d  ompnth=%(ompnth)d 
         verb=y free=n snap=%(snap)s jsnap=%(jsnap)d nbell=%(nbell)d
         nb=%(nb)d  dabc=%(dabc)s
         ccc=${SOURCES[1]}
         den=${SOURCES[2]}
         sou=${SOURCES[3]}
         rec=${SOURCES[4]}
         wfl=${TARGETS[1]}
         %(fdcustom)s
         ''' % par)

def ewefd3d(odat,owfl,idat,cccc,dens,sou,rec,custom,par):
    par['fdcustom'] = custom
    
    Flow( [odat,owfl],[idat,cccc,dens,sou,rec],
         '''
         ewefd3dtti
         ompchunk=%(ompchunk)d  ompnth=%(ompnth)d 
         verb=y free=n snap=%(snap)s jsnap=%(jsnap)d nbell=%(nbell)d
         nb=%(nb)d  dabc=%(dabc)s
         ccc=${SOURCES[1]}
         den=${SOURCES[2]}
         sou=${SOURCES[3]}
         rec=${SOURCES[4]}
         wfl=${TARGETS[1]}
         %(fdcustom)s
         ''' % par)

def ewefd3(odat,owfl,idat,cccc,dens,sou,rec,custom,par):
    par['fdcustom'] = custom
    
    Flow( [odat,owfl],[idat,cccc,dens,sou,rec,'Code/ewe3c.x'],
         '''
         Code/ewe3c.x         
         ompchunk=1 ompnth=0 verb=y free=n snap=y 
           jsnap=%(jsnap)d nb=%(nb)d ssou=n opot=n  nbell=5
         nb=%(nb)d  dabc=%(dabc)s
         ccc=${SOURCES[1]}
         den=${SOURCES[2]}
         sou=${SOURCES[3]}
         rec=${SOURCES[4]}
         wfl=${TARGETS[1]}
         %(fdcustom)s
         ''' % par)
# ------------------------------------------------------------
# heat modeling
def hdefd(dat,wfl,  wav,con,sou,rec,custom,par):
    par['fdcustom'] = custom
    
    Flow( [dat,wfl],[wav,con,sou,rec],
          '''
          hdefd
          verb=y free=n snap=%(snap)s jsnap=%(jsnap)d nb=%(nb)d
          con=${SOURCES[1]}
          sou=${SOURCES[2]}
          rec=${SOURCES[3]}
          wfl=${TARGETS[1]}
          %(fdcustom)s
          ''' % par)

# ------------------------------------------------------------    
# F-D modeling from arbitrary source/receiver geometry
def awe(odat,wfld,idat,velo,dens,sou,rec,custom,par):
    par['fdcustom'] = custom
    
    Flow( [odat,wfld],[idat,velo,dens,sou,rec],
          '''
          awe ompchunk=%(ompchunk)d
          verb=y abc=y free=n dens=y
          snap=%(snap)s jsnap=%(jsnap)d
          nbz=%(nbz)d tz=%(tz)g
          nbx=%(nbx)d tx=%(tx)g
          vel=${SOURCES[1]}
          den=${SOURCES[2]}
          sou=${SOURCES[3]}
          rec=${SOURCES[4]}
          wfl=${TARGETS[1]}
          %(fdcustom)s
          ''' % par)

# ------------------------------------------------------------
# shot-record reverse-time migration
def rtm(imag,sdat,rdat,velo,dens,sacq,racq,iacq,custom,par):

    swfl = imag+'_us' #   source wavefield
    rwfl = imag+'_ur' # receiver wavefield
    sout = imag+'_ds' #   source data (not the input sdat!)
    rout = imag+'_dr' # receiver data (not the input rdat!)

    # source wavefield (z,x,t)
    awefd(sout,swfl,sdat,velo,dens,sacq,iacq,custom,par)

    # receiver wavefield (z,x,t)
    tdat = imag+'_tds'
    tout = imag+'_tdr'

    Flow(tdat,rdat,'reverse which=2 opt=i verb=y')
    awefd(tout,rwfl,tdat,velo,dens,racq,iacq,custom,par)
    Flow(rout,tout,'reverse which=2 opt=i verb=y')

    # conventional (cross-correlation zero-lag) imaging condition
    Flow(imag,[sout,rout],'xcor uu=${SOURCES[1]} axis=2 verb=y nbuf=500')
    
#    corr = imag+'_cor'
#    Flow(corr,[sout,rout],'paradd mode=p ${SOURCES[1]} memsize=%d' %mem)
#    Flow(imag,corr,'stack axis=2')

# ------------------------------------------------------------
# exploding-reflector reverse-time migration
def zom(imag,data,rdat,velo,dens,sacq,racq,custom,par):

    rwfl = imag+'_ur' # receiver wavefield
    rout = imag+'_dr' # receiver data (not the input rdat)

    # receiver wavefield (z,x,t)
    tdat = imag+'_tds'
    twfl = imag+'_tur'

    Flow(tdat,rdat,'reverse which=2 opt=i verb=y')
    awefd(data,twfl,tdat,velo,dens,sacq,racq,custom+' jsnap=%d' % (par['nt']-1),par)

    Flow(imag,twfl,'window n3=1 f3=1')

# ------------------------------------------------------------
# wavefield-over-model plot
def wom(wom,wfld,velo,vmean,par):

    if(('wweight') not in par): par['wweight']=10
    if(('wclip') not in par):   par['wclip']=1.0

    chop = wfld+'_chop'
    Flow(chop,wfld,
         '''
         window
         min1=%(zmin)g max1=%(zmax)g
         min2=%(xmin)g max2=%(xmax)g |
         scale axis=123 |
         clip clip=%(wclip)g
         ''' % par)

    Flow(wom,[velo,chop],
         '''
         add add=-%g |
         scale axis=123 |
         spray axis=3 n=%d o=%g d=%g |
         math w=${SOURCES[1]} output="input+%g*w"
         ''' % (vmean,
                par['nt']/par['jsnap']+1,
                par['ot'],
                par['dt']*par['jsnap'],
                par['wweight']))

# ------------------------------------------------------------
# (elastic) wavefield-over-model
def wem(wom,wfld,velo,vmean,par):

    if(('wweight') not in par): par['wweight']=10
    if(('wclip') not in par):   par['wclip']=1.0

    Flow(velo+'-spray',
         velo,
         '''
         add add=-%g |
         scale axis=123 |
         spray axis=3 n=%d o=%g d=%g
         ''' % (vmean,
                par['nt']/par['jsnap'],
                par['ot'],
                par['dt']*par['jsnap']))

    for i in range(2):
        Flow(wom+'-%d' %i ,wfld,'window n3=1 f3=%d' %i)
        Flow(wom+'-%d' %i +'chop',
             wom+'-%d' %i,
             '''
             window
             min1=%(zmin)g max1=%(zmax)g
             min2=%(xmin)g max2=%(xmax)g |
             scale axis=123 |
             clip clip=%(wclip)g
             ''' % par)

        Flow(wom+'-%d' %i +'-wom',
             [velo+'-spray',wom+'-%d' %i + 'chop'],
             'math w=${SOURCES[1]} output="input+%g*w" | transp plane=34' % par['wweight'])

    Flow(wom,[wom+'-0-wom',wom+'-1-wom'],'cat axis=3 space=n ${SOURCES[1]}')
    
        
# ------------------------------------------------------------
# image-over-model plot
def iom(iom,imag,velo,vmean,par):

    if(('iweight') not in par): par['iweight']=10
    if(('iclip') not in par):   par['iclip']=1.0

    chop = imag+'_chop'
    Flow(chop,imag,
         '''
         window
         min1=%(zmin)g max1=%(zmax)g
         min2=%(xmin)g max2=%(xmax)g |
         scale axis=123 |
         clip clip=%(iclip)g
         ''' % par)

    Flow(iom,[velo,chop],
         '''
         add add=-%g |
         scale axis=123 |
         math w=${SOURCES[1]} output="input+%g*w"
         ''' % (vmean,par['iweight']))

# ------------------------------------------------------------
# wavefield snapshot plots
def wframe(frame,movie,index,custom,par):

    Flow([movie+'_plt',movie+'_bar'],movie,
         'byte bar=${TARGETS[1]} gainpanel=a pclip=100 %s' % custom)
    Plot  (frame,[movie+'_plt',movie+'_bar'],
           'window n3=1 f3=%d bar=${SOURCES[1]} |'% index
           + wgrey(custom,par))
    
# ------------------------------------------------------------
# elastic wavefield movie
#def emovieold(wfld,custom,axis,par):
#
#    # loop over wavefield components
#    for i in range(2):
#        Flow(wfld+str(i+1),wfld,
#             '''
#             window n3=1 f3=%d |
#             window min1=%g max1=%g min2=%g max2=%g
#             ''' % (i,par['zmin'],par['zmax'],par['xmin'],par['xmax']))
#
#    # join component wavefields
#    Flow(wfld+'all',[wfld+'1',wfld+'2'],
#         'cat axis=%d space=n ${SOURCES[1:2]}' % axis)
#
#    if(axis==1):        
#        height=2*par['height']
#        if(height>10): height=10
#        ratio =2*par['ratio']
#    if(axis==2):
#        height=par['height']
#        if(height>10): height=10
#        ratio =0.5*par['ratio']
#    
#    Result(wfld,wfld+'all',
#           '''
#           grey title="" wantaxis=y screenratio=%f screenht=%f
#           gainpanel=a pclip=99 %s
#           %s
#           ''' % (ratio,height,par['labelattr'],custom) )

# ------------------------------------------------------------
# elastic wavefield movie frames
def eframe(frame,movie,index,custom,axis,par,xscale=0.75,yscale=0.75,shift=-8.25):

    Flow([movie+'-plt',movie+'-bar'],movie,
         'byte bar=${TARGETS[1]} gainpanel=a pclip=100 %s' % custom)

    for i in range(2):
        Plot(frame+'-'+str(i),movie+'-plt',
             'window n3=1 f3=%d n4=1 f4=%d |' % (i,index)
             + cgrey('',par))
#        Result(frame+'-'+str(i),movie+'-plt',
#             'window n3=1 f3=%d n4=1 f4=%d |' % (i,index)
#             + cgrey('',par))

    if(axis==1):
        pplot.p2x1(frame,frame+'-1',frame+'-0',yscale,xscale,shift)
    else:
        pplot.p1x2(frame,frame+'-0',frame+'-1',yscale,xscale,shift)

# ------------------------------------------------------------
# elastic wavefield movie
def emovie(movie,wfld,nframes,custom,axis,par,xscale=0.75,yscale=0.75,shift=-8.25):
    
    for iframe in range(nframes):
        tag = "-%02d" % iframe
        eframe(movie+tag,wfld,iframe,custom,axis,par,xscale,yscale,shift)
        
    allframes = map(lambda x: movie+'-%02d'  % x,range(nframes))
    Result(movie,allframes,'Movie')

# ------------------------------------------------------------
# elastic data
def edata(plot,data,custom,par):

    Flow([plot+'_plt',plot+'_bar'],data,
         'scale axis=123 | byte bar=${TARGETS[1]} gainpanel=a pclip=100 %s' % custom)
    
    for i in range(2):
        Plot(  plot+str(i+1),[plot+'_plt',plot+'_bar'],
               'window n2=1 f2=%d bar=${SOURCES[1]} | transp |' % i
               + dgrey('%s' % custom,par))
        Result(plot+str(i+1),[plot+'_plt',plot+'_bar'],
               'window n2=1 f2=%d bar=${SOURCES[1]} | transp |' % i
               + dgrey('%s' % custom,par))

# ------------------------------------------------------------
# elastic image
def eimage(plot,imag,custom,par):

    Flow([plot+'_plt',plot+'_bar'],imag,
         'scale axis=123 | byte bar=${TARGETS[1]} gainpanel=a pclip=100 %s ' % custom)        

    for i in range(4):        
        Plot  (plot+str(i+1),[plot+'_plt',plot+'_bar'],
               'window n3=1 f3=%d bar=${SOURCES[1]} |'% i
               + cgrey('%s'% custom,par))
        Result(plot+str(i+1),[plot+'_plt',plot+'_bar'],
               'window n3=1 f3=%d bar=${SOURCES[1]} |'% i
               + cgrey('%s'% custom,par))
        
# ------------------------------------------------------------
# plot elastic wavelet
def ewavelet(wavelet,custom,par):
    
    for i in range(2):
        Plot(wavelet+'-'+str(i+1),wavelet,
             'window n2=1 f2=%d | transp | window |'%i +
             waveplot('%d %s'% (i,custom) ,par))
    pplot.p1x2(wavelet,wavelet+'-1',wavelet+'-2',0.5,0.5,-11)

def ewavelet3d(wavelet,custom,par):
    
    for i in range(3):
        Plot(wavelet+'-'+str(i+1),wavelet,
             'window n2=1 f2=%d | transp | window |'%i +
             waveplot('%d %s'% (i,custom) ,par))
    pplot.p1x3(wavelet,wavelet+'-1',wavelet+'-2',wavelet+'-3',0.35,0.35,-11)


# ------------------------------------------------------------
# acoustic RTM
def artm(imag,sdat,rdat,velo,dens,sacq,racq,iacq,custom,par):
    
    swfl = imag+'_us' #   source wavefield
    rwfl = imag+'_ur' # receiver wavefield
    twfl = imag+'_ut' #     temp wavefield
    sout = imag+'_ds' #   source data (not the input sdat!)
    rout = imag+'_dr' # receiver data (not the input rdat!)
    
    iwindow = ' ' + \
              '''
              nqz=%(nqz)d oqz=%(oqz)g
              nqx=%(nqx)d oqx=%(oqx)g
              jsnap=%(jdata)d jdata=%(jdata)d
              ''' % par + ' '
    
    # source wavefield (z,x,t)
    awefd(sout,swfl,sdat,velo,dens,sacq,iacq,custom+iwindow,par)
    
    # receiver wavefield (z,x,t)
    tdat = imag+'_tds'
    Flow(tdat,rdat,'reverse which=2 opt=i verb=y')
    awefd(rout,twfl,tdat,velo,dens,racq,iacq,custom+iwindow,par)
    Flow(rwfl,twfl,'reverse which=4 opt=i verb=y')
    
    # conventional (cross-correlation zero-lag) imaging condition
    Flow(imag,[swfl,rwfl],'xcor uu=${SOURCES[1]} axis=3 verb=y nbuf=100')
    
# ------------------------------------------------------------
# elastic RTM
def ertm(imag,sdat,rdat,cccc,dens,sacq,racq,iacq,custom,par):
    
    swfl = imag+'_us' #   source wavefield
    rwfl = imag+'_ur' # receiver wavefield
    twfl = imag+'_ut' #     temp wavefield
    sout = imag+'_ds' #   source data (not the input sdat!)
    rout = imag+'_dr' # receiver data (not the input rdat!)
    
    iwindow = ' ' + \
              '''
              nqz=%(nqz)d oqz=%(oqz)g
              nqx=%(nqx)d oqx=%(oqx)g
              jsnap=%(jdata)d jdata=%(jdata)d
              ''' % par + ' '
    
    # source wavefield (z,x,t)
    ewefd2(sout,swfl,sdat,cccc,dens,sacq,iacq," ssou=y opot=n" + custom+iwindow,par)
    
    # receiver wavefield (z,x,t)
    tdat = imag+'_tds'
    Flow(tdat,rdat,'reverse which=4 opt=i verb=y')
    ewefd2(rout,twfl,tdat,cccc,dens,racq,iacq," ssou=n opot=n" + custom+iwindow,par)

    # ------------------------------------------------------------
    Flow(swfl+'1',swfl,'window n3=1 f3=0')
    Flow(swfl+'2',swfl,'window n3=1 f3=1')

    Flow(rwfl+'1',twfl,'window n3=1 f3=0 | reverse which=4 opt=i verb=y')
    Flow(rwfl+'2',twfl,'window n3=1 f3=1 | reverse which=4 opt=i verb=y')

    # conventional (cross-correlation zero-lag) imaging condition
    for     i in ('1','2'):
        for j in ('1','2'):
            Flow(imag+i+j,[swfl+i,rwfl+j],
                 'xcor uu=${SOURCES[1]} axis=3 verb=y nbuf=100')
            
    Flow(imag,[imag+'11',imag+'12',imag+'21',imag+'22'],
         'cat axis=3 space=n ${SOURCES[1:4]}')

