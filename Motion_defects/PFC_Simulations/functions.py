import numpy as np
import matplotlib.pyplot as plt
import random

def En_NL_u(u,theta,kx,ky,GRADu_x,GRADu_y,gamma,Bx,By):
    nx=np.cos(theta);ny=np.sin(theta);
    T1=-ddFF(GRADu_x*(-ny*Bx+nx*By),kx,ky,1,0)-ddFF(GRADu_y*(-ny*Bx+nx*By),kx,ky,0,1)
    T2=+ddFF(-ny*(-GRADu_x*ny+GRADu_y*nx),kx,ky,1,0)+ddFF(nx*(-GRADu_x*ny+GRADu_y*nx),kx,ky,0,1)
    return gamma*T2+T1

def En_NL_theta(u,theta,GRADu_x,GRADu_y,gamma,Bx,By):
    nx=np.cos(theta);ny=np.sin(theta);
    return gamma*((nx*GRADu_x+ny*GRADu_y))*(-ny*GRADu_x+nx*GRADu_y)-0.5*(nx*Bx+ny*By)*(GRADu_x**2+GRADu_y**2)

def SH2D_FT_lin(u,theta,NL_u,NL_theta,kx2,ky2,dt,eps,alpha,sigma):
    L_u=eps-alpha*(1.0-2*(kx2+ky2)+(kx2+ky2)**2)
    res_u = (u+dt*NL_u)/(1-dt*L_u)
    scal_t=1
    L_theta=-sigma*(kx2+ky2)
    res_theta = (theta+dt*scal_t*NL_theta)/(1-dt*scal_t*L_theta)
    return res_u,res_theta

def u_dislo_single(xmi,ymi,x0,y0,b):
    nu=1./3.
    ux=np.zeros(xmi.shape)
    uy=np.zeros(xmi.shape)
    xm=xmi-x0;ym=ymi-y0;
    ux=ux+(b/(2*np.pi))*(np.arctan2(ym,xm)+xm*ym/(2*(1-nu)*(xm**2+ym**2+0.001)))
    uy=uy-(b/(2*np.pi))*(((1-2*nu)/(4-4*nu))*np.log(xm**2+ym**2)+(xm**2-ym**2)/(4*(1-nu)*(xm**2+ym**2+0.001)))
    return ux,uy

def u_dipole(xmi,ymi,L,n=8):
    uuxt,uuyt=u_dislo_single(xmi,ymi,0,L/n,(2*np.pi))
    uux,uuy=u_dislo_single(xmi,ymi,0,-L/n,(-2*np.pi))
    uux=uux+uuxt
    uuy=uuy+uuyt
    return uux,uuy

def generate_random_points(count,min_distance,L):
    points = np.zeros((count**2,2))
    dx=L/count
    i=0
    for ix in range(count):
        for iy in range(count):
            points[i,0]=dx*ix+np.random.uniform(-1,1)*dx/3-L/2
            points[i,1]=dx*iy+np.random.uniform(-1,1)*dx/3-L/2
            if(points[i,0]>L/2):
                points[i,0]=points[i,0]-L
            if(points[i,1]>L/2):
                points[i,1]=points[i,1]-L
            if(points[i,0]<-L/2):
                points[i,0]=points[i,0]+L
            if(points[i,1]<-L/2):
                points[i,1]=points[i,1]+L
            i=i+1
    return np.array(points),i-1

def u_ndefect(xmi,ymi,ndefects,L,dist):
    np.random.seed(42)
    uux=np.zeros(xmi.shape)
    uuy=np.zeros(xmi.shape)
    points,ndd=generate_random_points(ndefects,dist,L)
    for i in range(0,ndd):
        uuxt,uuyt=u_dislo_single(xmi,ymi,points[i,0],points[i,1],((-1)**int(np.random.rand()>0.5))*(2*np.pi))
        uux=uux+uuxt
        uuy=uuy+uuyt
    return uux,uuy

def ddFF(u,kx,ky,ex,ey):
    I = complex(0,1)
    uhat=np.fft.fft2(u)
    derder=np.real(np.fft.ifft2( ( (I*ky)**ex ) * ( (I*kx)**ey )*uhat) )
    return derder

def ddFF2(u,kx,ky):
    I = complex(0,1)
    uhat=np.fft.fft2(u)
    derder=np.real(np.fft.ifft2( ( ( ((I*kx)**2) + ((I*ky)**2) )**2) * uhat ) )
    return derder

def energyFT(u,theta,kx,ky,eps,alpha,gamma,sigma,Bx,By,GRADu_x,GRADu_y,GRADth_x,GRADth_y,dx):
    lap=ddFF(u,kx,ky,2,0)+ddFF(u,kx,ky,0,2)
    lap2=ddFF2(u,kx,ky)
    nx=np.cos(theta);ny=np.sin(theta)
    res1=-1.0*eps*(u**2)/2+alpha*((u**2)/2+0.25*(u**4)+u*(2*lap+lap2))
    res2=-0.5*(GRADu_x**2+GRADu_y**2)*(-ny*Bx+nx*By)
    res3=(gamma/2)*(-ny*GRADu_x+nx*GRADu_y)**2
    res4=sigma*(ddFF(theta,kx,ky,2,0)+ddFF(theta,kx,ky,0,2))
    return np.array(res1+res2+res3+res4)
