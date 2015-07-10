# -*- coding: utf-8 -*-
"""
Created on Sun Jul  5 16:15:10 2015

@author: Walid Hassani 
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Jul  5 10:59:01 2015

@author: walid
"""

from numpy import *
from matplotlib.pyplot import *

N = 1000

m=0.468
g=9.81

dt=0.001 

l=0.225 

k=2.980*10**(-6)
b=1.140*10**(-7)

Im = 3.357*10**(-5)
u=1000

phi = zeros(N)
psi = zeros(N)
theta = zeros(N)

phi_dot = zeros(N)
psi_dot = zeros(N)
theta_dot = zeros(N)

phi_ddot = zeros(N)
psi_ddot = zeros(N)
theta_ddot = zeros(N)

x = zeros(N)
y = zeros(N)
z = zeros(N)

x_dot = zeros(N)
y_dot = zeros(N)
z_dot = zeros(N)

x_ddot = zeros(N)
y_ddot = zeros(N)
z_ddot = zeros(N)

ez = zeros(N)
ez_dot = zeros(N)

ephi = zeros(N)
ephi_dot = zeros(N)

etheta = zeros(N)
etheta_dot = zeros(N)

epsi = zeros(N)
epsi_dot = zeros(N)



phi[0] = 0
psi[0] = 0
theta[0] = 0

phi_dot[0] = 0
psi_dot[0] = 0
theta_dot[0] = 0

x[0] = 0
y[0] = 0
z[0] = 0

x_dot[0] = 0
y_dot[0] = 0
z_dot[0] = 0
    
Ixx = 4.856*10**(-3)
Iyy = 4.856*10**(-3)
Izz = 8.801*10**(-3)
 
tau = array([0,0,0])
tau_tilde = array([0,0,0])

I = array([[Ixx, 0, 0], [0, Iyy, 0], [0, 0, Izz]])

phi_dot[0] = 0; 
theta_dot[0] = 0; 
psi_dot[0] = 0; 
eta_dot = array([phi_dot[0], theta_dot[0], psi_dot[0]])

phi[0] = 0; 
theta[0] = 0; 
psi[0] = 0;
eta = array([phi[0], theta[0], psi[0]])

u=1000
w1=u*ones(N); 
w2=u*ones(N);
w3=u*ones(N);
w4=u*ones(N);

T=zeros(N)
tau_phi=zeros(N)
tau_theta=zeros(N)
tau_psi=zeros(N)


zd=20*sin(np.array((range(N))) * 4*pi / 180. -pi/2)+20
thetad=zeros(N)
phid=zeros(N)
psid=sin(np.array((range(N))) * 0.004*pi / 180. -pi/2) + sin(np.array((range(N))) * 1*pi / 180.)
thetad=0.01*ones(N)

for i in range(N-1):

    J = array([[Ixx, 0, -Ixx*sin(theta[i])],[0, Iyy*(cos(phi[i])**2)+ Izz*(sin(phi[i]))**2, (Iyy-Izz)*cos(phi[i])*sin(phi[i])*cos(theta[i])],[-Ixx*sin(theta[i]), (Iyy-Izz)*cos(phi[i])*sin(phi[i])*cos(theta[i]), Ixx*sin(theta[i])**2+Iyy*sin(phi[i])**2*cos(theta[i])**2+Izz*cos(phi[i])**2*cos(theta[i])**2]])
    
    c = array([[0, (Iyy-Izz)*(theta_dot[i]*cos(phi[i])*sin(phi[i])+psi_dot[i]*(sin(phi[i]))**2*cos(theta[i]))+(Izz-Iyy)*psi_dot[i]*(cos(phi[i]))**2*cos(theta[i])-Ixx*psi_dot[i]*cos(theta[i]), (Izz-Iyy)*psi_dot[i]*cos(phi[i])*sin(phi[i])*(cos(theta[i]))**2],
               [(Izz-Iyy)*(theta_dot[i]*cos(phi[i])*sin(phi[i])+psi_dot[i]*(sin(phi[i]))**2*cos(theta[i]))+(Iyy-Izz)*psi_dot[i]*(cos(phi[i]))**2*cos(theta[i])-Ixx*psi_dot[i]*cos(theta[i]), (Izz-Iyy)*phi_dot[i]*cos(phi[i])*sin(phi[i]), -Ixx*psi_dot[i]*sin(theta[i])*cos(theta[i])+Iyy*psi_dot[i]*(sin(phi[i]))**2*sin(theta[i])*cos(theta[i])+Izz*psi_dot[i]*(cos(phi[i]))**2*sin(theta[i])*cos(theta[i])],
               [(Iyy-Izz)*psi_dot[i]*(cos(phi[i]))**2*sin(phi[i])*cos(phi[i])-Ixx*theta_dot[i]*cos(theta[i]), (Izz-Iyy)*(theta_dot[i]*cos(phi[i])*sin(phi[i])*sin(theta[i])+phi_dot[i]*(sin(phi[i]))**2*cos(theta[i]))+(Iyy-Izz)*phi_dot[i]*(cos(phi[i]))**2*cos(theta[i])+Ixx*psi_dot[i]*sin(theta[i])*cos(theta[i])-Iyy*psi_dot[i]*(sin(phi[i]))**2*sin(theta[i])*cos(theta[i])-Izz*psi_dot[i]*(cos(phi[i]))**2*sin(theta[i])*cos(theta[i]), (Iyy-Izz)*(phi_dot[i]*cos(phi[i])*sin(phi[i])*(cos(theta[i]))**2)-Iyy*theta_dot[i]*(sin(phi[i]))**2*cos(theta[i])*sin(theta[i])-Izz*theta_dot[i]*(cos(phi[i]))**2*cos(theta[i])*sin(theta[i])+Ixx*theta_dot[i]*cos(theta[i])*sin(theta[i])]             
                ])
                
    
    tau_tilde = dot(linalg.inv(J),(tau.T - dot(c, eta_dot)))

    
#    T = k*(w1[i]**2+w2[i]**2+w3[i]**2+w4[i]**2)
    x_ddot[i+1] = 1/m * T[i] * ( sin(phi[i]) * sin(psi[i]) + cos(phi[i]) * cos(psi[i]) * sin(theta[i]))
    y_ddot[i+1] = 1/m * T[i] * ( cos(phi[i]) * sin(psi[i]) * sin(theta[i]) - sin(phi[i]) * cos(psi[i]) )
    z_ddot[i+1] = T[i] * cos(theta[i]) * cos(phi[i]) - m*g
    
    x_dot[i+1] = x_ddot[i] * dt + x_dot[i]
    y_dot[i+1] = y_ddot[i] * dt + y_dot[i]
    z_dot[i+1] = z_ddot[i] * dt + z_dot[i]

    x[i+1] = x_dot[i] * dt + x[i]
    y[i+1] = y_dot[i] * dt + y[i]
    z[i+1] = z_dot[i] * dt + z[i]

    if z[i+1] < 0:
        z[i+1] = 0
        z_dot[i+1]=0
           
    phi_ddot[i+1] = tau_tilde[0]
    theta_ddot[i+1] = tau_tilde[1]
    psi_ddot[i+1] = tau_tilde[2]
    
    phi_dot[i+1] = phi_ddot[i] * dt + phi_dot[i]; 
    theta_dot[i+1] = theta_ddot[i] * dt + theta_dot[i]; 
    psi_dot[i+1] = psi_ddot[i] * dt + psi_dot[i]; 
    
    eta_dot = array([phi_dot[i], theta_dot[i], psi_dot[i]])
    
    phi[i+1] = phi_dot[i] * dt + phi[i]
    theta[i+1] = theta_dot[i] * dt + theta[i]
    psi[i+1] = psi_dot[i] * dt + psi[i]
    
    eta = array([phi[i+1], theta[i+1], psi[i+1]])
    
#    tau = array([l*k*(w2[i]**2-w4[i]**2),l*k*(w3[i]**2-w1[i]**2), b*(-w1[i]**2+w2[i]**2-w3[i]**2+w4[i]**2)])
    ez[i+1]=(zd[i]-z[i+1])
    ez_dot[i+1]=(ez[i+1]-ez[i])/dt
    T[i+1]=100*ez[i]+200*ez_dot[i]
    
    ephi[i+1]=(phid[i]-phi[i+1])
    ephi_dot[i+1]=(ephi[i+1]-ephi[i])/dt
    tau_phi[i+1]=0*ephi[i]+0*ephi_dot[i]
    
    etheta[i+1]=(thetad[i]-theta[i+1])
    etheta_dot[i+1]=(etheta[i+1]-etheta[i])/dt
    tau_theta[i+1]=10*etheta[i]+etheta_dot[i]
    
    epsi[i+1]=(psid[i]-psi[i+1])
    epsi_dot[i+1]=(epsi[i+1]-epsi[i])/dt
    tau_psi[i+1]=100*epsi[i]+2*epsi_dot[i]
    
    tau = array([tau_phi[i], tau_theta[i], tau_psi[i]])
    print x[i]
    
plot(sin(theta))
plot(sin(thetad))

figure()
plot(y)
