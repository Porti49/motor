import numpy as np
import matplotlib.pyplot as plt

rg = 287

def work12(T1,P1,rc,y):
    work12_s = (1-rc**(y-1))/(y-1)
    T2 = T1 * rc**(y-1)
    P2 = P1 * rc**y
    return work12_s, T2, P2

def heat23(T2,P2,rc,rp,y,y2,ypc):
    q23_s = rc**(y-1)*(rp/(ypc-1)-1/(y2-1)*mp)
    T3 = T2 * (rp)
    P3 = P2 * rp
    return q23_s,T3,P3

def work34(T3,P3,rc,rp,rv,y):
    work34_s = rc**(y-1)*rp*(rv-1)
    T4 = T3 * rv
    P4 = P3
    return work34_s,T4,P4

def heat34(rc,rp,rv,y,ypc):
    q34_s = rc**(y-1)*(rp*(rv-1)*ypc/(ypc-1))
    return q34_s

def work45(T4, P4, rc, re, rp, rv, yp, y):
    w45_s = (rc**(y-1)*rv*rp/(yp-1))*(1-(rv/re)**(yp-1))
    T5 = T4 * (rv/re)**(yp-1)
    P5 = P4 * (rv/re)**yp
    return w45_s,T5,P5

def work67(re,rc,rad_p,pi_re):
    w67_s = (re-1)/(rc*rad_p)
    return w67_s
    
def work81(rc,pi_re):
    w81_s = (rc-1)/rc
    return w81_s 

def ciclo(m,T1,P1,rc,re,rp,rv,y,y2,ypc,yp,yr,rad_p,pi_re,nu,q):

    #Volumenes
    V2 = q / (re-1)
    V1 = V2 * rc
    V3 = V2
    V4 = V3 * rv
    V5 = V3 * re
    V6 = V5
    V7 = V2
    V8 = V2
    
    #Trabajo 1 a 2
    w12 = work12(T1,P1,rc,y)
    w12_r = w12[0]#*(m*rg*T1)
    T2 = w12[1]
    P2 = w12[2]
    
    #Calor 2 a 3
    q23 = heat23(T2,P2,rc,rp,y, y2, ypc)
    q23_r = q23[0]#*(m*rg*T1)
    T3 = q23[1]
    P3 = q23[2]
    
    #Trabajo 3 a 4
    w34 = work34(T3,P3,rc,rp,rv,y)
    w34_r = w34[0]#*(m*rg*T1)
    T4 = w34[1]
    P4 = w34[2]
    
    #Calor 3 a 4
    q34 = heat34(rc,rp,rv,y,ypc)
    q34_r = q34#*(m*rg*T1)
    
    #Trabajo 4 a 5
    w45 = work45(T4, P4, rc, re, rp, rv, yp, y)
    w45_r = w45[0]#*(m*rg*T1)
    T5 = w45[1]
    P5 = w45[2]
    
    #5 a 6
    Pes_p = P1*rad_p
    P6 = Pes_p
    T6 = T5*(Pes_p/P5)**((yr-1)/yr)
    
    
    #Trabajo 6 a 7
    w67 = work67(re,rc,rad_p,pi_re)
    w67_r = w67#*(m*rg*T1)
    
    P7 = P6
    Pad_p = (1/rad_p) * Pes_p
    P8 = Pad_p
    
    #Trabajo 8 a 1
    w81 = work81(rc,pi_re)
    w81_r = w81#*(m*rg*T1)
    
    #Temperatura gases de escape
    Tr = T1 * nu*((rc**(y/yr-1)*(rp*rv**yp)**(1/yr)) / (re**(yp/yr-1)*mp*(rad_p*pi_re)**(1-1/yr)))
    T7 = Tr
    
    
    V = V1, V2, V3, V4, V5, V6, V7, V8
    P = P1, P2, P3, P4, P5, P6, P7, P8
    T = T1, T2, T3, T4, T5, T6, T7
    W = w12_r + w34_r + w45_r - (w67_r - w81_r)*(1/pi_re)
    Q = q23_r + q34_r
    return W, Q, T, P, V

#Cilindrada
q = 1.4*10**-3

#Ciclo = f(m, T1, P1, rc, re, rp, rv, y, y2, ypc, yp, yr, rad_p, pi_re, nu, q)
m = 0.183
#mp/m
mp = 1.0
Tad_re = 273
Patm = 1
#Parametros de presiones internos 
rc = 8
re = 10
rp = 1
rv = 2
y = 1.3
ypc = 1.3
yp = 1.3
yr = 1.3
y2 = y
yf = 1.3
y1 = y
nu = 1

#Parametros de presion previos y turbo
pi_def = 1
pi_fi = 1
pi_ctt = 2
pi_in = 1
mariposa = 1
pi_vad = 1.0
rad_p = 1.05
pi_re = 1.05
P1 = Patm*pi_def*pi_fi**-1*pi_ctt*pi_in**-1*mariposa*pi_vad**-1*pi_re
T1 = Tad_re*(1-1/(nu*rv*(rc**y*rad_p*pi_re)**(1/yp)))*yf/(1+(y1-1)/pi_re*(1-1/rc)*(y1-1)/(yr-1)*(1/(rc*rad_p*pi_re)))
Fr = 1
re_op = rv*(rc**y*rp*rad_p*pi_re)**(1/yp)
f = (1/(nu*re_op)*mp)
EGRa = 0
yh = 0

#rendimientos mecanicos


Ciclo = ciclo(m, T1, P1, rc, re, rp, rv, y, y2, ypc, yp, yr, rad_p, pi_re, nu, q)
W = Ciclo[0]
Q = Ciclo[1]
T = Ciclo[2]
P = Ciclo[3]
V = Ciclo[4]

#Lamdas en funcion de T
T_p = np.array([(T[1]+T[3])/2, (T[3]+T[4])/2]) #ypc, yp
T_r = np.array([(T[0]+T[1])/2, T[5]/2]) #y = y2, yr 
y_p = 1.411-0.03*Fr-0.067*T_p/1000
y_r = 1.43 - 0.095*T_r/1000-0.045*Fr-0.01*(f+EGRa+yh)  

while abs(y_p[0]-ypc)>0.0001 and abs(y_p[1]-yp)>0.0001 and abs(y_r[1]-yr)>0.0001 and (y_r[0]-y)>0.0001:
    ypc = y_p[0]
    yp = y_p[1]
    y = y_r[0]
    yr = y_r[1]
    y1 = y
    y2 = y
    Ciclo = ciclo(m, T1, P1, rc, re, rp, rv, y, y2, ypc, yp, yr, rad_p, pi_re, nu, q)
W = Ciclo[0]
Q = Ciclo[1]
T = Ciclo[2]
P = Ciclo[3]
V = Ciclo[4]

W_t = W * T1 * m * rg

#re a rendimiento maximo ceteris paribus
rendimiento_p = 0
rendimiento = W / Q
re_v = ()
rend_p = ()
rend_i = 0

while rend_i<=rendimiento_p:
    re = re + 0.001
    rend_i = rendimiento_p
    Ciclo = ciclo(m, T1, P1, rc, re, rp, rv, y, y2, ypc, yp, yr, rad_p, pi_re, nu, q)
    rendimiento_p = Ciclo[0]/Ciclo[1]
    re_v = np.append(re_v,re)
    rend_p = np.append(rend_p, rendimiento_p)

pmt = W_t/q


plt.plot(V,P,'bo-',linewidth = 1, markersize = 2.5)
plt.yscale('log')
plt.xscale('log')
plt.xlabel("Volumen")
plt.ylabel("Presion")
plt.title("Ciclo teórico de dos composiciones")



"""
#Esta hecho con re/2 para que se pueda representar los valores de 0.5 en 0.5
W = ()
re_v = ()
rend = ()
for rp in range(2,7,2):
    for re in range(10, 41):
        Ciclo = ciclo(m, T1, P1, rc, re/2, rp, rv, y, y2, ypc, yp, yr, rad_p, pi_re, nu, q)
        rendimiento = Ciclo[0]/Ciclo[1]
        W = np.append(W,Ciclo[0])
        re_v = np.append(re_v,re)
        rend = np.append(rend,rendimiento)
  
    plt.plot(re_v/2,rend,"bo",markersize = 2.5)
    plt.title("Rendimiento: rc fija")
    plt.xlabel("r.e")
    plt.ylabel("rendimiento")

 """  


