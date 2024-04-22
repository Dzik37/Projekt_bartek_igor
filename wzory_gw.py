from math import *
import numpy as np

a = 6378137.000
e2 = 0.00669438002290

a_kras = 6378245.000 
e2_kras = 0.00669342162296
#promień krzywizny przekroju poprzecznego N 
def N_promien(f, a, e2):
    N = a / np.sqrt(1 - e2 * np.sin(f)**2)
    return(N)

def hirvonen(X, Y, Z, a, e2):
    l = np.arctan2(Y, X)
    p = np.sqrt(X**2 + Y**2)
    f = np.arctan(Z / (p * (1 -e2)))
    while True:
        
        N = N_promien(f, a, e2)
        h = p / np.cos(f) - N
        fs = f
        f = np.arctan(Z / (p * (1 - (e2 * (N / (N + h))))))
        if np.abs(fs - f) < (0.000001/206265):
            break
    return(f, l, h)


# def dms(x):
#     sig = ' '
#     if x<0:
#         sig = '-'
#         x = abs(x)
#     x = x * 180 / pi
#     d = int(x)
#     m = int(60 * (x - d))
#     s = (x - d - m/60) * 3600
#     return(f'{sig}{d:3d}{chr(176)}{abs(m):2d}\'{abs(s):7.5f}\"')
    
    
    
# def dms_azymut(x):
#     sig = ' '
#     if x<0:
#         sig = '-'
#         x = abs(x)
#     x = x * 180 / pi
#     d = int(x)
#     m = int(60 * (x - d))
#     s = (x - d - m/60) * 3600
#     return(f'{sig}{d:3d}{chr(176)}{abs(m):2d}\'{abs(s):7.3f}\"')
   
def dms_azymut(x):
    if x<0:
        x = x + 2 * np.pi
    if x>2*np.pi:
        x = x - 2 * np.pi
    x = x * 180 / pi
    d = int(x)
    m = int(60 * (x - d))
    s = (x - d - m/60) * 3600
    return(f'{d:3d}{chr(176)}{abs(m):2d}\'{abs(s):7.3f}\"')    
def dms(x):
    if x<0:
        x = x + 2 * np.pi
    if x>2*np.pi:
        x = x - 2 * np.pi
    x = x * 180 / pi
    d = int(x)
    m = int(60 * (x - d))
    s = (x - d - m/60) * 3600
    return(f'{d:3d}{chr(176)}{abs(m):2d}\'{abs(s):7.5f}\"')    

def flh2xyz(f, l, h, a, e2):
    X = (N_promien(f, a, e2) + h) * np.cos(f) * np.cos(l)
    Y = (N_promien(f, a, e2) + h) * np.cos(f) * np.sin(l)
    Z = ((N_promien(f, a, e2) * (1 - e2)) + h) * np.sin(f)
    return(X, Y, Z)

#Promień krzywizny przekroju południkowego M 
def Mp(f,a,e2):
    M =(a*(1-e2)) / ((1-e2*(np.sin(f)**2))**(3/2))
    return(M)
    
def stop2dzies(stp,min=0,sek=0):
    stop=stp + min/60 + sek/3600
    return(stop)



def kivioj(f,l,A,s,a,e2):
    n = int(s/1000)
    ds = s/n
    for i in range(n):
        M = Mp(f, a, e2)
        N = N_promien(f, a, e2)
        df = np.cos(A) * ds / M
        dA = (np.sin(A) * np.tan(f) *ds) / N
        fm = f + df/2
        Am = A + dA/2
        Mm = Mp(fm, a, e2)
        Nm = N_promien(fm, a, e2)
        df = (np. cos(Am) * ds )/ Mm
        dA = (np.sin(Am) * np.tan(fm) * ds) / Nm
        dl = (np.sin(Am) * ds ) / (Nm * np.cos(fm))
        f = f + df
        l = l + dl
        A = A + dA
    A = A + pi
    if A > 2*pi:
        A = A - 2*pi
    return(f,l,A)
    


def clairaut(f,A,a,e2):
    N = N_promien(f, a, e2)
    C = N * np.cos(f)* np.sin(A)
    return(C)

def vincenty(fa, la, fb, lb, a, e2):
    b = a * sqrt((1 - e2))
    fl = 1- (b/a)
    dL = lb - la
    Ua = atan((1 - fl) * tan(fa))
    Ub = atan((1 - fl) * tan(fb))
    L = dL
    while True:
        sin_sig = sqrt((cos(Ub) * sin(L)) ** 2 + (cos(Ua) * sin(Ub) - sin(Ua) * cos(Ub) * cos(L)) ** 2)
        cos_sig = sin(Ua) * sin(Ub) + cos(Ua) * cos(Ub) * cos(L)
        sigma = atan2(sin_sig, cos_sig)
        sin_a = cos(Ua) * cos(Ub) * sin(L) / sin_sig
        cos2_a = 1 - sin_a ** 2
        cos2sm = cos_sig - 2 * sin(Ua) * sin(Ub) / cos2_a
        C = fl / 16 * cos2_a * (4 + fl * (4 - 3 * cos2_a))
        Ls = L
        L = dL + (1 - C) * fl* sin_a * (sigma + C * sin_sig * (cos2sm + C * cos_sig * (-1 + 2 * cos2sm ** 2)))
        if abs(Ls - L) < (0.000001 / 206265):
            break
    u2 = (a ** 2 - b ** 2) / b ** 2 * cos2_a 
    A = 1 + u2 / 16384 * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
    B = u2 / 1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
    ds = B * sin_sig * (cos2sm + B / 4 * (cos_sig * (-1 + 2 * cos2sm ** 2) - B / 6 * cos2sm * (-3 + 4 * sin_sig ** 2) * (-3 + 4 * cos2sm ** 2)))
    s = b * A * (sigma - ds)
    Aab = atan2(cos(Ub) * sin(L), cos(Ua) * sin(Ub) - sin(Ua) * cos(Ub) * cos(L))
    Aba = (atan2(cos(Ua) * sin(L), -sin(Ua) * cos(Ub) + cos(Ua) * sin(Ub) * cos(L))) + pi
    return (s, Aab, Aba)

def saz2neu(s, alfa, z):
    n = s* np.sin(z) * np.cos(alfa)
    e = s* np.sin(z) * np.sin(alfa)
    u = s* np.cos(z)
    dx = np.array([n,e,u])
    return(dx)


def neu2saz(dx):
    s = np.sqrt(dx[0]**2 + dx[1]**2 + dx[2]**2)
    alfa = np.arctan2(dx[1], dx[0])
    z = np.arccos(dx[2] / s)
    return(s,alfa,z)

def Rneu(fa,la):
    R = np.array([[-np.sin(fa)*np.cos(la) , -np.sin(la) , np.cos(fa)*np.cos(la)],
                  [-np.sin(fa)*np.sin(la) , np.cos(la) , np.cos(fa)*np.sin(la)],
                  [np.cos(fa) , 0 , np.sin(fa)]])
    return(R)

def neu2XYZ(dx,f, l):
    R=Rneu(f, l)
    dX= R @ dx
    return (dX)

def XYZ2neu(dX,f,l):
    R=Rneu(f, l)
    dx = R. T @ dX
    return(dx)
    

def stopxyz(h, f, l,a, e2):#zamiana ze wsp geo na wsp prost
    N = Np(f, a, e2)
    A = (N + h) * np.cos(f) * np.cos(l)
    B = (N + h) * np.cos(f) * np.sin(l)
    C = (N * (1 - e2)+h)*np.sin(f)
    return A, B, C

def sigma(phi,a,e2):
    A_0 = 1 - e2/4 - (3 * e2 * e2)/64 - (5 * e2 * e2 *e2)/256
    A_2 = 3/8 * (e2 + (e2 *e2)/4 + (15 * e2 * e2 * e2)/128)
    A_4 = 15/256 * (e2 * e2 + (3 * e2 * e2 * e2)/4)
    A_6 = (35 * e2 * e2 * e2)/3072

    sig = a * (A_0 * phi - A_2 * np.sin(2 * phi) + A_4 * np.sin(4 * phi) - A_6 * np.sin(6 * phi))
    return(sig)



def fl2gk(phi,lam,lam_0,a,e2):

    b2 = a**2 * (1 - e2)
    e22 = (a**2 - b2)/b2


    dlam = lam - lam_0
    t = np.tan(phi)
    eta2 = e22 * (np.cos(phi)**2)
    N = N_promien(phi, a, e2)
    sig = sigma(phi, a, e2)
   
    X_gk = sig+(dlam**2/2)*N*np.sin(phi)*np.cos(phi)*((1+(dlam**2/12)*(np.cos(phi))**2*(5-t**2+9*eta2+4*eta2**2)+(dlam**4/360)*np.cos(phi)**4*(61-58*t**2+t**4+270*eta2-330*eta2*t**2)))
    Y_gk = dlam * N * np.cos(phi) * (1 + (dlam**2/6) * (np.cos(phi)**2) * (1 - t**2 + eta2) + (dlam**4/120) * (np.cos(phi)**4) * (5 - 18 * t**2 + t**4 + 14 * eta2 - 58 * eta2 * t**2))
    return(X_gk,Y_gk)

def f1(X_gk,a,e2):
    A_0 = 1 - e2/4 - 3 * e2 **2/64 - 5 * e2 * e2**2/256
    phi_1 = X_gk/(a *A_0)
    while True:
        phi_i = phi_1
        s = sigma(phi_1,a,e2)
        phi_1 = phi_1 + (X_gk - s)/(a * A_0)
        if abs(phi_1 - phi_i) < (0.000001/206265):
                break
    return(phi_1)








def gk2fl(X_gk, Y_gk, lam_0, a, e2):
    b2 = a**2 * (1 - e2)
    e22 = (a**2 - b2)/b2
    f_1 = f1(X_gk, a, e2)
    eta2 = e22 * (np.cos(f_1)**2)
   
    M = Mp(f_1, a, e2)
    N = N_promien(f_1, a, e2)
    t = np.tan(f_1)
    b2 = a**2 * (1-e2)
    e2p = (a**2 - b2)/b2
    n2 = e2p * ((np.cos(f_1))**2)
    phi = f_1 - ((Y_gk**2) * t)/(2 * M * N) * (1 - ((Y_gk)**2)/(12 * N**2) * (5 + 3 * (t**2) + eta2 - 9 * eta2 * (t**2) - 4 * (eta2**2)) + ((Y_gk)**4)/(360 * N**4) * (61 + 90 * (t**2) + 45 * (t**4)))
    lam = lam_0 + Y_gk/(N * np.cos(f_1)) * (1 - (Y_gk**2)/(6 * N**2) * (1 + 2 * (t**2) + eta2) + (Y_gk**4)/(120 * (N**4)) * (5 + 28 * (t**2) + 24 * (t**4) + 6 * eta2 + 8 * eta2 * (t**2)))


    return(phi, lam)
    
def m_gk(x,y,a,e2):
    f_1 = f1(x, a, e2)   
    R1 = np.sqrt(Mp(f_1, a, e2)*N_promien(f_1, a, e2))
    m_gk = 1 + (y**2 /( 2*(R1**2))) + (y**4 / (24*R1**4))
    return(m_gk)


def prom_m(fa, fb):
    fm= (fa + fb) / 2
    M = Mp(fm, a, e2)
    N = N_promien(fm, a, e2)
    Rm = np.sqrt(M * N)
    return (Rm)




def dlug_s(spom, m0, yagk, ha, ybgk, hb, Rm, a, e2):
    s0 = np.sqrt(( (spom**2) - ((hb - ha)**2) )/( (1 + (ha/Rm)) * (1 + (hb/Rm)) ))
    selip = 2 * Rm * np.arcsin(s0/(2*Rm))
    r = spom * ((yagk**2 + yagk * ybgk + ybgk**2) / (6* Rm**2))
    sgk = selip + r 
    s = sgk * m0
    return(r, selip, sgk, s0 ,s)
    
def gamma(xgk, ygk, l0, a, e2):
    f = f1(xgk, a, e2)
    N = N_promien(f, a, e2)
    M = Mp(f, a, e2)
    t = np.tan(f)
    b2 = (a * np.sqrt(1-e2))**2
    e22 = (a**2 - b2)/b2
    ni2 = e22 * (np.cos(f))**2
    g = ((ygk/N)  * t) * (1 - (ygk**2/(3 * N**2)) * (1 + t**2 - ni2 - 2 * ni2**2) + (ygk**4/(15 * N**4)) * (2 + 5 * t**2 + 3 * t**4))
    return(g)

def delta(xa, ya, xb, yb, l0, a, e2):
    xr = (xa+xb)/2
    yr = (ya + yb)/2
    f, l = gk2fl(xr, yr, l0, a, e2)
    M = Mp(f, a, e2)
    N = N_promien(f, a, e2)
    R = np.sqrt(M * N)
    dab = (xb - xa) * (2 * ya + yb)/(6 * R**2)
    dba = (xa - xb) * (2 * yb + ya)/(6 * R**2)
    return(dab, dba)
    


def trans(alfa=0, beta=0, gam=0, kx=0, ky=0, kz=0, X0=0, Y0=0, Z0=0, X=0, Y=0, Z=0):
    B = np.array([[cos(gam)*cos(beta), sin(gam), -sin(beta)],
                  [-sin(gam), cos(alfa)*cos(gam), sin(alfa)],
                  [sin(beta), -sin(alfa), cos(beta)*cos(gam)]])
    #B = np.array([[1, gam, -beta],  #DLA MAŁYCH KĄTÓW (<5")
                  #[-gam, 1, alfa],
                  #[beta, -alfa, 1]])
    macierz_zniek = np.array([[kx, 0, 0],
                              [0, ky, 0],
                              [0, 0, kz]])
    MB = np.array([[kx, gam, -beta],
                   [-gam, ky, alfa],
                   [beta, -alfa, kz]])
    E = np.array([[1,0,0],
                  [0,1,0],
                  [0,0,1]])
    M = E + macierz_zniek
    rp = np.array([[X],
                   [Y],
                   [Z]])
    skala_ukladu = M @ rp
    r0 = np.array([[X0],
                   [Y0],
                   [Z0]])
    przesuniecie_ukladu = rp + r0
    rpp_quasi_afini_male_katy = rp + MB @ rp + r0
    return(rpp_quasi_afini_male_katy)

def od2000 (x2000,y2000, nr):
    m2000 = 0.999923
    xgk = x2000/m2000
    ygk = ((y2000 - 500000)-(nr*1000000))/m2000
    return(xgk, ygk)  
    
def od1992 (X92,Y92):
    M1992=0.9993
    Xgk = (X92+5300000)/M1992
    Ygk = (Y92-500000)/M1992
    return(Xgk,Ygk)   
    
def prom_xy(xa, ya, xb, yb, l_0): 
    xm = (xa + xb)/2
    ym = (ya + yb)/2
    fm, lm = gk2fl(xm, ym, l_0, a, e2)
    M = Mp(fm, a,e2)
    N = N_promien(fm, a ,e2)
    Rm = np.sqrt(M * N)
    return (Rm)   


def spom(X1,Y1,X2,Y2):
    dlugosc = np.sqrt(((X2-X1)**2) + ((Y2-Y1)**2))
    return(dlugosc)


def uklad_2000 (xgk,ygk,nr):
    m2000=0.999923
    x2000=xgk*m2000
    y2000=ygk*m2000+nr*1000000+500000
    return(x2000,y2000)

def uklad_1992 (xgk,ygk):
    m1992=0.9993
    x1992=xgk * m1992- 5300000
    y1992=ygk * m1992 + 500000
    return(x1992,y1992)




