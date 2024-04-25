import numpy as np
from math import sin, cos, sqrt, atan, degrees, pi, tan
import sys

class Transformacje:
    def __init__(self, model: str = "wgs84"):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            ecc2 - mimośród^2
        + WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        + Inne powierzchnie odniesienia: https://en.wikibooks.org/wiki/PROJ.4#Spheroid
        + Parametry planet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html
        """
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "mars":
            self.a = 3396900.0
            self.b = 3376097.80585952
        elif model == "krasowski":
            self.a = 6378245.000
            self.b = 6356863.019            
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flat = (self.a - self.b) / self.a
        self.ecc = sqrt(2 * self.flat - self.flat ** 2) # eccentricity  WGS84:0.0818191910428 
        self.ecc2 = (2 * self.flat - self.flat ** 2) # eccentricity**2


    
    def xyz2plh(self, X, Y, Z, output = 'dec_degree'):
        """
        Przeliczenie współrzędnych X, Y, Z do phi, lambda, h:
        
        Parameters
        ----------
        X, Y, Z : FLOAT
             współrzędne w układzie orto-kartezjańskim, 

        Returns
        -------
        lat
            [stopnie dziesiętne] - szerokość geodezyjna
        lon
            [stopnie dziesiętne] - długośc geodezyjna.
        h : TYPE
            [metry] - wysokość elipsoidalna
        output [STR] - optional, defoulf 
            dec_degree - decimal degree
            dms - degree, minutes, sec
        """
        r   = sqrt(X**2 + Y**2)           # promień
        lat_prev = atan(Z / (r * (1 - self.ecc2)))    # pierwsze przybliilizenie
        lat = 0
        while abs(lat_prev - lat) > 0.000001/206265:    
            lat_prev = lat
            N = self.a / sqrt(1 - self.ecc2 * sin(lat_prev)**2)
            h = r / cos(lat_prev) - N
            lat = atan((Z/r) * (((1 - self.ecc2 * N/(N + h))**(-1))))
        lon = atan(Y/X)
        N = self.a / sqrt(1 - self.ecc2 * (sin(lat))**2);
        h = r / cos(lat) - N       
        if output == "dec_degree":
            return degrees(lat), degrees(lon), h




    def plh2xyz(self, phi, lam, h):
        """ 
        Przeliczenie współrzędnych phi, lambda, h do X, Y, Z:
        
        Parameters:
        -----------
        FI, LAM, H : FLOAT
            współrzędne geodezyjne i wysokosc elipsoidalna,
            współrzędne - [stopnie], wysokosc [metry]
        
        
        Returns:
        --------
        Współrzędne ortokartezjańskie:
        X : FLOAT
            [metry]
        Y : FLOAT
            [metry]
        Z : FLOAT
            [metry]
            
        
        """
        N =  self.a / sqrt(1 - self.ecc2 * sin(phi)**2)
        X = (N + h) * cos(phi) * cos(lam)
        Y = (N + h) * cos(phi) * sin(lam)
        Z = (N * (1-self.ecc2) + h) * sin(phi) 
        return X, Z, Y



    def xyz2neu(self, X, Y, Z, X0, Y0, Z0):
        """
        Wyznaczenie współrzędnych topocentrycznych N,E,U:
            
        Parameters
        ----------
        X, Y, Z : FLOAT
            współrzędne geocentryczne satelitów
            [metry]
        X0, Y0, Z0 : FLOAT
            współrzędne geocentryczne anteny
            [metry]
        
        Returns
        -------
        Współrzędne topocentryczne satelitów:
        n : FLOAT
            [metry]
        e : FLOAT
            [metry]
        u : FLOAT
            [metry]
        

        """
        dX = np.array([[X - X0],
                       [Y - Y0],
                       [Z - Z0]])

        def xyz2plh(self, X, Y, Z, output = 'dec_degree'):
            r   = sqrt(X**2 + Y**2)           # promień
            lat_prev = atan(Z / (r * (1 - self.ecc2)))    # pierwsze przybliilizenie
            lat = 0
            while abs(lat_prev - lat) > 0.000001/206265:    
                lat_prev = lat
                N = self.a / sqrt(1 - self.ecc2 * sin(lat_prev)**2)
                h = r / cos(lat_prev) - N
                lat = atan((Z/r) * (((1 - self.ecc2 * N/(N + h))**(-1))))
            lon = atan(Y/X)
            N = self.a / sqrt(1 - self.ecc2 * (sin(lat))**2);
            h = r / cos(lat) - N       
            
            return lat, lon, h
        
        phi, lam, h = xyz2plh(self,X,Y,Z)

        def RneuT (phi,lam):
            R = np.array([[-np.sin(phi) * np.cos(lam),-np.sin(phi) * np.sin(lam), np.cos(phi)],
                          [-np.sin(lam), np.cos(lam), 0],
                          [np.cos(phi) * np.cos(lam), np.cos(phi) * np.sin(lam), np.sin(phi)]])
            return R
            
        n,e,u = RneuT(phi,lam) @ dX 
        
        return (n, e, u)
       



    def uklad2000(self, fi,lam,h):
        """
        Wyznaczenie współrzędnych w układzie 2000 
        ze współrzędnych geodezyjnych phi, lambda, h.

        Parameters
        ----------
        FI, LAM, H : FLOAT
            współrzędne geodezyjne i wysokosc elipsoidalna,
            współrzędne - [stopnie], wysokosc [metry]

        Returns
        -------
        Współrzędne w układzie 2000:
        X2000 : FLOAT
            [metry]
        Y2000 : FLOAT
            [metry]


        """
        
        if lam < 16.5 :
                L0 = 15
                nrst = 5
        if lam > 16.5 and lam < 19.5 :
                    L0 = 18
                    nrst = 6
        if lam > 19.5 and lam < 22.5:
                L0 = 21
                nrst = 7
        if lam > 22.5:
                L0 = 24
                nrst = 8
        #delta lambda
        l = lam-L0
        #zmiana jednostek
        fi = fi * pi/180
        lam = lam * pi/180
        l = l * pi/180
    
        b2 = self.a**2 *(1-self.ecc2)
        e_2 = ((self.a)**2-(b2))/((b2))
    
        t = tan(fi)
    
        eta2 = (e_2)*((cos(fi))**2)
    
        N = self.a/(sqrt( 1 -(self.ecc2) * (sin(fi))**2))
    
    
        A0 = 1 - (self.ecc2/4) - (3 * ((self.ecc2)**2)/64) - (5*((self.ecc2)**3)/256)
        A2 = (3/8) * ((self.ecc2) + (((self.ecc2)**2)/4) + (15 * ((self.ecc2)**3)/128))
        A4 = (15/256) * (((self.ecc2)**2)+(3*((self.ecc2)**3)/4))
        A6 = (35*(self.ecc2)**3)/3072
    
        sigma = self.a * ((A0) * (fi) - (A2) * (sin(2 * fi)) + (A4) * (sin(4 * fi)) - (A6) * (sin(6*fi)))
    
    
        Xgk = sigma + ((l**2)/2) * N * (sin(fi)) * (cos(fi)) * (1 + ((l**2)/12) * ((cos(fi))**2) * (5 - (t**2) + 9*(eta2) + 4*(eta2)**2) + ((l**4)/360) * ((cos(fi))**4) * (61 - 58*(t**2) + (t**4) + 270*(eta2) - 330*(eta2) * (t**2)))
        Ygk = l * N * (cos(fi)) * (1 + ((l**2)/6) * ((cos(fi))**2) * (1 - t**2 + eta2) + ((l**4)/120) * ((cos(fi))**4) * (5 - 18*(t**2) + (t**4) + 14*eta2 - 58*(eta2) * (t**2)))
    
        mo=0.999923
        X2000 = Xgk * mo
        #Y2000 = Ygk * mo + nr strefy * 1 000 000 + 500 000
        
        Y2000 = Ygk * mo + nrst * 1000000 + 500000
        return X2000, Y2000



    def uklad1992(self,fi,lam,h):
        """
        Wyznaczenie współrzędnych w układzie 1992 
        ze współrzędnych geodezyjnych phi, lambda, h.

        Parameters
        ----------
        FI, LAM, H : FLOAT
            współrzędne geodezyjne i wysokosc elipsoidalna,
            współrzędne - [stopnie], wysokosc [metry]
        

        Returns
        -------
        Współrzędne w układzie 1992:
        X92 : FLOAT
            [metry]
        Y92 : FLOAT
            [metry]

        """
    
        #delta lambda
        L0 = 19
        l = lam-L0
        #zmiana jednostek
        fi = fi * pi/180
        lam = lam * pi/180
        l = l * pi/180
    
        b2 = self.a**2 *(1-self.ecc2)
        e_2 = ((self.a)**2-(b2))/((b2))
    
        t = tan(fi)
    
        eta2 = (e_2)*((cos(fi))**2)
    
        N = self.a/(sqrt( 1 -(self.ecc2) * (sin(fi))**2))
    
    
        A0 = 1 - (self.ecc2/4) - (3 * ((self.ecc2)**2)/64) - (5*((self.ecc2)**3)/256)
        A2 = (3/8) * ((self.ecc2) + (((self.ecc2)**2)/4) + (15 * ((self.ecc2)**3)/128))
        A4 = (15/256) * (((self.ecc2)**2)+(3*((self.ecc2)**3)/4))
        A6 = (35*(self.ecc2)**3)/3072
    
        sigma = self.a * ((A0) * (fi) - (A2) * (sin(2 * fi)) + (A4) * (sin(4 * fi)) - (A6) * (sin(6*fi)))
    

        Xgk = sigma + ((l**2)/2) * N * (sin(fi)) * (cos(fi)) * (1 + ((l**2)/12) * ((cos(fi))**2) * (5 - (t**2) + 9*(eta2) + 4*(eta2)**2) + ((l**4)/360) * ((cos(fi))**4) * (61 - 58*(t**2) + (t**4) + 270*(eta2) - 330*(eta2) * (t**2)))
        Ygk = l * N * (cos(fi)) * (1 + ((l**2)/6) * ((cos(fi))**2) * (1 - t**2 + eta2) + ((l**4)/120) * ((cos(fi))**4) * (5 - 18*(t**2) + (t**4) + 14*eta2 - 58*(eta2) * (t**2)))
    
        mo=0.9993
        X92 = Xgk * mo - 5300000
        Y92 = Ygk * mo + 500000
        return X92, Y92


    
if __name__ == "__main__":
    # # utworzenie obiektu
    # geo = Transformacje(model = "mars")
    # # dane XYZ geocentryczne
    # X = 3664940.500; Y = 1409153.590; Z = 5009571.170
    # X2 = 3664940.520; Y2 =1409153.570; Z2 =5009571.167 
    # phi, lam, h = geo.xyz2plh(X, Y, Z)
    # print(phi, lam, h)
    # x,y,z = geo.plh2xyz(phi, lam, h)
    # print(x,y,z)

    #URUCHAMIANIE PROGRAMU Z KONSOLI, w konsole trzeba wpisać: python projekt_Bartek_igor.py xyz2plh wsp_xyz_inp.txt

    geo = Transformacje(model = "wgs84")
    
    sys.argv[0] = 'nazwa programu'
    wybrana_funkcja = sys.argv[1]
    plik = sys.argv[2]
    
    print(sys.argv[0], sys.argv[1], sys.argv[2])

    if wybrana_funkcja == 'xyz2plh':
        with open(plik, 'r') as f:
            lines = f.readlines()
            coord_lines = lines
        
            coords_plh = []
            for coord_line in coord_lines:
                coord_line = coord_line.strip('\n')
                x_str, y_str, z_str = coord_line.split(',')
                x, y, z = (float(x_str), float(y_str), float(z_str))
                phi, lam, h = geo.xyz2plh(x, y, z)
                coords_plh.append([phi, lam, h])
        
        file_out = 'wsp_plh_out.txt'
        f1 = open(file_out, 'w')
        
        for plh in coords_plh:
            s = f'{plh[0]:.5f},{plh[1]:.5f},{plh[2]:.3f} \n'
            f1.write(s)
            
        f1.close()    












    
# class Transformacje:
#     def __init__(self, model: str = "wgs84"):
#         """
#         Parametry elipsoid:
#             a - duża półoś elipsoidy - promień równikowy
#             e2 - mimośród^2
#         """
#         if model == "wgs84":
#             self.a = 6378137.0 
#             self.e2 = 0.006694379990141124
#         elif model == "grs80":
#             self.a = 6378137.0
#             self.e2 = 0.00669438002290
#         elif model == "krasowski":
#             self.a = 6378245.000 
#             self.e2 = 0.00669342162296
#         else:
#             raise NotImplementedError(f"{model} model not implemented")
       
#     def xyz2blh(self,X, Y, Z, output = 'f,l,h'):
#         l = np.arctan2(Y, X)
#         p = np.sqrt(X**2 + Y**2)
#         f = np.arctan(Z / (p * (1 -self.e2)))
#         while True:
        
#             N = self.a / np.sqrt(1 - self.e2 * np.sin(f)**2)
#             h = p / np.cos(f) - N
#             fs = f
#             f = np.arctan(Z / (p * (1 - (self.e2 * (N / (N + h))))))
#             if np.abs(fs - f) < (0.000001/206265):
#                 break
#             if output == "dec_degree":
#                 return f,l,h 
#             elif output == "dms":
#                 def dms(x):
#                     if x<0:
#                         x = x + 2 * np.pi
#                     if x>2*np.pi:
#                         x = x - 2 * np.pi
#                         x = x * 180 / pi
#                         d = int(x)
#                         m = int(60 * (x - d))
#                         s = (x - d - m/60) * 3600
#                         return(f'{d:3d}{chr(176)}{abs(m):2d}\'{abs(s):7.5f}\"')  
#                     return(dms(f),dms(l),h)
        
#         else:
#             raise NotImplementedError(f"{output} - output format not defined")
# if __name__ == "__main__":
#     # utworzenie obiektu
#     geo = Transformacje(model = "wgs84")
#     # dane XYZ geocentryczne
#     X = 3664940.500; Y = 1409153.590; Z = 5009571.170
#     phi, lam, h = geo.xyz2blh(X, Y, Z)
#     print(phi, lam, h)