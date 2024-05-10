import numpy as np
from math import sin, cos, sqrt, atan, degrees, pi, tan,radians
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
        elif model == "krasowski":
            self.a = 6378245.000
            self.b = 6356863.019            
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flat = (self.a - self.b) / self.a
        self.ecc = sqrt(2 * self.flat - self.flat ** 2) # eccentricity  WGS84:0.0818191910428 
        self.ecc2 = (2 * self.flat - self.flat ** 2) # eccentricity**2


    
    def xyz2plh(self, X, Y, Z):
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
        return X, Y, Z



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
            
        [[n],[e],[u]] = RneuT(phi,lam) @ dX 
        
        return n, e, u
       



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
    #URUCHAMIANIE PROGRAMU Z KONSOLI, w konsole trzeba wpisać: python transformacje.py --xyz2plh wsp_xyz_inp.txt --wgs84
    #DLA FUNKCJI xyz2neu: python transformacje.py --xyz2neu wsp_xyz_inp.txt --wgs84 3664940.500 1409153.590 5009571.170

   
    
    sys.argv[0] = 'transformacje'
    wybrana_funkcja = sys.argv[1]
    plik = sys.argv[2]
    elipsoida = sys.argv[3]
    
    if elipsoida == "--wgs84":
        geo = Transformacje(model = "wgs84" )
    elif elipsoida == "--grs80":
        geo = Transformacje(model =  "grs80")
    elif elipsoida == "--krasowski":
        geo = Transformacje(model = "krasowski")
    else:
        raise NameError('nie ma takiej elipsoidy, wybierz wgs84, grs80, mars lub krasowski')
     
    print(f'program: {sys.argv[0]}\nwybrana funkcja: {sys.argv[1]}\nplik wejsciowy: {sys.argv[2]}\nwybrana elipsoida: {sys.argv[3]}')

    if wybrana_funkcja == '--xyz2plh':
        
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

            
        print(f'plik wyjsciowy: wsp_plh_out.txt')
            

        
    elif wybrana_funkcja == '--plh2xyz':
        with open(plik, 'r') as f:
            lines = f.readlines()
            coord_lines = lines
        
            coords_xyz = []
            for coord_line in coord_lines:
                coord_line = coord_line.strip('\n')
                p_str, l_str, h_str = coord_line.split(',')
                p, l, h = (radians(float(p_str)), radians(float(l_str)), float(h_str))
                x, y, z = geo.plh2xyz(p, l, h)
                coords_xyz.append([x, y, z])
        
        file_out = 'wsp_xyz_out.txt'
        f1 = open(file_out, 'w')
        
        for xyz in coords_xyz:
            s = f'{xyz[0]:.3f},{xyz[1]:.3f},{xyz[2]:.3f} \n'
            f1.write(s)
            
        f1.close() 
        
        print(f'plik wyjsciowy: wsp_xyz_out.txt')
        
    
    elif wybrana_funkcja == '--xyz2neu':
        try:
            x0 = sys.argv[4]
            y0 = sys.argv[5]
            z0 = sys.argv[6]
            x0, y0, z0 = float(x0), float(y0), float(z0)
            
            with open(plik, 'r') as f:
                lines = f.readlines()
                coord_lines = lines
            
                coords_neu = []
                for coord_line in coord_lines:
                    coord_line = coord_line.strip('\n')
                    x_str, y_str, z_str = coord_line.split(',')
                    x, y, z = (float(x_str), float(y_str), float(z_str))
                    n, e, u = geo.xyz2neu(x, y, z, x0, y0, z0)
                    coords_neu.append([n, e, u])
            
            file_out = 'wsp_neu_out.txt'
            f1 = open(file_out, 'w')
            
            for neu in coords_neu:
                s = f'{neu[0]:.3f},{neu[1]:.3f},{neu[2]:.3f} \n'
                f1.write(s)
                
            f1.close()  
            
            print(f'plik wyjsciowy: wsp_neu_out.txt')
        except IndexError:
            raise AttributeError('nie podano współrzędnych geocentrycznych: x0, y0, z0')
        
    elif wybrana_funkcja == '--uklad2000':
        with open(plik, 'r') as f:
            lines = f.readlines()
            coord_lines = lines
        
            coords_xy2000 = []
            for coord_line in coord_lines:
                coord_line = coord_line.strip('\n')
                p_str, l_str, h_str = coord_line.split(',')
                p, l, h = (float(p_str), float(l_str), float(h_str))
                x, y = geo.uklad2000(p, l, h)
                coords_xy2000.append([x, y])
        
        file_out = 'wsp_xy2000_out.txt'
        f1 = open(file_out, 'w')
        
        for xy in coords_xy2000:
            s = f'{xy[0]:.3f},{xy[1]:.3f} \n'
            f1.write(s)
            
        f1.close() 
        
        print(f'plik wyjsciowy: wsp_xy2000_out.txt')
        
    elif wybrana_funkcja == '--uklad1992':
        with open(plik, 'r') as f:
            lines = f.readlines()
            coord_lines = lines
        
            coords_xy1992 = []
            for coord_line in coord_lines:
                coord_line = coord_line.strip('\n')
                p_str, l_str, h_str = coord_line.split(',')
                p, l, h = float(p_str),float(l_str), float(h_str)
                x, y = geo.uklad1992(p, l, h)
                coords_xy1992.append([x, y])
        
        file_out = 'wsp_xy1992_out.txt'
        f1 = open(file_out, 'w')
        
        for xy in coords_xy1992:
            s = f'{xy[0]:.5f},{xy[1]:.5f} \n'
            f1.write(s)
            
        f1.close()
        
        print(f'plik wyjsciowy: wsp_xy1992_out.txt')
        
    else:
        raise NameError('nie ma takiej funkcji')
        
         












