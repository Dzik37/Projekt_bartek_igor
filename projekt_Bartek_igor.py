import numpy as np
from math import *
class Transformacje:
    def __init__(self, model: str = "wgs84"):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            e2 - mimośród^2
        """
        if model == "wgs84":
            self.a = 6378137.0 
            self.e2 = 0.006694379990141124
        elif model == "grs80":
            self.a = 6378137.0
            self.e2 = 0.00669438002290
        elif model == "krasowski":
            self.a = 6378245.000 
            self.e2 = 0.00669342162296
        else:
            raise NotImplementedError(f"{model} model not implemented")
       
    def xyz2blh(self,X, Y, Z, output = 'f,l,h'):
        l = np.arctan2(Y, X)
        p = np.sqrt(X**2 + Y**2)
        f = np.arctan(Z / (p * (1 -self.e2)))
        while True:
        
            N = self.a / np.sqrt(1 - self.e2 * np.sin(f)**2)
            h = p / np.cos(f) - N
            fs = f
            f = np.arctan(Z / (p * (1 - (self.e2 * (N / (N + h))))))
            if np.abs(fs - f) < (0.000001/206265):
                break
            if output == "dec_degree":
                return f,l,h 
            elif output == "dms":
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
                    return(dms(f),dms(l),h)
        
        else:
            raise NotImplementedError(f"{output} - output format not defined")
if __name__ == "__main__":
    # utworzenie obiektu
    geo = Transformacje(model = "wgs84")
    # dane XYZ geocentryczne
    X = 3664940.500; Y = 1409153.590; Z = 5009571.170
    phi, lam, h = geo.xyz2blh(X, Y, Z)
    print(phi, lam, h)