import numpy as np
from wzory_gw import *

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
       
