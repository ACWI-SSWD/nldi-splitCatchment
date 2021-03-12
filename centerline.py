import geopandas as gpd
import numpy as np
from tspline import tspline

class centerline:
    """
    class to take a centerline as a shapefile and create a tension spline to be used
    with inherited class curvilinear grid
    """

    def __init__(self, center_shp, nx, tension):
        self.center_shp = center_shp
        self.gpdata = center_shp
        self.x = np.empty(0, dtype=np.double)
        self.y = np.empty(0, dtype=np.double)
        self.tension = 0.5
        self.numclpts = 0
        self.numInterpPts = nx
        self.si = np.empty(0, dtype=np.double)
        self.temp = np.empty(0, dtype=np.double)
        self.yp = np.empty(0, dtype=np.double)
        self.xo_interp = np.empty(0, dtype=np.double)
        self.yo_interp = np.empty(0, dtype=np.double)
        self.phi_interp = np.empty(0, dtype=np.double)
        self.r_interp = np.empty(0, dtype=np.double)
        self.stot = np.empty(0, dtype=np.double)
        self.sout = np.empty(0, dtype=np.double)
        self.__initialize()
        self.__getspline(nx, self.tension)

    def __initialize(self):
        tx = []
        ty = []
        for index, row in self.gpdata.iterrows():
            for pt in list(row['geometry'].coords):
                # print(pt, type(pt))
                tx.append(pt[0])
                ty.append(pt[1])
        self.x = np.array(tx)
        self.y = np.array(ty)
        self.numclpts = self.x.size

    # def description(self):
    #     return "{} used the shapfile {}".format("centerline", self.center_shp)

    # def getCLShapeFile(self):
    #     return self.center_shp

    def getlength(self):
        return

    def getpoints(self):
        return self.x, self.y

    def getinterppts(self):
        return self.xo_interp, self.yo_interp

    def getinterppts_dyn(self, numinterppts, tension):
        print(tension)
        self.__getspline(numinterppts, tension)
        return self.xo_interp, self.yo_interp

    def getphiinterp(self, index):
        return self.phi_interp[index]

    def __getspline(self, numinterppts, tension):
        self.numInterpPts = numinterppts
        self.tension = tension
        self.xo_interp = np.zeros(self.numInterpPts)
        self.yo_interp = np.zeros(self.numInterpPts)
        self.phi_interp = np.zeros(self.numInterpPts)
        self.r_interp = np.zeros(self.numInterpPts)
        self.stot = np.zeros(self.numInterpPts)
        self.sout = np.zeros(self.numInterpPts)

        self.si = np.zeros(self.numclpts)
        self.temp = np.zeros(self.numclpts)
        self.yp = np.zeros(self.numclpts)

        for i in range(0, self.numclpts):
            if i == 0:
                self.si[i] = 0.0
            else:
                ds = np.power((np.power((self.x[i] - self.x[i - 1]), 2) +
                               np.power((self.y[i] - self.y[i - 1]), 2)), 0.5)
                self.si[i] = self.si[i-1] + ds

        self.stot = self.si[self.numclpts - 1]
        for i in range(0, self.numInterpPts):
            self.sout[i] = i * self.stot / (self.numInterpPts - 1)
            # if i == 0:
            #     nsInc = self.stot / (self.numInterpPts - 1)
        if self.numclpts < 3:
            self.yp = np.zeros(3)
            self.temp = np.zeros(3)
            sitmp = np.zeros(3)
            sitmp = np.append(sitmp, [self.si[0]])
            sitmp = np.append(sitmp, [self.si[1] / 2])
            sitmp = np.append(sitmp, self.si[1])

            txctmp = np.zeros(3)
            txctmp = np.append(txctmp, self.x[0])
            txctmp = np.append(txctmp, self.x[0] + (self.x[1] - self.x[0]) / 2.)
            txctmp = np.append(txctmp, self.x[1])

            tyctmp = np.zeros(3)
            tyctmp = np.append(tyctmp, self.y[0])
            tyctmp = np.append(tyctmp, self.y[0] + (self.y[1] - self.y[0]) / 2.)
            tyctmp = np.append(tyctmp, self.y[1])

            tspline(sitmp, txctmp, 3, self.sout, self.xo_interp,
                    self.numInterpPts, self.tension, self.yp, self.temp)
            tspline(sitmp, tyctmp, 3, self.sout, self.yo_interp,
                    self.numInterpPts, self.tension, self.yp, self.temp)
        else:
            tspline(self.si, self.x, self.numclpts, self.sout, self.xo_interp,
                    self.numInterpPts, self.tension, self.yp, self.temp)
            tspline(self.si, self.y, self.numclpts, self.sout, self.yo_interp,
                    self.numInterpPts, self.tension, self.yp, self.temp)
        self.__calcCurvature()

    def __calcCurvature(self):
        for i in range(1, self.numInterpPts):
            dx = self.xo_interp[i] - self.xo_interp[i - 1]
            dy = self.yo_interp[i] - self.yo_interp[i - 1]
            if dx == 0.:
                if dy > 0.:
                    self.phi_interp[i] = np.pi/2.
                elif dy < 0.:
                    self.phi_interp[i] = -1.*np.pi/2.
            else:
                self.phi_interp[i] = np.arctan2(dy, dx)
        self.phi_interp[0] = (2. * self.phi_interp[1]) - self.phi_interp[2]
        scals = self.stot/(self.numInterpPts - 1)
        for i in range(1, self.numInterpPts):
            dx = self.xo_interp[i] - self.xo_interp[i - 1]
            dphi = np.fabs(self.phi_interp[i]) - np.fabs(self.phi_interp[i - 1])
            if dphi <= 0.0001:
                self.r_interp[i] = 100000000.
            else:
                self.r_interp[i] = scals/dphi