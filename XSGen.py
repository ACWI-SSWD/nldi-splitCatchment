from centerline import centerline
import numpy as np
from shapely.geometry import LineString, Point
import geopandas as gpd
import pandas as pd


class XSGen:
    def __init__(self, point, cl_geom, ny, width) -> None:
        self.cl_geom = cl_geom

        self.point = point
        self.tension = 0.5
        self.width = width
        self.cl_length = self.cl_geom.geometry[0].length
        self.ny = ny
        if self.cl_length > 20.0:
            self.nx = np.int(self.cl_length / 10)
        else:
            self.nx = np.int(self.cl_length / 1)
        self.cl = centerline(cl_geom, self.nx, self.tension)
        self.x = np.zeros(self.ny, dtype=np.double)
        self.y = np.zeros(self.ny, dtype=np.double)
        self._buildxs()

    def _get_perp_index(self, clx, cly):
        mind = 1e6
        id = -1
        for index, p in enumerate(zip(clx, cly)):
            dist = np.sqrt(np.power(self.point.geometry[0].x - p[0], 2) + np.power(self.point.geometry[0].y - p[1], 2))
            if dist < mind:
                mind = dist
                id = index

        return id

    def _buildxs(self):
        clx, cly = self.cl.getinterppts()
        delt = np.double(self.width / (self.ny - 1))
        nm = np.int((self.ny + 1) / 2)
        index = self._get_perp_index(clx, cly)

        for id, j in enumerate(range(0, self.ny)):
            self.x[id] = clx[index] + delt * (nm - j - 1) * np.sin(self.cl.getphiinterp(index))
            self.y[id] = cly[index] - delt * (nm - j - 1) * np.cos(self.cl.getphiinterp(index))

    def get_xs(self):
        points = gpd.GeoSeries(map(Point, zip(self.x, self.y)))
        ls = LineString((points.to_list()))
        d = {0: {'name': 'cross-section', 'geometry': ls}}
        df = pd.DataFrame.from_dict(d, orient='index')
        gdf = gpd.GeoDataFrame(df, geometry=df.geometry)
        return gdf


    def get_xs_points(self):
        return self.x, self.y


    def get_strm_seg_spine(self):
        x,y = self.cl.getinterppts()
        points = gpd.GeoSeries(map(Point, zip(x,y)))
        ls = LineString((points.to_list()))
        d = {0: {'name': 'strm_seg_spline', 'geometry': ls}}
        df = pd.DataFrame.from_dict(d, orient='index')
        gdf = gpd.GeoDataFrame(df, geometry=df.geometry)
        return gdf