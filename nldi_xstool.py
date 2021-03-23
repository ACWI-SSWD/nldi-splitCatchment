"""Main module."""
# import aiodns
# print('ok: aiodns')
# import aiohttp
# print('ok: aiohttp')
# # import brotlipy
# # print('ok: brotlipy')
# import chardet
# print('ok: chardet')
# import cytoolz
# print('ok: cytoolz')
# import defusedxml
# print('ok: defusedxml')
# import geopandas
# print('ok: geopandas')
# # import nest-asyncio
# # print('ok: nest-asyncio')
# import netCDF4
# print('ok: netCDF4')
# import orjson
# print('ok: orjson')
# import owslib
# print('ok: owslib')
# import pygeoogc
# print('ok: pygeoogc')
# import pygeoutils
# print('ok: pygeoutils')
# import pyproj
# print('ok: pyproj')
# import rasterio
# print('ok: rasterio')
# import requests
# print('ok: requests')
# import setuptools
# print('ok: setuptools')
# import shapely
# print('ok: shapely')
# import simplejson
# print('ok: simplejson')
# import xarray
# print('ok: xarray')


from XSGen import XSGen
import requests
import json
import py3dep
# from pynhd import NLDI
import xarray as xr
from matplotlib import pyplot as plt
from shapely.geometry import Point
import geopandas as gpd
import pandas as pd
#import os.path as path

class HPoint(Point):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __hash__(self):
       return hash(tuple(self.coords))


def dataframe_to_geodataframe(df):
    geometry = [HPoint(xy) for xy in zip(df.x, df.y)]
    df = df.drop(['x','y'], axis=1)
    gdf = gpd.GeoDataFrame(df, geometry=geometry)
    return gdf

def get_strm_seg(comid):
    url = f'https://labs.waterdata.usgs.gov/api/nldi/linked-data/comid/{comid}?f=json'
    print('url ', url)
    r = requests.get(url)
    r = r.json()
    strmseg = gpd.GeoDataFrame.from_features(r, crs="EPSG:4326")
    print("response: ", type(strmseg), strmseg)
    return strmseg

def getXSAtPoint(point, numpoints, width, file=None):
    tpoint = f'POINT({point[1]} {point[0]})'
    df = pd.DataFrame({'pointofinterest':['this'],
                        'Lat':[point[0]],
                        'Lon':[point[1]]})
    gpd_pt = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.Lon, df.Lat))
    gpd_pt.set_crs(epsg=4326, inplace=True)
    gpd_pt.to_crs(epsg=3857, inplace=True)
    comid = getCIDFromLatLon(point)
    print(f'comid = {comid}')
    # strm_seg = NLDI().getfeature_byid("comid", comid)#.to_crs('epsg:3857')  # "3561878" , basin=False
    # strm_seg = gpd.GeoDataFrame.from_features({"type":"FeatureCollection","features":[{"type":"Feature","geometry":{"type":"LineString","coordinates":[[-103.802528001368,40.2687375992537],[-103.801507800817,40.2685499936342],[-103.797674000263,40.2685364931822],[-103.795592702925,40.2677875980735]]},"properties":{"source":"comid","sourceName":"NHDPlus comid","identifier":"3561878","name":"","uri":"","comid":"3561878","navigation":"https://labs.waterdata.usgs.gov/api/nldi/linked-data/comid/3561878/navigation"}}]}, crs="EPSG:4326").to_crs('epsg:3857')
    strm_seg = get_strm_seg(comid).to_crs('epsg:3857')
    print('strm_seg: ', strm_seg)
    xs = XSGen(point=gpd_pt, cl_geom=strm_seg, ny=100, width=1000)
    xs_line = xs.get_xs()
    # get topo polygon with buffer to ensure there is enough topography to interpolate xs line
    # With coarsest DEM (30m) 100. m should
    bb = xs_line.total_bounds - ((100., 100., -100., -100.))
    dem = py3dep.get_map("DEM", tuple(bb), resolution=10, geo_crs="EPSG:3857", crs="epsg:3857")  
    x,y = xs.get_xs_points()
    dsi = dem.interp(x=('z', x), y=('z', y))
    pdsi = dsi.to_dataframe()

    # gpdsi = gpd.GeoDataFrame(pdsi, gpd.points_from_xy(pdsi.x.values, pdsi.y.values))
    gpdsi = dataframe_to_geodataframe(pdsi)
    gpdsi.set_crs(epsg=3857, inplace=True)
    gpdsi.to_crs(epsg=4326, inplace=True)
    if(file):
        print("Output file: ", type(file), file)
        if not isinstance(file, str):
        # with open(file, "w") as f:
            file.write(gpdsi.to_json())
            file.close()
            return 0
        else:
            with open(file, "w") as f:
                f.write(gpdsi.to_json())
                f.close()
        # gpdsi.to_file(file, driver="GeoJSON")
            return 0
    else:
        return gpdsi

def latlonToPoint(lat, lon):
    return Point(lat, lon)

def getCIDFromLatLon(point):
    print(point)
    pt = latlonToPoint(point[1], point[0])
    location = pt.wkt
    location = f'POINT({point[1]} {point[0]})'
    baseURL = 'https://labs.waterdata.usgs.gov/api/nldi/linked-data/comid/position?f=json&coords='
    url = baseURL+location
    print(url)
    response = requests.get(url)
    jres = response.json()
    comid = jres['features'][0]['properties']['comid']
    return comid


