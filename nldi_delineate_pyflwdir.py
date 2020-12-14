## NLDI split catchment delineation script

# -----------------------------------------------------
# Martyn Smith USGS
# 10/29/2020
# NLDI Delineation script
# -----------------------------------------------------

# list of required python packages:
# gdal, pysheds, requests

###### CONDA CREATE ENVIRONMENT COMMAND
#conda create -n delineate python=3.6.8 gdal pysheds requests
###### CONDA CREATE ENVIRONMENT COMMAND

import requests
import rasterio
import rasterio.mask
import pyflwdir
import pyproj
from shapely.ops import transform
from shapely.geometry import shape, mapping, Point, Polygon, GeometryCollection
import json
import time
import numpy as np

#arguments
NLDI_URL = 'https://labs.waterdata.usgs.gov/api/nldi/linked-data/comid/'
NLDI_GEOSERVER_URL = 'https://labs.waterdata.usgs.gov/geoserver/wmadata/ows'
NHDPLUS_FLOWLINES_QUERY_URL = 'https://hydro.nationalmap.gov/arcgis/rest/services/nhd/MapServer/6/query'
OUT_PATH = 'C:/NYBackup/GitHub/nldi-splitCatchment/data/'
IN_FDR = 'C:/NYBackup/GitHub/nldi-splitCatchment/data/nhdplus/NHDPlusMA/NHDPlus02/NHDPlusFdrFac02b/fdr'

class Watershed:
    """Define inputs and outputs for the main Watershed class"""

    def __init__(self, x=None, y=None):

        self.x = x
        self.y = y
        self.catchmentIdentifier = None

        #geoms
        self.catchmentGeom = None
        self.splitCatchmentGeom = None
        self.upstreamBasinGeom = None
        self.mergedCatchmentGeom = None    

        #outputs
        self.catchment = None
        self.splitCatchment = None
        self.upstreamBasin = None
        self.mergedCatchment = None

        #create transform
        self.transformToRaster = None
        self.transformToWGS84 = None

        #kick off
        self.run()

    def serialize(self):
        return {
            'catchment': self.catchment,
            'splitCatchment': self.splitCatchment, 
            'upstreamBasin': self.upstreamBasin,
            'mergedCatchment': self.mergedCatchment
        }

## helper functions
    def geom_to_geojson(self, geom, name, simplify_tolerance=10, write_output=False):
        """Return a geojson from an OGR geom object"""

        geojson_dict = mapping(geom)

        if write_output:
            f = open(OUT_PATH + name + '.geojson','w')
            f.write(json.dumps(geojson_dict))
            f.close()
            print('Exported geojson:', name)
        
        return geojson_dict

## main functions
    def run(self):
        self.catchmentIdentifier, self.catchmentGeom = self.get_local_catchment(self.x,self.y)
        self.splitCatchmentGeom = self.split_catchment(self.catchmentGeom, self.x, self.y)

        self.upstreamBasinGeom = self.get_upstream_basin(self.catchmentIdentifier)
        self.mergedCatchmentGeom = self.merge_geometry(self.catchmentGeom, self.splitCatchmentGeom, self.upstreamBasinGeom)

        #outputs
        self.catchment = self.geom_to_geojson(self.catchmentGeom, 'catchment')
        self.splitCatchment = self.geom_to_geojson(self.splitCatchmentGeom, 'splitCatchment')

        #print('test', self.splitCatchment)
        self.upstreamBasin = self.geom_to_geojson(self.upstreamBasinGeom, 'upstreamBasin')
        self.mergedCatchment = self.geom_to_geojson(self.mergedCatchmentGeom, 'mergedCatchment')

    def transform_geom(self, proj, geom):
        """Transform geometry"""

        projected_geom = transform(proj, geom)

        return projected_geom

    def get_local_catchment(self, x, y):
        """Perform point in polygon query to NLDI geoserver to get local catchment geometry"""

        print('requesting local catchment...')

        wkt_point = "POINT(%f %f)" %  (x , y)
        cql_filter = "INTERSECTS(the_geom, %s)" % (wkt_point)

        payload = {
            'service': 'wfs', 
            'version': '1.0.0', 
            'request': 'GetFeature', 
            'typeName': 'wmadata:catchmentsp', 
            'outputFormat': 'application/json',
            'srsName': 'EPSG:4326',
            'CQL_FILTER': cql_filter
        }

        #request catchment geometry from point in polygon query from NLDI geoserver
        # https://labs.waterdata.usgs.gov/geoserver/wmadata/ows?service=wfs&version=1.0.0&request=GetFeature&typeName=wmadata%3Acatchmentsp&outputFormat=application%2Fjson&srsName=EPSG%3A4326&CQL_FILTER=INTERSECTS%28the_geom%2C+POINT%28-73.745860+44.006830%29%29
        r = requests.get(NLDI_GEOSERVER_URL, params=payload)
        resp = r.json()

        #get catchment id
        catchmentIdentifier = json.dumps(resp['features'][0]['properties']['featureid'])

        #get main catchment geometry polygon
        features = resp['features']
        catchmentGeom = GeometryCollection([shape(feature["geometry"]).buffer(0) for feature in features])

        print('got local catchment')
        return catchmentIdentifier, catchmentGeom

    def get_upstream_basin(self, catchmentIdentifier):
        """Use local catchment identifier to get upstream basin geometry from NLDI"""

        print('getting upstream basin...')

        #request upstream basin
        payload = {'f': 'json', 'simplified': 'false'}
        
        #request upstream basin from NLDI using comid of catchment point is in
        r = requests.get(NLDI_URL + catchmentIdentifier + '/basin', params=payload)

        #print('upstream basin', r.text)
        resp = r.json()

        #convert geojson to ogr geom
        features = resp['features']
        upstreamBasinGeom = GeometryCollection([shape(feature["geometry"]).buffer(0) for feature in features])

        print('finished getting upstream basin')
        return upstreamBasinGeom

    def merge_geometry(self, catchment, splitCatchment, upstreamBasin):
        """Attempt at merging geometries"""

        #Test if point is on a flowline we have an upstream basin and need to do some geometry merging
        if 1==1:

            print('merging geometries...')
            d = 0.00045
            #d2 = 0.00015 # distance
            cf = 1.3  # cofactor

            splitCatchment = splitCatchment.simplify(d)

            diff = catchment.difference(splitCatchment).buffer(-d).buffer(d*cf).simplify(d)
            mergedCatchmentGeom = upstreamBasin.difference(diff).buffer(-d).buffer(d*cf).simplify(d)

            print('finished merging geometries')
            #write out
            return mergedCatchmentGeom

        #otherwise, we can just return the split catchment
        else:
            mergedCatchmentGeom = splitCatchment

        return mergedCatchmentGeom

    def split_catchment(self, catchment_geom, x, y): 
        """Use catchment bounding box to clip NHD Plus v2 flow direction raster, and product split catchment delienation from X,Y"""

        print('start clip raster')
        with rasterio.open(IN_FDR, 'r') as ds:
            #get raster crs
            dest_crs = ds.crs

            #create wgs84 crs
            wgs84 = pyproj.CRS('EPSG:4326')

            #check to see if raster is already wgs84
            latlon = dest_crs == wgs84
            
            self.transformToRaster = pyproj.Transformer.from_crs(wgs84, dest_crs, always_xy=True).transform
            self.transformToWGS84 = pyproj.Transformer.from_crs(dest_crs, wgs84, always_xy=True).transform

            #transform catchment geometry to use for clip
            projected_catchment_geom = self.transform_geom(self.transformToRaster, catchment_geom)

            #clip input fd
            flwdir, flwdir_transform = rasterio.mask.mask(ds, projected_catchment_geom, crop=True)
            print('finish clip raster')
            
        #import clipped fdr into pyflwdir
        flw = pyflwdir.from_array(flwdir[0], ftype='d8', transform=flwdir_transform, latlon=latlon)

        point_geom = Point(self.x,self.y)
        print('original point:',point_geom)

        projected_point = transform(self.transformToRaster, point_geom)
        print('projected point:',projected_point)

        xy = projected_point.coords[:][0]

        #used for snapping click point
        stream_order = flw.stream_order()

        print('start split catchment...')

        # delineate subbasins
        subbasins = flw.basins(xy=xy, streams=stream_order>4)

        #convert subbasins from uint32
        subbasins = subbasins.astype(np.int32)

        #convert raster to features
        mask = subbasins != 0
        polys = rasterio.features.shapes(subbasins, transform=flwdir_transform, mask=mask)

        #just get one we want [not sure why we need to grab this]
        poly = next(polys)

        #project back to wgs84
        split_geom = transform(self.transformToWGS84, shape(poly[0]))

        print('finish split catchment...')
        return split_geom       

if __name__=='__main__':

    timeBefore = time.perf_counter()  

    #test site
    point = (-73.82705, 43.29139)
    #point = (-73.72435569763185, 43.17261895666325)

    #start main program
    delineation = Watershed(point[0],point[1])

    timeAfter = time.perf_counter() 
    totalTime = timeAfter - timeBefore
    print("Total Time:",totalTime)