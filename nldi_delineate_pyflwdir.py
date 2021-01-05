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
from pyproj import Geod
from shapely.ops import transform, split, snap
import shapely.geometry
from shapely.geometry import shape, mapping, Point, GeometryCollection
import json
import time
import numpy as np
import geopandas
from osgeo import ogr, osr, gdal
from functools import partial


#arguments
NLDI_URL = 'https://labs.waterdata.usgs.gov/api/nldi/linked-data/comid/'
NLDI_GEOSERVER_URL = 'https://labs.waterdata.usgs.gov/geoserver/wmadata/ows'
NHDPLUS_FLOWLINES_QUERY_URL = 'https://hydro.nationalmap.gov/arcgis/rest/services/nhd/MapServer/6/query'
NHD_FLOWLINES_URL = 'https://labs.waterdata.usgs.gov/geoserver/wmadata/ows?service=wfs&version=1.0.0&request=GetFeature&typeName=wmadata%3Anhdflowline_network&maxFeatures=500&outputFormat=application%2Fjson&srsName=EPSG%3A4326&CQL_FILTER=comid%3D'
ROOT_PATH = 'C:/Users/ahopkins/nldi/ACWI-SSWD/'
OUT_PATH = ROOT_PATH + 'nldi-splitCatchment/data/'
#IN_FDR = ROOT_PATH + 'nldi-splitCatchment/data/NHDPlusMA/NHDPlus02/NHDPlusFdrFac02b/fdr'
IN_FDR_COG = '/vsicurl/https://prod-is-usgs-sb-prod-publish.s3.amazonaws.com/5fe0d98dd34e30b9123eedb0/fdr.tif'

class Watershed:
    """Define inputs and outputs for the main Watershed class"""

    def __init__(self, x=None, y=None, fullBasin=True, export=True):

        self.x = x
        self.y = y
        self.fullBasin = fullBasin
        self.export = export
        self.catchmentIdentifier = None
        self.flowlines = None
        self.flw = None
        self.xy = None

        #geoms
        self.catchmentGeom = None
        self.splitCatchmentGeom = None
        self.upstreamBasinGeom = None
        self.mergedCatchmentGeom = None 
        self.intersectionPointGeom = None
        self.downstreamPathGeom = None   
        self.nhdFlowlinesGeom = None

        #outputs
        self.catchment = None
        self.splitCatchment = None
        self.upstreamBasin = None
        self.mergedCatchment = None
        self.intersectionPoint = None
        self.downstreamPath = None
        self.nhdFlowlines = None
        self.streamInfo = None

        #create transform
        self.transformToRaster = None
        self.transformToWGS84 = None
        self.transformToEqualDistance = None

        #kick off
        self.run()

    def serialize(self):
        return {
            'catchment': self.catchment,
            'splitCatchment': self.splitCatchment, 
            'upstreamBasin': self.upstreamBasin,
            'mergedCatchment': self.mergedCatchment,
            'intersectionPoint': self.intersectionPoint,
            'downstreamPath': self.downstreamPath,
            'nhdFlowlines': self.nhdFlowlines,
            'streamInfo': self.streamInfo
        }

## helper functions
    def geom_to_geojson(self, geom, name, write_output=False):
        """Return a geojson from an OGR geom object"""

        from pyproj import Geod
        from shapely import wkt

        # specify a named ellipsoid
        geod = Geod(ellps="WGS84")

        poly = geom
        #area = abs(geod.geometry_area_perimeter(poly)[0])

        #print(name, 'area: {:12.3f} mi^2'.format(area*0.00000038610))

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
        self.downstreamPathGeom, self.intersectionPointGeom, self.nhdFlowlinesGeom, self.streamInfo = self.get_downstreamPath(self.catchmentGeom, self.catchmentIdentifier, self.flw ,self.xy)

        #outputs
        self.catchment = self.geom_to_geojson(self.catchmentGeom, 'catchment', self.export)
        self.splitCatchment = self.geom_to_geojson(self.splitCatchmentGeom, 'splitCatchment', self.export)

        if self.fullBasin:
            self.upstreamBasinGeom = self.get_upstream_basin(self.catchmentIdentifier)
            self.mergedCatchmentGeom = self.merge_geometry(self.catchmentGeom, self.splitCatchmentGeom, self.upstreamBasinGeom)
            self.upstreamBasin = self.geom_to_geojson(self.upstreamBasinGeom, 'upstreamBasin', self.export)
            self.mergedCatchment = self.geom_to_geojson(self.mergedCatchmentGeom, 'mergedCatchment', self.export)

        self.upstreamBasin = self.geom_to_geojson(self.upstreamBasinGeom, 'upstreamBasin')
        self.mergedCatchment = self.geom_to_geojson(self.mergedCatchmentGeom, 'mergedCatchment')
        self.intersectionPoint = self.geom_to_geojson(self.intersectionPointGeom, 'intersectionPoint')
        self.downstreamPath = self.geom_to_geojson(self.downstreamPathGeom, 'downstreamPath')
        self.nhdFlowlines = self.geom_to_geojson(self.nhdFlowlinesGeom, 'nhdFlowlines')

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

        print('got local catchment:', catchmentIdentifier)
        return catchmentIdentifier, catchmentGeom

    def get_local_flowlines(self, catchmentIdentifier):
        # Request NDH Flowlines

        cql_filter = "comid=%s" % (catchmentIdentifier) 

        payload = {
            'service': 'wfs', 
            'version': '1.0.0', 
            'request': 'GetFeature', 
            'typeName': 'wmadata:nhdflowline_network', 
            'maxFeatures': '500',
            'outputFormat': 'application/json',
            'srsName': 'EPSG:4326',
            'CQL_FILTER': cql_filter
        }

        #request  flowline geometry from point in polygon query from NLDI geoserver
        r = requests.get(NLDI_GEOSERVER_URL, params=payload)

        print('request url: ', r.url)
        flowlines = r.json()

        print('got local flowlines')
        return flowlines

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
        with rasterio.open(IN_FDR_COG, 'r') as ds:
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
        self.flw = flw

        point_geom = Point(self.x,self.y)
        print('original point:',point_geom)

        projected_point = transform(self.transformToRaster, point_geom)
        print('projected point:',projected_point)

        xy = projected_point.coords[:][0]
        self.xy = xy

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


    def get_downstreamPath(self, catchment_geom, catchmentIdentifier, flw, xy):
        """Use catchment bounding box to clip NHD Plus v2 flow direction raster, and trace a flowpath from X,Y"""

        # Request the NHD Flowline from the catchment
        r = requests.get(NHD_FLOWLINES_URL + catchmentIdentifier)
        flowlines = r.json()
        print('got flowlines   ', flowlines)

        # Convert the flowline to a geometry colelction to be exported
        nhdGeom = flowlines['features'][0]['geometry']
        nhdFlowlines = GeometryCollection([shape(nhdGeom)])

        streamname = flowlines['features'][0]['properties']['gnis_name']
        if streamname == ' ':
            streamname = 'none'

        streamInfo = {'gnis_name' : streamname,
                      'comid' : flowlines['features'][0]['properties']['comid'],
                      'lengthkm' : flowlines['features'][0]['properties']['lengthkm']}

        # Convert the flowlines to a geopandas dataframe
        dfNHD = geopandas.GeoDataFrame.from_features(flowlines, crs="EPSG:4326")

        # Project the flowlines to the same crs as the flw raster
        projectedNHD = self.transform_geom(self.transformToRaster, dfNHD.geometry[0][0])
           
        # Convert the flowline coordinates to a format that can be iterated
        line = list(projectedNHD.coords)
        print('created list of nhd coords  ')

        # loop thru the flowline coordinates, grab the xy coordinantes and put them in separate lists.
        # Use these lists in the index function of pyflwdir to grap the ids of the cell in which these points fall
        xlist = []
        ylist = []
        cellIndexList = []
        for i in line:
            if i == line[0]:
                xlist = []
            if i != 0:
                xlist = (i[0])
                ylist = (i[1])
            cellIndex = flw.index(xlist, ylist)
            cellIndexList.append(cellIndex)
        print('nhd converted to raster  ')

        # create mask from in the same of the flw raster
        nhdMask = np.zeros(flw.shape, dtype=np.bool)

        # Set the flowline cells to true
        nhdMask.flat[cellIndexList] = True

        # trace downstream
        path, dist = flw.path( 
            xy=xy,
            mask=nhdMask
            )
        print('traced downstreampath   ')

        # get points on downstreamPath 
        pathPoints = flw.xy(path)

        # loop thru the downstream path points and create a dict of coords
        lastPointID = pathPoints[0][0].size - 1
        i = 0
        coordlist = {'type': 'LineString', 'coordinates': []}
        while i <= lastPointID:
            x = pathPoints[0][0][i]
            y = pathPoints[1][0][i]
            coordlist['coordinates'].append([x,y])
            i+=1

        # Convert the dict of coords to ogr geom
        pathGeom = GeometryCollection([shape(coordlist)])

        # Project the ogr geom to WGS84
        projectedPathGeom = transform(self.transformToWGS84, pathGeom)

        # Snap downstreamPath points to the flowline within a buffer
        snapPath = snap(projectedPathGeom[0], nhdFlowlines[0], .00045)
        
        # Conver snapPath to a geometry collection
        snapPath = GeometryCollection([snapPath])
        
        # Grap all the points of intersection between the downstreamPath and the flowline
        intersectionpoints = snapPath.intersection(nhdFlowlines[0])
        print('intersectionpoints', type(intersectionpoints), intersectionpoints)

        # Use if statements to filter the intersecting points bt geometry type. The downstream path 
        # will then be split by each point in the intersectionpoints geom. 
        if type(intersectionpoints) == shapely.geometry.multipoint.MultiPoint: 
            for i in intersectionpoints:
                splitPoint = snap(Point(i.coords), snapPath, .0002)
                snapPath = split(snapPath[0], splitPoint)
        if type(intersectionpoints) == shapely.geometry.linestring.LineString: 
            for i in intersectionpoints.coords:
                splitPoint = snap(Point(i), snapPath, .0002)
                snapPath = split(snapPath[0], splitPoint)
        if type(intersectionpoints) == shapely.geometry.point.Point:
            splitPoint = snap(intersectionpoints, snapPath, .0002)
            snapPath = split(snapPath[0], splitPoint)
        if type(intersectionpoints) == shapely.geometry.multilinestring.MultiLineString:
            for i in intersectionpoints:
                for j in i.coords:
                    splitPoint = snap(Point(j), snapPath, .0002)
                    snapPath = split(snapPath[0], splitPoint)
        if type(intersectionpoints) == shapely.geometry.collection.GeometryCollection:
            for i in intersectionpoints:
                for j in i.coords:
                    splitPoint = snap(Point(j), snapPath, .0002)
                    snapPath = split(snapPath[0], splitPoint)
           
        # The first linestring in the snapPath geometry collection in the downstreamPath
        downstreamPath = snapPath[0]

        # The last point on the downstreamPath is the intersectionPoint
        coordID = len(downstreamPath.coords) - 1
        intersectionPoint = Point(downstreamPath.coords[coordID])

        # Get lengths of NHD Flowline before and after the intersectionpoint
        geod = Geod(ellps="WGS84")
        #snapIntersectionPoint = snap(intersectionPoint, nhdFlowlines[0], .0004)
        print('Does the intersectionPoint intersect the nhdFlowlines? ', nhdFlowlines[0].intersects(intersectionPoint))
        NHDFlowlinesCut = split(nhdFlowlines[0], intersectionPoint)
        print('NHDFlowlinesCut', type(NHDFlowlinesCut), NHDFlowlinesCut)
        streamInfo['downstreamPathDist'] = round(geod.geometry_length(downstreamPath), 2)
        streamInfo['flowlineLength'] = round(geod.geometry_length(nhdFlowlines[0]) / 1000, 3)
        try:
             NHDFlowlinesCut[1]
        except: 
            streamInfo['distFromStart'] = round(geod.geometry_length(NHDFlowlinesCut), 2)
            print('except')
        else:
            streamInfo['distFromStart'] = round(geod.geometry_length(NHDFlowlinesCut[0]), 2)
            streamInfo['distToEnd'] = round(geod.geometry_length(NHDFlowlinesCut[1]), 2)
            print('else')
        
        print('intersectionPoint', type(intersectionPoint), intersectionPoint)
        print('downstreamPath', type(downstreamPath), downstreamPath)
        print('streamInfo: ', type(streamInfo), streamInfo)
        print('nhdFlowlines: ', type(nhdFlowlines), nhdFlowlines[0].length)

        return downstreamPath, intersectionPoint, nhdFlowlines, streamInfo

if __name__=='__main__':

    timeBefore = time.perf_counter()  

    #test site
    #point = (-73.82705, 43.29139)
    point = (-73.72435569763185, 43.17261895666325)

    #start main program
    delineation = Watershed(point[0],point[1], True, False)

    timeAfter = time.perf_counter() 
    totalTime = timeAfter - timeBefore
    print("Total Time:",totalTime)