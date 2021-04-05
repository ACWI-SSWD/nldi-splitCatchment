import requests
import json
from shapely.geometry import Point, LineString, Polygon

##### Resolution types and their respective IDs for the Rest Service
res_types = {'res_1m': 19, 'res_3m': 19, 'res_5m': 20, 'res_10m': 21, 'res_30m': 22, 'res_60m': 23}

###### Create a bounding box from any geo type ##########
# 'width' is in meters. It is the width of the buffer to place around the input geometry
def make_bbox(shape_type, coords, width):
    if shape_type == 'point':
        shape = Point(coords)
    
    if shape_type == 'line':
        shape = LineString(coords)
    
    if shape_type == 'polygon':
        shape = Polygon(coords)
        
    converted_width = convert_width(width)
    buff_shape = shape.buffer(converted_width)
    return buff_shape.bounds

##### Convert the width of the buffer from meters to 'degree' ######
# This is NOT a precise conversion, just a quick overestimation
# Maybe fix this later
def convert_width(width):
    factor = 1/70000
    dist = factor*width
    return dist

###### Get request from arcgis Rest Services ###########
def get_dem(bbox, res_type):
    miny = str(bbox[0])
    minx = str(bbox[1])
    maxy = str(bbox[2])
    maxx = str(bbox[3])
    res_id = res_types[res_type]

    url = f'https://index.nationalmap.gov/arcgis/rest/services/3DEPElevationIndex/MapServer/{res_id}/query'
    payload = {
        "where":"",
        "text":"",
        "objectIds":"",
        "time":"",
        "geometry":"{xmin:\""+minx+"\",ymin:\""+miny+"\",xmax:\""+maxx+"\",ymax:\""+maxy+"\",spatialReference:{wkid:4326}}",
        "geometryType":"esriGeometryEnvelope",
        "inSR":"EPSG:4326",
        "spatialRel":"esriSpatialRelWithin",
        "relationParam":"",
        "outFields":"",
        "returnGeometry":"false",
        "returnTrueCurves":"false",
        "maxAllowableOffset":"",
        "geometryPrecision":"",
        "outSR":"EPSG:4326",
        "having":"",
        "returnIdsOnly":"false",
        "returnCountOnly":"false",
        "orderByFields":"",
        "groupByFieldsForStatistics":"",
        "outStatistics":"",
        "returnZ":"false",
        "returnM":"false",
        "gdbVersion":"",
        "historicMoment":"",
        "returnDistinctValues":"false",
        "resultOffset":"",
        "resultRecordCount":"",
        "queryByDistance":"",
        "returnExtentOnly":"false",
        "datumTransformation":"",
        "parameterValues":"",
        "rangeValues":"",
        "quantizationParameters":"",
        "featureEncoding":"esriDefault",
        "f":"geojson"
        }
    
    r = requests.get(url, params=payload)
    resp = r.json()

    #### If the Rest Services has a 200 response
    if 'features' in resp:
        if resp['features'] != []:  # If the features are not empty, then the DEM exist
            return True
        if resp['features'] == []: # If features are empty, then no DEM
            return  False
    if 'error' in resp: # If 'error', then there was an unsuccessful request
        return  'Error!'

###### The function to loop thru all resolutions and submit queries ############
def query_dems(shape_type, coords, width=100):
    resp = {}  # Create an empty dictionary for the response
    bbox = make_bbox(shape_type, coords, width)  # Make the bbox
    for res_type in res_types:   # Loop thru all resolutions and submit a query for each one
        resp[res_type] = (get_dem(bbox, res_type))  # Add the response to the dictionary

    # print(resp)
    return(resp)

###### Tests #########
# query_dems('point', [(39,-96)])
# query_dems('line', [(39,-96.0001), (39,-96)])
# query_dems('polygon', [(39,-96.0008), (39,-96), (39.0008, -96)])