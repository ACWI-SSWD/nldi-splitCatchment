import requests
import json
from shapely.geometry import Point, LineString, Polygon

res_types = {'res_1m': 19, 'res_3m': 19, 'res_5m': 20, 'res_10m': 21, 'res_30m': 22, 'res_60m': 23}

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

def convert_width(width):

    factor = 1/70000
    dist = factor*width
    return dist

def query_dem(bbox, res_type):
    # print(bbox)
    miny = bbox[0]
    minx = bbox[1]
    maxy = bbox[2]
    maxx = bbox[3]
    # print(minx)
    # print(miny)
    # print(maxx)
    # print(maxy)

    res_id = res_types[res_type]

    url = f'https://index.nationalmap.gov/arcgis/rest/services/3DEPElevationIndex/MapServer/{res_id}/query?where=&text=&objectIds=&time=&geometry=%7B%22xmin%22%3A{minx}%2C%22ymin%22%3A{miny}%2C%22xmax%22%3A{maxx}%2C%22ymax%22%3A{maxy}%2C%22spatialReference%22%3A%7B%22wkid%22%3A4326%7D%7D&geometryType=esriGeometryEnvelope&inSR=EPSG%3A4326&spatialRel=esriSpatialRelWithin&relationParam=&outFields=&returnGeometry=false&returnTrueCurves=false&maxAllowableOffset=&geometryPrecision=&outSR=EPSG%3A4326&having=&returnIdsOnly=false&returnCountOnly=false&orderByFields=&groupByFieldsForStatistics=&outStatistics=&returnZ=false&returnM=false&gdbVersion=&historicMoment=&returnDistinctValues=false&resultOffset=&resultRecordCount=&queryByDistance=&returnExtentOnly=false&datumTransformation=&parameterValues=&rangeValues=&quantizationParameters=&featureEncoding=esriDefault&f=geojson'
    # print(url)
    r = requests.get(url)
    resp = r.json()
    # print(resp)
    if resp['features'] != []:  # [0]['type']
        # print(f'{res_type} DEM exists at this point')
        return True
    if resp['features'] == []:
        # print('False')
        # print(f"{res_type} not available")
        return  False

#query_dem([44.99857142857143, -93.10142857142857, 45.10142857142857, -92.99857142857142], '10m')

def query_dems(shape_type, coords, width=100):
    resp = {}
    bbox = make_bbox(shape_type, coords, width)
    # print(bbox)
    for res_type in res_types:
        resp[res_type] = (query_dem(bbox, res_type))

    print(resp)
    return(resp)

query_dems('point', [(39,-96)])
# query_dems('line', [(39,-96.0001), (39,-96)])
# query_dems('polygon', [(39,-96.0008), (39,-96), (39.0008, -96)])





# def get_flowtrace():

#     payload = {
#         "inputs": [
#             {
#             "id": "lat",
#             "type": "text/plain",
#             "value": "45"
#             },
#             {
#             "id": "lng",
#             "type": "text/plain",
#             "value": "-93"
#             },
#             {
#             "id": "raindroptrace",
#             "type": "text/plain",
#             "value": "True"
#             },
#             {
#             "id": "direction",
#             "type": "text/plain",
#             "value": "down"
#             }
#         ]
#         }
#     flowtrace_url = 'https://nhgf.wim.usgs.gov/processes/nldi-flowtrace/jobs?response=document'
    
#     r = requests.post(flowtrace_url, data=payload)
#     resp = r.json()
#     print(r)

# get_flowtrace()



# https://index.nationalmap.gov/arcgis/rest/services/3DEPElevationIndex/MapServer/21/query?where=&text=&objectIds=&time=&geometry=%7B%22xmin%22%3A-93.10142857142857%2C%22ymin%22%3A44.99857142857143%2C%22xmax%22%3A-92.99857142857142%2C%22ymax%22%3A45.10142857142857%2C%22spatialReference%22%3A%7B%22wkid%22%3A4326%7D%7D&geometryType=esriGeometryEnvelope&inSR=EPSG%3A4326&spatialRel=esriSpatialRelWithin&relationParam=&outFields=&returnGeometry=false&returnTrueCurves=false&maxAllowableOffset=&geometryPrecision=&outSR=EPSG%3A4326&having=&returnIdsOnly=false&returnCountOnly=false&orderByFields=&groupByFieldsForStatistics=&outStatistics=&returnZ=false&returnM=false&gdbVersion=&historicMoment=&returnDistinctValues=false&resultOffset=&resultRecordCount=&queryByDistance=&returnExtentOnly=false&datumTransformation=&parameterValues=&rangeValues=&quantizationParameters=&featureEncoding=esriDefault&f=geojson
# https://index.nationalmap.gov/arcgis/rest/services/3DEPElevationIndex/MapServer/21/query?where=&text=&objectIds=&time=&geometry=%7B%22xmin%22%3A-93.00142857142858%2C%22ymin%22%3A44.99857142857143%2C%22xmax%22%3A-92.99857142857142%2C%22ymax%22%3A45.00142857142857%2C%22spatialReference%22%3A%7B%22wkid%22%3A4326%7D%7D&geometryType=esriGeometryEnvelope&inSR=EPSG%3A4326&spatialRel=esriSpatialRelWithin&relationParam=&outFields=&returnGeometry=false&returnTrueCurves=false&maxAllowableOffset=&geometryPrecision=&outSR=EPSG%3A4326&having=&returnIdsOnly=false&returnCountOnly=false&orderByFields=&groupByFieldsForStatistics=&outStatistics=&returnZ=false&returnM=false&gdbVersion=&historicMoment=&returnDistinctValues=false&resultOffset=&resultRecordCount=&queryByDistance=&returnExtentOnly=false&datumTransformation=&parameterValues=&rangeValues=&quantizationParameters=&featureEncoding=esriDefault&f=geojson