## StreamStats delineation script flask app server wrapper

# -----------------------------------------------------
# Martyn Smith USGS
# 09/30/2019
# StreamStats Delineation script flask server
#
# run with: "python -m flask run"
#
# -----------------------------------------------------

from flask import Flask, request, jsonify
from flask_cors import CORS, cross_origin
from datetime import datetime
import splitcatchment
import flowtrace
import nldi_xstool 
import dem_query
import time
import json
from shapely.geometry import LineString, mapping
import geojson

app = Flask(__name__)
cors = CORS(app)
app.config['CORS_HEADERS'] = 'Content-Type'

@app.route("/")

def home():
    return "sample delineation server"

@app.route("/delineate")
@cross_origin(origin='*')
def main():

    lat = float(request.args.get('lat'))
    lng = float(request.args.get('lng'))
    runsplitcatchment = request.args.get('runsplitcatchment')
    truefalse = bool(request.args.get('truefalse'))
    direction = request.args.get('direction')
    query = request.args.get('query')
    res = request.args.get('res')
    xstool = request.args.get('xstool')

    #start main program
    timeBefore = time.perf_counter()  
    
############# Splitcatchment ##############
    if runsplitcatchment == 'true':
        results = splitcatchment.SplitCatchment(lng, lat, truefalse)
        results = jsonify(results.serialize())

############### Flowtrace ###############
    if runsplitcatchment == 'false' and xstool == 'false' and query == 'false':
        results = flowtrace.Flowtrace(lng, lat, truefalse, direction)
        results = jsonify(results.serialize())

############# Dem_query ##################
    if runsplitcatchment == 'false' and xstool == 'false' and query == 'true':
        output = dem_query.query_dems('point', [(lat,lng)])
        results = output
        results = json.dumps(results)
        results = json.loads(results)

############# XSTool ###################
    if runsplitcatchment == 'false' and xstool == 'true':
        output = nldi_xstool.getXSAtPoint((lat, lng), 100, 1000, res, None)
        output = output.to_json()
        output = json.loads(output)

        # print('output',type(output), output)

        line = []
        x = 0
        while x < 100:
            p = (output['features'][x]['geometry']['coordinates'][0], output['features'][x]['geometry']['coordinates'][1])
            line.append(p)
            x += 1

        linestring = LineString(line)
        geojson_dict = mapping(linestring)
        feature1 = geojson.Feature(geometry=geojson_dict, id='streamCrossSection')
        results = geojson.FeatureCollection([feature1])


    
    timeAfter = time.perf_counter() 
    totalTime = timeAfter - timeBefore
    print("Total Time:",totalTime)
    print("results: ", type(results) , results)
    return results