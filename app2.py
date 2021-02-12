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
import time

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
    # print('runsplitcatchment: ', type(runsplitcatchment), runsplitcatchment, 'truefalse: ', type(truefalse) , truefalse, 'direction: ', direction)

    print(lat,lng)

    #start main program
    timeBefore = time.perf_counter()  
    
    if runsplitcatchment == 'true':
        results = splitcatchment.SplitCatchment(lng, lat, truefalse)

    if runsplitcatchment == 'false':
        results = flowtrace.Flowtrace(lng, lat, truefalse, direction)

    # print("results: ", type(results) , results, results.serialize())
    
    timeAfter = time.perf_counter() 
    totalTime = timeAfter - timeBefore
    print("Total Time:",totalTime)
    return jsonify(results.serialize())