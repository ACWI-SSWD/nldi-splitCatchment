from utils import geom_to_geojson, get_local_catchment, get_local_flowlines, get_coordsys, project_point
from utils import get_flowgrid, get_onFlowline, get_raindropPath, get_intersectionPoint, get_reachMeasure, split_flowline, merge_downstreamPath
import geojson

class Flowtrace:
    """Define inputs and outputs for the main Flowtrace class"""

    def __init__(self, x=None, y=None, raindropTrace=bool, direction=str):

        self.x = x
        self.y = y
        self.raindropTrace = raindropTrace
        self.direction = direction
        self.catchmentIdentifier = None
        self.flowlines = None
        self.flw = None
        self.flwdir_transform = None
        self.projected_xy = None
        self.onFlowline = bool

        #geoms
        self.catchmentGeom = None
        self.splitCatchmentGeom = None
        self.upstreamBasinGeom = None
        self.mergedCatchmentGeom = None 
        self.intersectionPointGeom = None
        self.raindropPathGeom = None   
        self.nhdFlowlineGeom = None
        self.upstreamFlowlineGeom = None
        self.downstreamFlowlineGeom = None
        self.downstreamPathGeom = None

        #outputs
        self.catchment = None
        self.splitCatchment = None
        self.upstreamBasin = None
        self.mergedCatchment = None
        self.intersectionPoint = None
        self.raindropPath = None
        self.nhdFlowline = None
        self.streamInfo = None
        self.upstreamFlowline = None
        self.downstreamFlowline = None
        self.downstreamPath = None

        #create transform
        self.transformToRaster = None  
        self.transformToWGS84 = None 

        #kick off
        self.run()

    def serialize(self):
        # return {
        #     'catchment': self.catchment,
        #     'intersectionPoint': self.intersectionPoint,
        #     'raindropPath': self.raindropPath,
        #     'nhdFlowline': self.nhdFlowline,
        #     'streamInfo': self.streamInfo, 
        #     'upstreamFlowline': self.upstreamFlowline, 
        #     'downstreamFlowline': self.downstreamFlowline,
        #     'downstreamPath': self.downstreamPath
        # }
        if self.onFlowline == True:
            if self.direction == 'up':
                feature1 = geojson.Feature(geometry=self.upstreamFlowline, id='upstreamFlowline', property=self.streamInfo)
                featurecollection = geojson.FeatureCollection([feature1])
                return featurecollection
                # return {
                #     'streamInfo': self.streamInfo,
                #     'upstreamFlowline': self.upstreamFlowline
                # }
            if self.direction == 'down':
                feature1 = geojson.Feature(geometry=self.downstreamFlowline, id='downstreamFlowline', property=self.streamInfo)
                featurecollection = geojson.FeatureCollection([feature1])
                return featurecollection
                # return {                    
                #     'streamInfo': self.streamInfo, 
                #     'downstreamFlowline': self.downstreamFlowline
                # }
            if self.direction == 'none':
                feature1 = geojson.Feature(geometry=self.nhdFlowline, id='nhdFlowline', property=self.streamInfo)
                featurecollection = geojson.FeatureCollection([feature1])
                return featurecollection
                # return {
                #     'nhdFlowline': self.nhdFlowline,
                #     'streamInfo': self.streamInfo
                # }
        
        if self.onFlowline == False:
            if self.direction == 'up' and self.raindropTrace == True:
                feature1 = geojson.Feature(geometry=self.upstreamFlowline, id='upstreamFlowline', property=self.streamInfo)
                feature2 = geojson.Feature(geometry=self.raindropPath, id='raindropPath')
                featurecollection = geojson.FeatureCollection([feature1, feature2])
                return featurecollection
                # return {
                #     'raindropPath': self.raindropPath,
                #     'streamInfo': self.streamInfo, 
                #     'upstreamFlowline': self.upstreamFlowline
                # }
            if self.direction == 'down' and self.raindropTrace == True:
                feature1 = geojson.Feature(geometry=self.downstreamPath, id='downstreamPath', property=self.streamInfo)
                featurecollection = geojson.FeatureCollection([feature1])
                return featurecollection
                # return {
                #     'streamInfo': self.streamInfo,  
                #     'downstreamPath': self.downstreamPath
                # }
            if self.direction == 'none' and self.raindropTrace == True:
                feature1 = geojson.Feature(geometry=self.nhdFlowline, id='nhdFlowline', property=self.streamInfo)
                feature2 = geojson.Feature(geometry=self.raindropPath, id='raindropPath')
                featurecollection = geojson.FeatureCollection([feature1, feature2])
                return featurecollection
                # return {
                #     'raindropPath': self.raindropPath,
                #     'nhdFlowline': self.nhdFlowline,
                #     'streamInfo': self.streamInfo
                # }
            if self.direction == 'up' and self.raindropTrace == False:
                feature1 = geojson.Feature(geometry=self.upstreamPath, id='upstreamPath', property=self.streamInfo)
                featurecollection = geojson.FeatureCollection([feature1])
                return featurecollection
                # return {
                #     'streamInfo': self.streamInfo, 
                #     'upstreamFlowline': self.upstreamFlowline
                # }
            if self.direction == 'down' and self.raindropTrace == False:
                feature1 = geojson.Feature(geometry=self.downstreamFlowline, id='downstreamFlowline', property=self.streamInfo)
                featurecollection = geojson.FeatureCollection([feature1])
                return featurecollection
                # return {
                #     'streamInfo': self.streamInfo,  
                #     'downstreamFlowline': self.downstreamFlowline
                # }
            if self.direction == 'none' and self.raindropTrace == False:
                feature1 = geojson.Feature(geometry=self.nhdFlowline, id='nhdFlowline', property=self.streamInfo)
                featurecollection = geojson.FeatureCollection([feature1])
                return featurecollection
                # return {
                #     'nhdFlowline': self.nhdFlowline,
                #     'streamInfo': self.streamInfo
                # }
                

## main functions
    def run(self):
        # Order of these functions is important!
        self.catchmentIdentifier, self.catchmentGeom = get_local_catchment(self.x, self.y)
        self.flowlines, self.nhdFlowlineGeom = get_local_flowlines(self.catchmentIdentifier)
        self.transformToRaster, self.transformToWGS84 = get_coordsys()
        self.projected_xy = project_point(self.x, self.y, self.transformToRaster)
        self.flw, self.flwdir_transform = get_flowgrid( self.catchmentGeom, self.transformToRaster, self.transformToWGS84 )
        self.onFlowline = get_onFlowline( self.projected_xy, self.flowlines, self.transformToRaster, self.transformToWGS84)
        self.catchment = geom_to_geojson(self.catchmentGeom, 'catchment')


        if self.onFlowline == True:
            self.intersectionPointGeom = get_intersectionPoint(self.x, self.y, self.onFlowline)
            self.streamInfo = get_reachMeasure(self.intersectionPointGeom, self.flowlines)
            self.upstreamFlowlineGeom, self.downstreamFlowlineGeom = split_flowline(self.intersectionPointGeom, self.flowlines)

            # Outputs
            if self.direction == 'up':
                self.upstreamFlowline = geom_to_geojson(self.upstreamFlowlineGeom, 'upstreamFlowline')
            
            if self.direction == 'down':
                self.downstreamFlowline = geom_to_geojson(self.downstreamFlowlineGeom, 'downstreamFlowline')
            
            if self.direction == 'none':
                self.nhdFlowline = geom_to_geojson(self.nhdFlowlineGeom, 'nhdFlowline')

        if self.onFlowline == False:
            self.raindropPathGeom = get_raindropPath(self.flw, self.projected_xy,  self.nhdFlowlineGeom, self.flowlines, self.transformToRaster, self.transformToWGS84)
            self.intersectionPointGeom = get_intersectionPoint(self.x, self.y, self.onFlowline, self.raindropPathGeom)
            self.streamInfo = get_reachMeasure(self.intersectionPointGeom, self.flowlines)
            self.upstreamFlowlineGeom, self.downstreamFlowlineGeom = split_flowline(self.intersectionPointGeom, self.flowlines)            

            # Outputs
            if self.direction == 'up' and self.raindropTrace == True:
                self.upstreamFlowline = geom_to_geojson(self.upstreamFlowlineGeom, 'upstreamFlowline')
                self.raindropPath = geom_to_geojson(self.raindropPathGeom, 'raindropPath')
            
            if self.direction == 'down' and self.raindropTrace == True:
                # self.downstreamFlowline = geom_to_geojson(self.downstreamFlowlineGeom, 'downstreamFlowline')
                # self.raindropPath = geom_to_geojson(self.raindropPathGeom, 'raindropPath')
                self.downstreamPathGeom = merge_downstreamPath(self.raindropPathGeom, self.downstreamFlowlineGeom)
                self.downstreamPath = geom_to_geojson(self.downstreamPathGeom, 'downstreamPath')
            
            if self.direction == 'none' and self.raindropTrace == True:
                self.nhdFlowline = geom_to_geojson(self.nhdFlowlineGeom, 'nhdFlowline')
                self.raindropPath = geom_to_geojson(self.raindropPathGeom, 'raindropPath')
            
            if self.direction == 'up' and self.raindropTrace == False:
                 self.upstreamFlowline = geom_to_geojson(self.upstreamFlowlineGeom, 'upstreamFlowline')
            
            if self.direction == 'down' and self.raindropTrace == False:
                self.downstreamFlowline = geom_to_geojson(self.downstreamFlowlineGeom, 'downstreamFlowline')
            
            if self.direction == 'none' and self.raindropTrace == False:
                self.nhdFlowline = geom_to_geojson(self.nhdFlowlineGeom, 'nhdFlowline')

        # self.intersectionPoint = geom_to_geojson(self.intersectionPointGeom, 'intersectionPoint')