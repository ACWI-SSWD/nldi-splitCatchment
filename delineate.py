from utils import geom_to_geojson, transform_geom, get_local_catchment, get_local_flowlines
from utils import split_catchment, get_onFlowline, get_upstream_basin, merge_geometry

class Delineate:
    """Define inputs and outputs for the main Delineate class"""

    def __init__(self, x=None, y=None):

        self.x = x
        self.y = y
        self.catchmentIdentifier = None
        self.flowlines = None
        self.flw = None
        self.xy = None
        self.clickonFlowline = bool
        self.nhdCellList = None


        #geoms
        self.catchmentGeom = None
        self.splitCatchmentGeom = None
        self.upstreamBasinGeom = None
        self.mergedCatchmentGeom = None 
        self.intersectionPointGeom = None
        self.downstreamPathGeom = None   
        self.nhdFlowlineGeom = None
        self.upstreamFlowlineGeom = None
        self.downstreamFlowlineGeom = None

        #outputs
        self.catchment = None
        self.splitCatchment = None
        self.upstreamBasin = None
        self.mergedCatchment = None
        self.intersectionPoint = None
        self.downstreamPath = None
        self.nhdFlowline = None
        self.streamInfo = None
        self.upstreamFlowline = None
        self.downstreamFlowline = None

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
            'nhdFlowline': self.nhdFlowline,
            'streamInfo': self.streamInfo, 
            'upstreamFlowline': self.upstreamFlowline, 
            'downstreamFlowline': self.downstreamFlowline
        }

## main functions
    def run(self):
        # Order of these functions is important!
        self.catchmentIdentifier, self.catchmentGeom = get_local_catchment(self.x, self.y)
        self.flowlines, self.nhdFlowlineGeom = get_local_flowlines(self.catchmentIdentifier)
        self.splitCatchmentGeom, self.flw, self.xy, self.transformToRaster, self.transformToWGS84 = split_catchment(self.catchmentGeom, self.x, self.y)
        self.clickonFlowline, self.intersectionPointGeom  = get_onFlowline( self.flw, self.xy, self.flowlines, self.transformToRaster, self.transformToWGS84)
        self.catchment = geom_to_geojson(self.catchmentGeom, 'catchment')

        if self.clickonFlowline == False:
            self.splitCatchment = geom_to_geojson(self.splitCatchmentGeom, 'splitCatchment')
            
        if self.clickonFlowline == True:
            self.upstreamBasinGeom = get_upstream_basin(self.catchmentIdentifier)
            self.mergedCatchmentGeom = merge_geometry(self.catchmentGeom, self.splitCatchmentGeom, self.upstreamBasinGeom)
            # self.upstreamBasin = self.geom_to_geojson(self.upstreamBasinGeom, 'upstreamBasin')
            self.mergedCatchment = geom_to_geojson(self.mergedCatchmentGeom, 'mergedCatchment')
            

    