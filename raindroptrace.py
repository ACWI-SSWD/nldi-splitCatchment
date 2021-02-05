from utils import geom_to_geojson, get_local_catchment, get_local_flowlines, get_coordsys, project_point
from utils import get_flowgrid, split_catchment, get_onFlowline, get_raindropPath, get_reachMeasure


class RaindropTrace:
    """Define inputs and outputs for the main RaindropTrace class"""

    def __init__(self, x=None, y=None, ):

        self.x = x
        self.y = y
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
            'mergedCatchment': self.mergedCatchment,
            'intersectionPoint': self.intersectionPoint,
            'raindropPath': self.raindropPath,
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
        self.transformToRaster, self.transformToWGS84 = get_coordsys()
        self.projected_xy = project_point(self.x, self.y, self.transformToRaster)
        self.flw, self.flwdir_transform = get_flowgrid( self.catchmentGeom, self.transformToRaster, self.transformToWGS84 )
        self.splitCatchmentGeom = split_catchment(self.flw, self.flwdir_transform, self.projected_xy, self.transformToWGS84)
        self.onFlowline  = get_onFlowline( self.flw, self.projected_xy, self.flowlines, self.transformToRaster, self.transformToWGS84)
        
        # self.catchment = geom_to_geojson(self.catchmentGeom, 'catchment')

        if self.onFlowline == False:
            self.raindropPathGeom, self.intersectionPointGeom = get_raindropPath(self.flw, self.projected_xy,  self.nhdFlowlineGeom, self.flowlines, self.transformToRaster, self.transformToWGS84)
            self.raindropPath = geom_to_geojson(self.raindropPathGeom, 'raindropPath')
            self.streamInfo, self.upstreamFlowlineGeom, self.downstreamFlowlineGeom = get_reachMeasure(self.intersectionPointGeom, self.flowlines, self.raindropPathGeom)
        
        if self.onFlowline == True:
            self.streamInfo, self.upstreamFlowlineGeom, self.downstreamFlowlineGeom = get_reachMeasure(self.intersectionPointGeom, self.flowlines)

        # Data returned
        # self.intersectionPoint = geom_to_geojson(self.intersectionPointGeom, 'intersectionPoint')
        # self.upstreamFlowline = geom_to_geojson(self.upstreamFlowlineGeom, 'upstreamFlowline')
        self.downstreamFlowline = geom_to_geojson(self.downstreamFlowlineGeom, 'downstreamFlowline')
        # self.nhdFlowline = geom_to_geojson(self.nhdFlowlineGeom, 'nhdFlowline')
