from utils import geom_to_geojson, get_local_catchment, get_local_flowlines, get_coordsys, project_point
from utils import get_flowgrid, split_catchment, get_onFlowline, get_upstream_basin, merge_geometry

class Delineate:
    """Define inputs and outputs for the main Delineate class"""

    def __init__(self, x=None, y=None):

        self.x = x
        self.y = y
        self.catchmentIdentifier = None
        self.flowlines = None
        self.flw = None
        self.flwdir_transform = None
        self.projected_xy = None
        self.onFlowline = bool
        self.nhdCellList = None


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
        self.onFlowline = get_onFlowline(self.projected_xy, self.flowlines, self.transformToRaster, self.transformToWGS84)
        self.catchment = geom_to_geojson(self.catchmentGeom, 'catchment')

        if self.onFlowline == False:
            self.splitCatchment = geom_to_geojson(self.splitCatchmentGeom, 'splitCatchment')
            
        if self.onFlowline == True:
            self.upstreamBasinGeom = get_upstream_basin(self.catchmentIdentifier)
            self.mergedCatchmentGeom = merge_geometry(self.catchmentGeom, self.splitCatchmentGeom, self.upstreamBasinGeom)
            # self.upstreamBasin = self.geom_to_geojson(self.upstreamBasinGeom, 'upstreamBasin')
            self.mergedCatchment = geom_to_geojson(self.mergedCatchmentGeom, 'mergedCatchment')
            

    