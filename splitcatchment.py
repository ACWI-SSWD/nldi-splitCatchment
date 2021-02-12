from utils import geom_to_geojson, get_local_catchment, get_local_flowlines, get_coordsys, project_point
from utils import get_total_basin, split_catchment, get_onFlowline, get_upstream_basin, merge_geometry

class SplitCatchment:
    """Define inputs and outputs for the main SplitCatchment class"""

    def __init__(self, x=None, y=None, upstream=bool):

        self.x = x
        self.y = y
        self.catchmentIdentifier = None
        self.flowlines = None
        self.flw = None
        self.flwdir_transform = None
        self.projected_xy = None
        self.onFlowline = bool
        self.upstream = upstream

        #geoms
        self.catchmentGeom = None
        self.splitCatchmentGeom = None
        self.totalBasinGeom = None
        self.upstreamBasinGeom = None
        self.mergedCatchmentGeom = None 
        self.nhdFlowlineGeom = None

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
        if self.upstream == False:
            return {
                'catchment': self.catchment,
                'splitCatchment': self.splitCatchment, 
            }
        if self.upstream == True and self.onFlowline == True:
            return {
                'catchment': self.catchment,
                'mergedCatchment': self.mergedCatchment
            }
        if self.upstream == True and self.onFlowline == False:
            return {
                'catchment': self.catchment,
                'splitCatchment': self.splitCatchment,
                'upstreamBasin': self.upstreamBasin
            }
        # return {
        #     'catchment': self.catchment,
        #     'splitCatchment': self.splitCatchment,
        #     'upstreamBasin': self.upstreamBasin,
        #     'mergedCatchment': self.mergedCatchment
        # }

## main functions
    def run(self):
        # Order of these functions is important!
        self.catchmentIdentifier, self.catchmentGeom = get_local_catchment(self.x, self.y)
        self.flowlines, self.nhdFlowlineGeom = get_local_flowlines(self.catchmentIdentifier)
        self.transformToRaster, self.transformToWGS84 = get_coordsys()
        self.projected_xy = project_point(self.x, self.y, self.transformToRaster)
        self.splitCatchmentGeom = split_catchment(self.catchmentGeom, self.projected_xy, self.transformToRaster, self.transformToWGS84)
        self.onFlowline = get_onFlowline(self.projected_xy, self.flowlines, self.transformToRaster, self.transformToWGS84)
        self.catchment = geom_to_geojson(self.catchmentGeom, 'catchment')

        # outputs
        if self.upstream == False:
            self.splitCatchment = geom_to_geojson(self.splitCatchmentGeom, 'splitCatchment')

        if self.upstream == True and self.onFlowline == True:
            self.totalBasinGeom = get_total_basin(self.catchmentIdentifier, self.catchmentGeom)     
            self.mergedCatchmentGeom = merge_geometry(self.catchmentGeom, self.splitCatchmentGeom, self.totalBasinGeom)
            self.mergedCatchment = geom_to_geojson(self.mergedCatchmentGeom, 'mergedCatchment')
        
        if self.upstream == True and self.onFlowline == False: 
            self.splitCatchment = geom_to_geojson(self.splitCatchmentGeom, 'splitCatchment')
            self.totalBasinGeom = get_total_basin(self.catchmentIdentifier, self.catchmentGeom)
            self.upstreamBasinGeom = get_upstream_basin(self.catchmentGeom, self.totalBasinGeom)
            self.upstreamBasin = geom_to_geojson(self.upstreamBasinGeom, 'upstreamBasin')
    