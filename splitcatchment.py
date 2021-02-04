from utils import geom_to_geojson, get_local_catchment, split_catchment


class SplitCatchment:
    """Define inputs and outputs for the main SplitCatchment class"""

    def __init__(self, x=None, y=None):

        self.x = x
        self.y = y
        self.catchmentIdentifier = None
        self.flw = None
        self.xy = None

        #geoms
        self.catchmentGeom = None
        self.splitCatchmentGeom = None
        
        #outputs
        self.catchment = None
        self.splitCatchment = None

        #create transform
        self.transformToRaster = None
        self.transformToWGS84 = None

        #kick off
        self.run()

    def serialize(self):
        return {
            'catchment': self.catchment,
            'splitCatchment': self.splitCatchment, 
        }

## main function
    def run(self):
        # Order of these functions is important!
        self.catchmentIdentifier, self.catchmentGeom = get_local_catchment(self.x, self.y)
        self.splitCatchmentGeom, self.flw, self.xy, self.transformToRaster, self.transformToWGS84 = split_catchment(self.catchmentGeom, self.x, self.y)
        self.catchment = geom_to_geojson(self.catchmentGeom, 'catchment')
        self.splitCatchment = geom_to_geojson(self.splitCatchmentGeom, 'splitCatchment')
    
