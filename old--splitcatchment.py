from utils import geom_to_geojson, get_local_catchment, get_coordsys, project_point, get_flowgrid, split_catchment


class SplitCatchment:
    """Define inputs and outputs for the main SplitCatchment class"""

    def __init__(self, x=None, y=None):

        self.x = x
        self.y = y
        self.catchmentIdentifier = None
        self.flw = None
        self.flwdir_transform = None
        self.projected_xy = None

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
        self.transformToRaster, self.transformToWGS84 = get_coordsys()
        self.projected_xy = project_point(self.x, self.y, self.transformToRaster)
        self.flw, self.flwdir_transform = get_flowgrid(self.catchmentGeom, self.transformToRaster, self.transformToWGS84)
        self.splitCatchmentGeom  = split_catchment(self.flw, self.flwdir_transform, self.projected_xy, self.transformToWGS84)

        # Outputs
        self.catchment = geom_to_geojson(self.catchmentGeom, 'catchment')
        self.splitCatchment = geom_to_geojson(self.splitCatchmentGeom, 'splitCatchment')
    
