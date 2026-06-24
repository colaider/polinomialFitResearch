import numpy as np

class ShapeGenrator:
    def __init__(self, shape, range_x = (-20,20), range_y = (-20, 20), resolution=0.04):
        self.shape = shape
        self.range_x = range_x
        self.range_y = range_y
        self.resolution = resolution
        self.noise_std = 0.0        
        self.point = np.dtype([ ('x', np.float32), ('y', np.float32), ('z', np.float32) ])
        self.points = []

    def addGrooves(self):
        pass
    
    def interpolateShape(self):
        xs = np.arange(self.range_x[0], self.range_x[1], self.resolution)
        ys = np.arange(self.range_y[0], self.range_y[1], self.resolution)    
        XX, YY = np.meshgrid(xs, ys)
        ZZ = self.shape(XX, YY).astype(np.float32)
        self.points = np.zeros(XX.size, dtype=self.point)
        self.points['x'] = XX.ravel()
        self.points['y'] = YY.ravel()
        self.points['z'] = ZZ.ravel()
        
        
    def addGroove(self, center_x, width, depth, wall_width=1.0, center_y=None, length=None):
        x = self.points['x']
        y = self.points['y']
        half = width / 2.0

        # y mask
        if center_y is not None and length is not None:
            y_mask = np.abs(y - center_y) < length / 2.0
        else:
            y_mask = np.ones(len(self.points), dtype=bool)

        dist = np.abs(x - center_x)

        # flat bottom
        flat_mask = (dist < half) & y_mask

        # smooth wall using cosine blend
        wall_mask = (dist >= half) & (dist < half + wall_width) & y_mask

        # t goes 0 at groove edge → 1 at surface
        t = (dist[wall_mask] - half) / wall_width

        # cosine gives natural smooth rounded fillet
        smooth = depth * 0.5 * (1.0 + np.cos(np.pi * t))

        self.points['z'][flat_mask] -= depth
        self.points['z'][wall_mask] -= smooth
                
    def addNoise(self, std):
        """
        std : noise standard deviation in mm (e.g. 0.003 for 3 micron scanner noise)
        """
        noise = np.random.normal(0, std, len(self.points))
        self.points['z'] += noise.astype(np.float32)