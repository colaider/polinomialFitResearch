from shapeGenrator import ShapeGenrator
import matplotlib.pyplot as plt 


class PlyGenartor(ShapeGenrator):
    def __init__(self, shape, range_x = (-10,10), range_y = (-20, 20), resolution=0.04):
        super().__init__(shape, range_x, range_y, resolution)
        

    def save_PLY(self, path):
        """Save to PLY file compatible with your C++ reader"""
        if self.points is None:
            raise RuntimeError("Call generate() first")

        with open(path, 'w') as f:
            f.write("ply\n")
            f.write("format ascii 1.0\n")
            f.write(f"element vertex {len(self.points)}\n")
            f.write("property float x\n")
            f.write("property float y\n")
            f.write("property float z\n")
            f.write("end_header\n")
            for p in self.points:
                f.write(f"{p['x']:.4f} {p['y']:.4f} {p['z']:.4f} \n")
        print(f"Saved {len(self.points)} points to {path}")
    
    
    
    def plot(self):
        fig = plt.figure()
        ax  = fig.add_subplot(111, projection='3d')
        ax.scatter(
            self.points['x'],
            self.points['y'],
            self.points['z'],
            s=0.1, c=self.points['z'], cmap='viridis'
        )
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        plt.tight_layout()
        plt.show()