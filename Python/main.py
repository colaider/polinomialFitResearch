from plyGenrator import PlyGenartor

shape = lambda x, y: 0.000005*x**4 - 0.002*x**2 + 0.00005*y**4 - 0.002*y**2
s = PlyGenartor(shape)
s.interpolateShape()
s.addGroove(center_x=-4.5, width=2.0, depth=2.0, wall_width=1.0, center_y=0, length=10.0)
s.addGroove(center_x=4.5, width=2.0, depth=2.0, wall_width=1.0, center_y=0, length=10.0)
s.addNoise(std=0.04)   # 3 micron noise
s.save_PLY(r'C:\Users\ominaiev1\source\repos\PolinomilPlaneWithCurvatireDapation\x64\Debug\ply\output.ply')
