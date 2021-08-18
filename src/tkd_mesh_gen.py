'''TKD mesh generator using Netgen for alveolar simulations.
'''

__author__ = 'pzuritas'
__email__ = 'pzurita@uc.cl'

from netgen.csg import Vec, Plane, OrthoBrick, Pnt, CSGeometry
import meshio

def TKD_generator(
    R=100, h=6, maxh=2, save=True, name='./raw-data/TKD', convert=True
):
    '''Creates a TKD cuboid mesh as alveolar approximation using Netgen.

    R: alveolar radius in um. (int or float)
    h: blood-air barrier/alvelar wall thickness in um. (int or float)
    maxh: maximum element diameter h for mesher. (int or float)
    save: saves the mesh in Gmsh format. (bool)
    name: path for Gmsh file when saved. (string)
    convert: converts the Gmsh to VTK using meshio. (bool)
    '''

    # Plane directions for TKD

    n = {
        1: Vec(1, 1, 1),
        2: Vec(1, -1, 1),
        3: Vec(-1, 1, 1),
        4: Vec(-1, -1, 1),
        5: Vec(1, 1, -1),
        6: Vec(1, -1, -1),
        7: Vec(-1, 1, -1),
        8: Vec(-1, -1, -1)
    }

    # Points for geometry

    p = {
        1: Pnt(R, R, 0),
        2: Pnt(R, -R, 0),
        3: Pnt(-R, R, 0),
        4: Pnt(-R, -R, 0),
        5: Pnt(R-h, R-h, 0),
        6: Pnt(R-h, -R+h, 0),
        7: Pnt(-R+h, R-h, 0),
        8: Pnt(-R+h, -R+h, 0)
    }

    # Planes

    up_left = Plane(p[1], n[1])
    up_front = Plane(p[2], n[2])
    up_right = Plane(p[3], n[3])
    up_back = Plane(p[4], n[4])

    bot_left = Plane(p[1], n[5])
    bot_front = Plane(p[2], n[6])
    bot_right = Plane(p[3], n[7])
    bot_back = Plane(p[4], n[8])

    octahedron = up_left*up_right*up_front*up_back\
        *bot_left*bot_right*bot_front*bot_back

    cube = OrthoBrick(
        Pnt(-R*1.2, -R*1.2, -R*1.2),
        Pnt(R*1.2, R*1.2, R*1.2)
    )

    mini_up_left = Plane(p[5], n[1])
    mini_up_front = Plane(p[6], n[2])
    mini_up_right = Plane(p[7], n[3])
    mini_up_back = Plane(p[8], n[4])

    mini_bot_left = Plane(p[5], n[5])
    mini_bot_front = Plane(p[6], n[6])
    mini_bot_right = Plane(p[7], n[7])
    mini_bot_back = Plane(p[8], n[8])

    mini_octahedron = mini_up_left*mini_up_right*mini_up_front*mini_up_back\
        *mini_bot_left*mini_bot_right*mini_bot_front*mini_bot_back

    mini_cube = OrthoBrick(
        Pnt((-R+h)*1.2, (-R+h)*1.2, (-R+h)*1.2),
        Pnt((R-h)*1.2, (R-h)*1.2, (R-h)*1.2)
    )

    hole = mini_octahedron*mini_cube
    tkd = octahedron*cube
    alv = tkd - hole

    geo = CSGeometry()
    geo.Add(alv)
    mesh = geo.GenerateMesh(maxh=maxh)

    if save:
        mesh.Export(name+".gmsh","Gmsh2 Format")
        if convert:
            meshio_mesh = meshio.read(name+".gmsh", file_format="gmsh")
            meshio_mesh.write(name+".vtk")

    return mesh

if __name__ == '__main__':
    TKD_generator()