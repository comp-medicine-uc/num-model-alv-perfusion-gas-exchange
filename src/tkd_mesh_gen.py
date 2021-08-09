from ngsolve import *
from netgen.csg import *


n_1 = Vec(1, 1, 1)
n_2 = Vec(1, -1, 1)
n_3 = Vec(-1, 1, 1)
n_4 = Vec(-1, -1, 1)

n_5 = Vec(1, 1, -1)
n_6 = Vec(1, -1, -1)
n_7 = Vec(-1, 1, -1)
n_8 = Vec(-1, -1, -1)

p_1 = Pnt(1, 1, 0)
p_2 = Pnt(1, -1, 0)
p_3 = Pnt(-1, 1, 0)
p_4 = Pnt(-1, -1, 0)

up_left = Plane(p_1, n_1)
up_front = Plane(p_2, n_2)
up_right = Plane(p_3, n_3)
up_back = Plane(p_4, n_4)

bot_left = Plane(p_1, n_5)
bot_front = Plane(p_2, n_6)
bot_right = Plane(p_3, n_7)
bot_back = Plane(p_4, n_8)

octahedron = up_left*up_right*up_front*up_back\
    *bot_left*bot_right*bot_front*bot_back

cube = OrthoBrick(
    Pnt(-1, -1, -1),
    Pnt(1, 1, 1)
)

tkd = octahedron*cube

geo = CSGeometry()
geo.Add(tkd)
mesh = geo.GenerateMesh(maxh=0.1)
