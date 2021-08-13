from ngsolve import *
from netgen.csg import *


L = 1
h = 0.1

n_1 = Vec(1, 1, 1)
n_2 = Vec(1, -1, 1)
n_3 = Vec(-1, 1, 1)
n_4 = Vec(-1, -1, 1)

n_5 = Vec(1, 1, -1)
n_6 = Vec(1, -1, -1)
n_7 = Vec(-1, 1, -1)
n_8 = Vec(-1, -1, -1)

p_1 = Pnt(L, L, 0)
p_2 = Pnt(L, -L, 0)
p_3 = Pnt(-L, L, 0)
p_4 = Pnt(-L, -L, 0)

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
    Pnt(-L*1.2, -L*1.2, -L*1.2),
    Pnt(L*1.2, L*1.2, L*1.2)
)

p_5 = Pnt(L-h, L-h, 0)
p_6 = Pnt(L-h, -L+h, 0)
p_7 = Pnt(-L+h, L-h, 0)
p_8 = Pnt(-L+h, -L+h, 0)

mini_up_left = Plane(p_5, n_1)
mini_up_front = Plane(p_6, n_2)
mini_up_right = Plane(p_7, n_3)
mini_up_back = Plane(p_8, n_4)

mini_bot_left = Plane(p_5, n_5)
mini_bot_front = Plane(p_6, n_6)
mini_bot_right = Plane(p_7, n_7)
mini_bot_back = Plane(p_8, n_8)

mini_octahedron = mini_up_left*mini_up_right*mini_up_front*mini_up_back\
    *mini_bot_left*mini_bot_right*mini_bot_front*mini_bot_back

mini_cube = OrthoBrick(
    Pnt((-L+h)*1.2, (-L+h)*1.2, (-L+h)*1.2),
    Pnt((L-h)*1.2, (L-h)*1.2, (L-h)*1.2)
)

hole = mini_octahedron*mini_cube

tkd = octahedron*cube

alv = tkd - hole

geo = CSGeometry()
geo.Add(alv)
mesh = geo.GenerateMesh(maxh=0.1)
mesh.Export("./raw-data/TKD.gmsh","Gmsh2 Format")