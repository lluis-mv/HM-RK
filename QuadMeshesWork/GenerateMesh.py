# -*- coding: utf-8 -*-
import gmsh
from math import cos, sin, pi


gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)

gmsh.model.add("model")


# %% Variables asociadas a la geometría:
L = 1.0

# %% Se crea la geometría:

# La siguiente línea permite abreviar el código:
geom = gmsh.model.geo

# Puntos:
p1 = geom.addPoint(0, 0, 0)  # Punto central del círculo
p2 = geom.addPoint(0.2, 0, 0)
p3 = geom.addPoint(0.2, 1.0, 0)
p4 = geom.addPoint(0, 1.0, 0)


# Líneas
l1 = geom.addLine(p1, p2)
l2 = geom.addLine(p2, p3)
l3 = geom.addLine(p3, p4)
l4 = geom.addLine(p4, p1)

# Curve loops:
cl1 = geom.addCurveLoop([l1, l2, l3, l4])


# Dado que queremos una malla estructurada, definimos el número de elementos
# que queremos en ciertas líneas:
# Sintaxis: gmsh.model.geo.mesh.setTransfiniteCurve(tag, # nodos)

n = 15  # Líneas perpendiculares al círculo

for tag in [l1, l2, l3, l4]:
    gmsh.model.geo.mesh.setTransfiniteCurve(tag, n+1)

# Finalmente creamos la superficie:
s1 = geom.addPlaneSurface([cl1])
# Definimos la superficie física:
print(s1)
gmsh.model.addPhysicalGroup(2, [s1], 101)
gmsh.model.setPhysicalName(2, 101, 'my surface')


# %% Definimos que la superficie creada tenga una malla estructurada:
geom.mesh.setTransfiniteSurface(s1)   

# Y recombinamos para que use elementos cuadriláteros en lugar de triángulos
geom.mesh.setRecombine(2, s1)


geom.synchronize()

# Ver las "caras" de los elementos finitos 2D
gmsh.option.setNumber('Mesh.SurfaceFaces', 1)

# Ver los nodos de la malla
gmsh.option.setNumber('Mesh.Points', 1)

gmsh.model.mesh.generate(2)
gmsh.option.setNumber('Mesh.SecondOrderIncomplete', 1)
gmsh.model.mesh.setOrder(2)
# Se guarda la malla
filename = 'mesh.msh'
gmsh.write(filename)

#gmsh.fltk.run()

gmsh.finalize()

from leer_GMSH import plot_msh  # Funciones para leer y graficar la malla

plot_msh(filename, '2D', True)


