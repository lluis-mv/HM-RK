
"""
"""

import gmsh
from math import cos, sin, pi


gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)

gmsh.model.add("modelo_2")


# %% Variables asociadas a la geometría:
L = 6.5/2
r = 1

# %% Se crea la geometría:

# La siguiente línea permite abreviar el código:
geom = gmsh.model.geo

# Puntos:
pc = geom.addPoint(0, 0, 0)  # Punto central del círculo
p2 = geom.addPoint(L, 0, 0)
p3 = geom.addPoint(L, L, 0)
p4 = geom.addPoint(0, L, 0)



# Líneas
l1 = geom.addLine(pc, p2)
l2 = geom.addLine(p2, p3)
l3 = geom.addLine(p3, p4)
l4 = geom.addLine(p4, pc)


# Curve loops:
cl1 = geom.addCurveLoop([l1, l2, l3, l4])




s1 = geom.addPlaneSurface([cl1])

# Definimos la superficie física:
gmsh.model.addPhysicalGroup(2, [s1], 101)
gmsh.model.setPhysicalName(2, 101, 'mi superficie')


# %% Definimos que la superficie creada tenga una malla estructurada:
geom.mesh.setTransfiniteSurface(s1)   



# Y recombinamos para que use elementos cuadriláteros en lugar de triángulos
geom.mesh.setRecombine(2, s1)



# %% Finalmente se crea la malla y se guarda:

geom.synchronize()

# Ver las "caras" de los elementos finitos 2D
gmsh.option.setNumber('Mesh.SurfaceFaces', 1)

# Ver los nodos de la malla
gmsh.option.setNumber('Mesh.Points', 1)

gmsh.model.mesh.generate(2)

# Se guarda la malla
filename = 'ejm_2.msh'
gmsh.write(filename)

# Podemos visualizar el resultado en la interfaz gráfica de GMSH
gmsh.fltk.run()

# Se finaliza el programa
gmsh.finalize()

# %% Podemos graficar la malla para ver el resultado:

from leer_GMSH import plot_msh  # Funciones para leer y graficar la malla

#plot_msh(filename, '2D', True)