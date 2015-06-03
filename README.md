# PV-HighOrder

Original Author: Sébastien Blaise
http://perso.uclouvain.be/sebastien.blaise/index.html

Paraview plugin to visualise high order data

See http://perso.uclouvain.be/sebastien.blaise/tools.html for more info

Element node ordering see http://geuz.org/gmsh/doc/texinfo/gmsh.html#High-order-elements

The CellData arrays *\_HOsol\_i (where i is the solution point index) should be populated for all solution points in the cell according to the gmsh ordering above. (Note this is a change from the original code).

On startup if a PointData array that matches the *\_HOsol\_i is not found the high order solution is ignored.

