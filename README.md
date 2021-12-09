# Mesh Parameterization
## Setup
To run this project, download the folder mesh_parametrization and import the package by using 
```
from mesh_parametrization import mesh_deformation
```
## Usage
### Import mesh
Import the mesh by using the function `read_mesh`:
```
mesh = read_mesh("LProfilTetMsh.xdmf")
```

### Import nodes and displacements
The nodes and displacements can be imported by calling the function `read_input_file` 
```
nodes, displacements = read_input_file("example_input.txt")
```
### Mesh deformation
Lastly, the mesh can be deformed by calling the function `deform_mesh`
```
deform_mesh(mesh, nodes, displacements)
```

# Corner Detection

Slightly modified and adapted implementation of:  
_3D INTEREST POINT DETECTION USING LOCAL SURFACE CHARACTERISTICS  
WITH APPLICATION IN ACTION RECOGNITION,  
Michael B. Holte, DOI: 10.1109/ICIP.2014.7026160_

## Usage
Import the corner detection function using
```
from cornerdetection_DoN.py import cornerdetection
```
Detect corners by calling the function on a given mesh
```
cornerdetection(mesh, "testoutput.txt", mode="angleSum", max_search_depth=2)
```
Where `mesh` is the mesh on which corners should be detected and `"testoutput.txt"` is the file name where the output should be written. Additionally, the two parameters `mode` and `max_search_depth` can be specified if needed. The parameter `mode` can either be `"DoN"` or `"angleSum"`. These are two different methods of finding the maximum corner likelyness of a corner compared to its surrounding corners. The parameter `max_search_depth` specifies in a neighbourhood of how many nodes a maximum should be searched.