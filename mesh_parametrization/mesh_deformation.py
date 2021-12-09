from dolfin import *
from vtkplotter.dolfin import plot, screenshot
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d


# Optimization options for the form compiler
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "uflacs"

# Create mesh and define function space


def deform_mesh(mesh, nodes, displacements):
    V = VectorFunctionSpace(mesh, "Lagrange", 1)

    #bcl = DirichletBC(V, c, left)
    bcs = []
    for i in range(len(nodes)):
        value = Constant(displacements[i])
        bound = Boundary(nodes[i])
        print("creating BC at " + str(bound) + " with u = " + str(displacements[i]))
        bc = DirichletBC(V, value, bound, method='pointwise')
        bcs.append(bc)

    # Define variational problem
    du = TrialFunction(V)
    d = du.geometric_dimension()  # space dimension
    v = TestFunction(V)

    # Compute solution
    u = Function(V)
    F = inner(grad(v), sigma(u)) * dx

    # Compute Jacobian of F
    J = derivative(F, u, du)

    # Solve variational problem
    solve(F == 0, u, bcs, J=J)

    file_displacement = File("displacement.pvd")
    file_displacement << u

    # file_displacement = File("mesh.pvd")
    # file_displacement << mesh
    u_mesh = project(u, V)
    # plot(mesh)
    ALE.move(mesh, u_mesh)
    # plot(mesh,linewidth=0.3,color="black")
    #plot(u)
    #plot(mesh)
    plot(u, mode="displacement")
    # lgd = legend(plots_source=None)
    # lgd.set_entry(label="Blue Square", color="blue")
    #screenshot()
    # screenshot()

    V = FunctionSpace(mesh, 'P', 1)

    # Compute magutude of displacement
    u_magnitude = sqrt(dot(u, u))
    u_magnitude = project(u_magnitude, V)
    print('min/max u:', u_magnitude.vector().min(), u_magnitude.vector().max() / 1000)

def epsilon(u):
    return 0.5 * (nabla_grad(u) + nabla_grad(u).T)

    # return sym(nabla_grad(u))

def sigma(u):
    #todo optionally, change dimension based on mesh
    d = 3
    return epsilon(u) - (tr(epsilon(u)) * Identity(d) / 3)


# Define boundary condition class
class Boundary(SubDomain):
    # input argument x is the coordinates of the displacement in x[0], x[1], x[2]
    def __init__(self, x):
        super().__init__()
        self.x = x

    def __str__(self):
        return 'Boundary at (x = ' + str(self.x[0]) + ' , y = ' + str(self.x[1]) + ' , z = ' + str(self.x[1]) + ')'

    def inside(self, x, on_boundary):
        tol = 1e-12
        # todo: change from choosing based on coordinates to based on node IDX or similar (mark node maybe?)
        return near(x[0], self.x[0], tol) and near(x[1], self.x[1], tol) and near(x[2], self.x[2], tol)


def read_mesh(filename):
    #todo: add input check, if file is in the correct format etc.
    mesh = Mesh()
    with XDMFFile(filename) as f:
        f.read(mesh)
    return mesh

def read_position(line, position):
    return line.split(",")[position].strip()

def read_input_file(filename):
    with open(filename, "r") as f:
        input_from_CD = f.readlines()
        input_from_CD = input_from_CD[1:]
    nodes = [(float(read_position(line, 2)), float(read_position(line, 3)), float(read_position(line, 4))) for line in input_from_CD]
    displacements = [(float(read_position(line, 5)), float(read_position(line, 6)), float(read_position(line, 7))) for line in input_from_CD]
    print("nodes imported: " + str(nodes))
    print("displacements imported: " + str(displacements))
    return nodes, displacements
