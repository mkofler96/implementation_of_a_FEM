from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patheffects
from mpl_toolkits.mplot3d import axes3d
from itertools import combinations

# attemt to implement
# 3D INTEREST POINT DETECTION USING LOCAL SURFACE CHARACTERISTICS
# WITH APPLICATION IN ACTION RECOGNITION, Michael B. Holte, DOI: 10.1109/ICIP.2014.7026160

# Parameters for corner detection
tau_n = np.pi/9
tau_phi = np.pi/2.1

class boundaryVertex(object):
    """data class that stores boundary vertices and additional informations"""
    def __init__(self, boundaryindex, globalindex, vertex):
        self.boundaryindex = boundaryindex
        self.globalindex = globalindex
        self.vertex = vertex
        self.DoN = None

def calc_vertex_avg_normal_vector(vertexfaces):
    """calculate the average vertex normal vector from the surrounding face normal vectors 
    (normals calculation only works with original mesh faces and not with boundary faces)"""
    avgnormalvector = np.array([face.normal().array() for face in vertexfaces]).mean(axis=0)
    avgnormalvector = avgnormalvector / np.linalg.norm(avgnormalvector)
    return avgnormalvector

def get_vertex_faces(bvertex, meshfacesdict, fb_to_f):
    """get the faces from the original mesh of which the boundary vertex is part of"""
    vertexfaces = [meshfacesdict[fb_to_f[face.index()]] for face in faces(bvertex)]
    return vertexfaces

def get_connected_vertices(bvertex, bverticesdict):
    """get the vertices on the boundary which are connected to the input vertex"""
    edge_entity_index_list = [edge.entities(0) for edge in edges(bvertex)]
    connected_vertices_index_set = set(np.array(edge_entity_index_list).flatten())
    connected_vertices_index_set.remove(bvertex.index())
    connected_vertices = [bverticesdict[idx] for idx in connected_vertices_index_set]
    return connected_vertices

def check_curvature(DoNij_list, tau_n):
    """check the DoNij against tau_n. If more than one DoNij is grater than tau_n, the vertex is
    considered to have a significant change in curvature
    returns the index of the DoNij if more than one DoNij is grater than tau_n, else None"""
    idx_list = []
    for idx, DoNij in enumerate(DoNij_list):
        if DoNij > tau_n:
            idx_list.append(idx)
    if len(idx_list) > 1:
        return idx_list
    else:
        return None

def neighbours_to_neighbourObjs(neighbours, bvertexObjs):
    """converts neighbours (vertices) to neighbourObjs (boundaryVertex class)"""
    neighbourObjs = []
    for neighbour in neighbours:
        for bvertexObj in bvertexObjs:
            if neighbour.index() == bvertexObj.boundaryindex:
                neighbourObjs.append(bvertexObj)
    return neighbourObjs

def get_neighbours(vertex, bverticesdict, depth=1):
    """returns the vertex itself and the neighbours of the vertex with specified depth. 
    depth = 0 returns the vertex itself, depth = 1 returns the vertex and all vertices directly surrounding the given vertex and so on.
    TODO: optimize this function, needs too much computation for high depths"""
    """Test function:
        neigs = get_neighbours(cornervertex_candidate.vertex, bverticesdict, depth=7)
        plot(bmesh)
        for neig in neigs:
            ax.plot(neig.point().x(), neig.point().y(), neig.point().z(), "ro", zorder=3)
        plt.show()
    """
    all_neighbours = get_connected_vertices(vertex, bverticesdict)
    all_neighbours.append(vertex)
    initial_neighbours = get_connected_vertices(vertex, bverticesdict)
    if depth == 0:
        return [vertex]
    elif depth == 1:
        return all_neighbours
    else:
        for neighbour in initial_neighbours:
            neighbours = get_neighbours(neighbour, bverticesdict, depth-1)
            all_neighbours.extend(neighbours)
        return list(set(all_neighbours))

def  write_output(cornervertices, outputfile):
    with open(outputfile, "w") as f:
        f.write("LocalIndex, GlobalIndex, nodex, nodey, nodez, ux, uy, uz\n")
        for idx, vertex in enumerate(cornervertices):
            f.write(f"{idx}, {vertex.globalindex}, {vertex.vertex.point().x()}, {vertex.vertex.point().y()}, {vertex.vertex.point().z()}, 0, 0, 0\n")

def cornerdetection(mesh, outputfile = "output.txt", mode = "angleSum", max_search_depth = 2):
    """detect corners of given mesh, plot corners and output corners with coordinates and wanted displacements to file. 
    mode: set mode for non maximum suppression: either "DoN" or "angleSum" 
    max_search_depth: set depth in which the maximum interest point should be calculated"""
    # get boundary mesh
    bmesh = BoundaryMesh(mesh, "exterior")

    # maps between boundary and normal mesh
    vb_to_v = bmesh.entity_map(0)
    fb_to_f = bmesh.entity_map(2)

    # lookup dicts for boundary and normal faces and vertices
    meshfacesdict = dict((face.index(), face) for face in faces(mesh))
    bfacesdict = dict((face.index(), face) for face in faces(bmesh))
    bverticesdict = dict((vertex.index(), vertex) for vertex in vertices(bmesh))

    cornervertices_before_suppression = []
    bvertexObjs = []
    for bvertex in vertices(bmesh):
        # create boundary vertex object
        bvertexindex = bvertex.index()
        vertexindex = vb_to_v[bvertexindex]
        bvertexObj = boundaryVertex(bvertexindex, vertexindex, bvertex)

        # calculate normal vectors of connected vertices
        connected_vertices = get_connected_vertices(bvertex, bverticesdict)
        connected_vertices_normal_vectors = [calc_vertex_avg_normal_vector(get_vertex_faces(bvertex, meshfacesdict, fb_to_f)) for bvertex in connected_vertices]

        # calculate normal vector of current vertex
        vertexfaces = get_vertex_faces(bvertex, meshfacesdict, fb_to_f)
        avgnormalvector = calc_vertex_avg_normal_vector(vertexfaces)

        # calculate DoNij and DoN for current vertex (see paper for details), modified to calculate normalized DoN
        DoNij = [np.arccos(avgnormalvector.dot(connected_normal_vector)) for connected_normal_vector in connected_vertices_normal_vectors]
        DoN = sum(DoNij)/len(DoNij)
        bvertexObj.DoN = DoN

        # check if corner satisfies the curvature condition (see paper for details)
        idx_list = check_curvature(DoNij, tau_n)
        # if any corneres satisfy the condition calculate the angels between all vectors a, b
        if idx_list:
            combos = combinations(idx_list, 2)
            phis = []
            for combo in combos:
                a = connected_vertices[combo[0]].point() - bvertex.point()
                b = connected_vertices[combo[1]].point() - bvertex.point()
                phi = np.arccos(a.dot(b)/(a.norm()*b.norm()))
                phis.append(phi)
            # check if 3 of the angles satisfy the angle condition 
            # (because the paper doesn't specify how much should satisfy it, and 3 seemed to produce good results)
            relaxed_corner_condition_list = []
            for phi in phis:
                if tau_phi < phi < np.pi - tau_phi or np.pi + tau_phi < phi < 2*np.pi - tau_phi:
                    relaxed_corner_condition_list.append(phi)
            if len(relaxed_corner_condition_list) > 3:
                # calculate average angles, for second mode of non maximum supression
                bvertexObj.avgphi = sum(phis)/len(idx_list)
                cornervertices_before_suppression.append(bvertexObj)
        bvertexObjs.append(bvertexObj)

    # supress non maximum interest points with chosen parameters
    cornervertices_after_suppression = []
    for cornervertex_candidate in cornervertices_before_suppression:
        neighbours = get_neighbours(cornervertex_candidate.vertex, bverticesdict, max_search_depth)
        neighbourObjs = neighbours_to_neighbourObjs(neighbours, bvertexObjs)

        chosen_cornervertexObj = cornervertex_candidate
        if mode == "DoN":
            for neighbourObj in neighbourObjs:
                if neighbourObj.DoN > chosen_cornervertexObj.DoN:
                    chosen_cornervertexObj = neighbourObj
        elif mode == "angleSum":
            for neighbourObj in neighbourObjs:
                if hasattr(neighbourObj, "avgphi"):
                    if neighbourObj.avgphi < chosen_cornervertexObj.avgphi:
                        chosen_cornervertexObj = neighbourObj
        else:
            print("mode not supported")
        if chosen_cornervertexObj == cornervertex_candidate:
            cornervertices_after_suppression.append(chosen_cornervertexObj)

    # write output
    write_output(cornervertices_after_suppression, outputfile)

    # Plotting
    ax = plt.gca(projection='3d')
    # set correct aspect ratio (https://stackoverflow.com/a/64487277)
    ax.set_box_aspect((np.ptp([vertex.vertex.point().x() for vertex in bvertexObjs]), np.ptp([vertex.vertex.point().y() for vertex in bvertexObjs]), np.ptp([vertex.vertex.point().z() for vertex in bvertexObjs])))
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    # Plot boundary mesh
    plot(bmesh)
    # Plot detected corners and labels
    for idx, cornervertex in enumerate(cornervertices_after_suppression):
        ax.plot(cornervertex.vertex.point().x(), cornervertex.vertex.point().y(), cornervertex.vertex.point().z(), "ro", zorder=3)
        txt = ax.text(cornervertex.vertex.point().x(), cornervertex.vertex.point().y(), cornervertex.vertex.point().z(), idx, zorder=100)
        txt.set_path_effects([patheffects.withStroke(linewidth=2, foreground='w')])
    plt.show()

def main():
    mesh = UnitCubeMesh(10, 10, 10)

    mesh = Mesh()
    with XDMFFile("LProfilTetMsh.xdmf") as f:
        f.read(mesh)
    cornerdetection(mesh, "testoutput.txt", mode="angleSum", max_search_depth=2)

#TODO: code cleanup, optimize get_neighbours()

# call main function if profiling is false, else profile the code for performance analysis
profiling = False
if __name__ == "__main__":
    if profiling:
        import cProfile
        import pstats
        profile = cProfile.Profile()
        profile.runcall(main)
        ps = pstats.Stats(profile)
        #ps.print_stats()
        ps.sort_stats('ncalls').print_stats(20)
    else:
        main()