##! [unwarped coordinates]
import numpy as np
from sympy.integrals.quadrature import gauss_legendre
import matplotlib.pyplot as plt

degree = 3 # degree p of the polynomials
row_size = degree + 1 # computationally, it is usually more convenient to work in terms of the size of each row of quadrature points
nodes, weights = gauss_legendre(row_size, 20) # fetch the Gauss-Legendre nodes
nodes = (np.array(nodes).astype(np.float64) + 1)/2 # map from [-1, 1] (the standard convention) to [0, 1] (the Hexed convention)

def physical_qpoint_coords(vertex_coords):
    """ compute the physical coordinates of the element quadrature points from the coordinates of the element vertices """
    phys_coords = np.zeros((3, row_size**3))
    for i in range(row_size):
        for j in range(row_size):
            for k in range(row_size):
                # start with the coordinates of all vertices, and then perform linear interpolation one dimension at a time
                qpoint_coords = vertex_coords + 0
                for i_dimension in range(3):
                    node = nodes[(i, j, k)[i_dimension]]
                    # perform linear interpolation along `i_dimension`, reducing the number of points left to interpolate by a factor of 2
                    qpoint_coords = (1 - node)*qpoint_coords[:qpoint_coords.shape[0]//2, :] + node*qpoint_coords[qpoint_coords.shape[0]//2:, :]
                # having performed linear interpolation three times (trilinear interpolation),
                # we have now gone from 8 points to just 1 point with 3 values (the 3 coordinates)
                phys_coords[:, ((i*row_size) + j)*row_size + k] = qpoint_coords
    return phys_coords

# vertices for an example element which is a frustrum of a pyramid (the top/+z face is shrunk by 0.2 in the x and y directions)
example_vertices = np.array([
    [0., 0., 0.],
    [.2, .2, 1.],
    [0., 1., 0.],
    [.2, .8, 1.],
    [1., 0., 0.],
    [.8, .2, 1.],
    [1., 1., 0.],
    [.8, .8, 1.],
])
# plot the vertices
ax = plt.gcf().add_subplot(projection = "3d")
ax.scatter(example_vertices[:, 0], example_vertices[:, 1], example_vertices[:, 2])
# compute the quadrature point positions of this example element
example_qpoints = physical_qpoint_coords(example_vertices)
# plot the quadrature points
ax.scatter(example_qpoints[0, :], example_qpoints[1, :], example_qpoints[2, :])
ax.set_xlabel("$x_0$")
ax.set_ylabel("$x_1$")
ax.set_zlabel("$x_2$")
plt.show()
##! [unwarped coordinates]

##! [warped coordinates]
def warped_qpoint_coords(vertex_coords, face_warping):
    """ computes the physical quadrature point coordinates again, but this time with face warping """
    unwarped_coords = physical_qpoint_coords(vertex_coords) # start with the unwarped coordinates; face warping is a perturbation to these
    warped_coords = unwarped_coords + 0 # will adjust these as we go along
    for i_dimension in range(3): # start with warping for -/+ x faces, then y, etc
        # the interior quadrature points are a 3D array, whereas the face quadrature points are a 2D array
        # we need to compute some strides in order to map face quadrature points to the corresponding 1D slice of interior quadrature points
        interior_stride = row_size**(2 - i_dimension) # stride between two interior quadrature points in the interior slice
        # stride between two interior quadrature points corresponding to consecutive rows of face quadrature points
        if (i_dimension == 0): face_row_stride = row_size
        else: face_row_stride = row_size**2
        # stride between two interior quadrature points corresponding to consecutive columns of face quadrature points
        if (i_dimension == 2): face_col_stride = row_size
        else: face_col_stride = 1
        # loop through face quadrature points (do - and + faces for the same dimension simultaneously)
        for face_row in range(row_size):
            for face_col in range(row_size):
                # compute the slice of interior quadrature points
                slice_start = face_row*face_row_stride + face_col*face_col_stride
                interior_slice = slice(slice_start, slice_start + row_size*interior_stride, interior_stride)
                interior_coords = unwarped_coords[:, interior_slice]
                # extrapolate the unwarped coordinates of interior quadrature points to faces
                unwarped_face_coords = np.zeros((3, 2)) # column 0 is - face, column 1 is + face
                for i in range(row_size):
                    for face_sign in [0, 1]: # face_sign == 0 indicates negative-facing face whereas face_sign == 1 indicates positive-facing
                        unwarped_face_coords[:, face_sign] += interior_coords[:, i]*np.array([(face_sign - nodes[j])/(nodes[i] - nodes[j]) for j in range(row_size) if j != i]).prod()
                interp_coefs = [1 - nodes[np.newaxis, :], nodes[np.newaxis, :]]
                # compute the difference between coordinates of corresponding quadrature points on opposite faces (\xi_+ - \xi_-)
                diff = unwarped_face_coords[:, [1]] - unwarped_face_coords[:, [0]]
                # interpolate back to interior quadrature points
                for face_sign in [0, 1]:
                    warped_coords[:, interior_slice] += interp_coefs[face_sign]*face_warping[i_dimension, face_sign, face_row*row_size + face_col]*diff
    return warped_coords

# compute an example warping function which is a simple quadratic warping of the -y face
example_warping = np.zeros((3, 2, row_size**2)) # first index is face dimension, second index is face sign, last index is the quadrature point of the face
example_warping[1, 0, :] = 1.5*(nodes[:, np.newaxis]*(1 - nodes[:, np.newaxis])*np.ones(row_size)).flatten()
# plot the vertices
ax = plt.gcf().add_subplot(projection = "3d")
ax.scatter(example_vertices[:, 0], example_vertices[:, 1], example_vertices[:, 2])
# compute the quadrature point positions of this example element
example_qpoints = warped_qpoint_coords(example_vertices, example_warping)
# plot the quadrature points
ax.scatter(example_qpoints[0, :], example_qpoints[1, :], example_qpoints[2, :])
ax.set_xlabel("$x_0$")
ax.set_ylabel("$x_1$")
ax.set_zlabel("$x_2$")
plt.show()
##! [warped coordinates]
