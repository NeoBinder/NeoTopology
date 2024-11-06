import numpy as np


def apply_matrix_1d(vec, matrix):
    vec_dim = vec.shape[0]
    if vec_dim == 3:
        vec = np.asarray([*vec, 1])
    return matrix.dot(vec)[:, :vec_dim]


def apply_matrix_2d(vecs, matrix):
    vecs_dim = vecs.shape[1]
    if vecs_dim == 3:
        vecs = np.hstack([vecs, np.ones(vecs.shape[0]).reshape(-1, 1)])
    vecs = matrix.dot(vecs.T).T
    return vecs[:, :vecs_dim]


def rigid_transform_3D(A, B):
    assert A.shape == B.shape

    num_rows, num_cols = A.shape
    if num_rows != 3:
        raise Exception(f"matrix A is not 3xN, it is {num_rows}x{num_cols}")

    num_rows, num_cols = B.shape
    if num_rows != 3:
        raise Exception(f"matrix B is not 3xN, it is {num_rows}x{num_cols}")

    # find mean column wise
    centroid_A = np.mean(A, axis=1)
    centroid_B = np.mean(B, axis=1)

    # ensure centroids are 3x1
    centroid_A = centroid_A.reshape(-1, 1)
    centroid_B = centroid_B.reshape(-1, 1)

    # subtract mean
    Am = A - centroid_A
    Bm = B - centroid_B

    H = Am @ np.transpose(Bm)

    # sanity check
    # if linalg.matrix_rank(H) < 3:
    #    raise ValueError("rank of H = {}, expecting 3".format(linalg.matrix_rank(H)))

    # find rotation
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T

    # special reflection case
    if np.linalg.det(R) < 0:
        print("det(R) < R, reflection detected!, correcting for it ...")
        Vt[2, :] *= -1
        R = Vt.T @ U.T

    t = -R @ centroid_A + centroid_B

    matrix = np.eye(4)
    matrix[:3, :3] = R
    matrix[:3, 3] = t.reshape(-1)

    return matrix


class RotationMatrix:

    def __init__(self):
        self.matrix = np.eye(4)

    def __repr__(self):
        content = super().__repr__()
        content = content + "\nwith matrix:\n" + str(self.matrix)
        return content

    @classmethod
    def from_matrix(cls, matrix):
        rmatrix = cls()
        if matrix.shape[0] == 3 and matrix.shape[1] == 3:
            rmatrix.matrix[:3, :3] = matrix
        if matrix.shape[0] == 3 and matrix.shape[1] == 4:
            rmatrix.matrix[:3] = matrix
        else:
            rmatrix.matrix = matrix
        return rmatrix

    @classmethod
    def from_bivec(cls, vec1, vec2):
        # untested
        vec1 /= np.linalg.norm(vec1)
        vec2 /= np.linalg.norm(vec2)

        vec_z = np.cross(vec1, vec2)
        vec_z /= np.linalg.norm(vec_z)

        vec_y = np.cross(vec_z, vec1)
        vec_y /= np.linalg.norm(vec_y)

        rotm = cls()
        rotm.matrix[0, :3] = vec1
        rotm.matrix[1, :3] = vec_y
        rotm.matrix[2, :3] = vec_z
        return rotm

    def apply(self, vecs):
        if len(vecs.shape) == 2:
            vecs = apply_matrix_2d(vecs, self.matrix)
        elif len(vecs.shape) == 1:
            vecs = apply_matrix_1d(vecs, self.matrix)
        return vecs

    @staticmethod
    def generate_rotation_matrix(rmatrix, rotation_origin):
        matrix = np.eye(4)
        matrix[:3, 3] = -rotation_origin
        matrix = rmatrix.dot(matrix)
        matrix[:3, 3] += rotation_origin
        return matrix
