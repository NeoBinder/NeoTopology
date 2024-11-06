import math

import numpy as np

from ttk.core import unit

RAD_TO_DEG = 180 / math.pi


def parameters_to_box_matrix(a_length, b_length, c_length, alpha, beta, gamma):
    a = np.array([a_length.magnitude, 0, 0])
    b = np.array(
        [b_length.magnitude * math.cos(gamma), b_length.magnitude * math.sin(gamma), 0]
    )
    cx = c_length.magnitude * math.cos(beta)
    cy = (
        c_length.magnitude
        * (math.cos(alpha) - math.cos(beta) * math.cos(gamma))
        / math.sin(gamma)
    )
    cz = math.sqrt(c_length.magnitude * c_length.magnitude - cx * cx - cy * cy)
    c = np.array([cx, cy, cz])

    for i in range(3):
        if abs(a[i]) < 1e-6:
            a[i] = 0.0
        if abs(b[i]) < 1e-6:
            b[i] = 0.0
        if abs(c[i]) < 1e-6:
            c[i] = 0.0

    # Make sure they're in the reduced form required by OpenMM.

    c = c - b * round(c[1] / b[1])
    c = c - a * round(c[0] / a[0])
    b = b - a * round(b[0] / a[0])

    box_matrix = np.eye(4)
    box_matrix[0, :3] = a
    box_matrix[1, :3] = b
    box_matrix[2, :3] = c
    return box_matrix


class UnitCell:

    def __init__(self, matrix):
        self.matrix = matrix

    @classmethod
    def from_openmm(cls, topology):
        from openmm import unit as openmm_unit

        top_PeriodicBoxVectors = topology.getPeriodicBoxVectors()
        if top_PeriodicBoxVectors:
            box_matrix = np.array(
                top_PeriodicBoxVectors.value_in_unit(openmm_unit.nanometer)
            )
        else:
            raise NotImplementedError("openmm topology not defined")
        return cls(box_matrix)

    @classmethod
    def from_positions(cls, positions):
        from openmm import unit as openmm_unit

        if not hasattr(positions, "shape"):
            # if isnot np ndarray
            positions = np.asarray([_pos for _pos in positions])
        positions = positions.value_in_unit(openmm_unit.nanometer)
        box_vec = positions.max(0) - positions.min(0) + 1.0
        box_matrix = np.eye(4)
        box_matrix[:3, :3] *= box_vec

        return cls(box_matrix)

    @classmethod
    def from_parameter(
        cls,
        a_length,
        b_length,
        c_length,
        alpha,
        beta,
        gamma,
        threshold=(1e-6) * unit.nanometer,
    ):
        a_length = (float(a_length) * unit.angstrom).to(unit.nanometer)
        b_length = (float(b_length) * unit.angstrom).to(unit.nanometer)
        c_length = (float(c_length) * unit.angstrom).to(unit.nanometer)
        alpha = float(alpha) * math.pi / 180.0 * unit.radians
        beta = float(beta) * math.pi / 180.0 * unit.radians
        gamma = float(gamma) * math.pi / 180.0 * unit.radians
        if a_length < threshold or b_length < threshold or c_length < threshold:
            raise ValueError(
                "a: {}, b: {},c: {} not meet threshold:{}".format(
                    a_length, b_length, c_length, threshold
                )
            )

        box_matrix = parameters_to_box_matrix(
            a_length, b_length, c_length, alpha, beta, gamma
        )
        return cls(box_matrix)

    @classmethod
    def from_content(cls, content, threshold=(1e-6) * unit.nanometer):
        # unit in angstrom
        if not content.startswith("CRYST1"):
            raise NotImplementedError(
                "{}\n not supported for parsing unit cell".format(content)
            )
        a_length = (float(content[6:15]) * unit.angstrom).to(unit.nanometer)
        b_length = (float(content[15:24]) * unit.angstrom).to(unit.nanometer)
        c_length = (float(content[24:33]) * unit.angstrom).to(unit.nanometer)
        alpha = float(content[33:40]) * math.pi / 180.0 * unit.radians
        beta = float(content[40:47]) * math.pi / 180.0 * unit.radians
        gamma = float(content[47:54]) * math.pi / 180.0 * unit.radians

        if a_length < threshold or b_length < threshold or c_length < threshold:
            raise ValueError(
                "a: {}, b: {},c: {} not meet threshold:{}".format(
                    a_length, b_length, c_length, threshold
                )
            )

        box_matrix = parameters_to_box_matrix(
            a_length, b_length, c_length, alpha, beta, gamma
        )

        return cls(box_matrix)

    def to_openmm(self):
        """Set the vectors defining the periodic box."""
        from openmm import unit as openmm_unit

        vectors = self.matrix[:3, :3] * openmm_unit.nanometer
        return vectors

    def to_vector(self):
        vec = (self.matrix[0][0], self.matrix[1][1], self.matrix[2][2])
        return np.array(vec) * unit.nanometer

    def to_cryst1(self):
        #  a, b, c = self.matrix[0][0], self.matrix[1][1], self.matrix[2][2]
        a, b, c = self.matrix[0], self.matrix[1], self.matrix[2]
        a_length = np.linalg.norm(a)
        b_length = np.linalg.norm(b)
        c_length = np.linalg.norm(c)
        alpha = math.acos(np.dot(b, c) / (b_length * c_length))
        beta = math.acos(np.dot(c, a) / (c_length * a_length))
        gamma = math.acos(np.dot(a, b) / (a_length * b_length))
        return (
            a_length * 10,
            b_length * 10,
            c_length * 10,
            alpha * RAD_TO_DEG,
            beta * RAD_TO_DEG,
            gamma * RAD_TO_DEG,
        )
