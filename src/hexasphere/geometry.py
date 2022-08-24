
import numpy as np

# Radius of earth
R = 6371

# Golden Number
phi = (1 + np.sqrt(5)) / 2


class Icosahedron:
    def __init__(self, *args, **kwargs):

        # Face-to-center distance
        self.FtoC = np.sqrt(phi**2 - 1 / 3)

        # Vertex-to-center distance
        self.VtoC = np.sqrt(1 + phi**2)

        # Normal vector of all face triangles
        self.k = (np.sqrt(3) / (3 * (1 + phi))) * np.array(
            [
                [2 * phi + 1, phi, 0],
                [1 + phi, 1 + phi, 1 + phi],
                [phi, 0, 2 * phi + 1],
                [1 + phi, -(1 + phi), 1 + phi],
                [2 * phi + 1, -phi, 0],
                [1 + phi, 1 + phi, -(1 + phi)],
                [0, 2 * phi + 1, phi],
                [-phi, 0, 2 * phi + 1],
                [0, -(2 * phi + 1), phi],
                [1 + phi, -(1 + phi), -(1 + phi)],
                [-(2 * phi + 1), phi, 0],
                [-(1 + phi), 1 + phi, -(1 + phi)],
                [-phi, 0, -(2 * phi + 1)],
                [-(1 + phi), -(1 + phi), -(1 + phi)],
                [-(2 * phi + 1), -phi, 0],
                [-(1 + phi), 1 + phi, 1 + phi],
                [0, 2 * phi + 1, -phi],
                [phi, 0, -(2 * phi + 1)],
                [0, -(2 * phi + 1), -phi],
                [-(1 + phi), -(1 + phi), 1 + phi],
            ]
        )

        # Directing summit of all face triangles
        self.a = np.array(
            [
                [phi, 0, -1],
                [1, phi, 0],
                [0, 1, phi],
                [0, -1, phi],
                [1, -phi, 0],
                [1, phi, 0],
                [0, 1, phi],
                [0, -1, phi],
                [1, -phi, 0],
                [phi, 0, -1],
                [-phi, 0, 1],
                [-1, phi, 0],
                [0, 1, -phi],
                [0, -1, -phi],
                [-1, -phi, 0],
                [-1, phi, 0],
                [0, 1, -phi],
                [0, -1, -phi],
                [-1, -phi, 0],
                [-phi, 0, 1],
            ]
        )

        # Directing edge of all face triangles
        self.e1 = np.array(
            [
                [1 - phi, phi, 1],
                [-1, 1 - phi, phi],
                [0, -2, 0],
                [1, 1 - phi, -phi],
                [phi - 1, phi, -1],
                [1 - phi, phi, 1],
                [-1, 1 - phi, phi],
                [0, -2, 0],
                [1, 1 - phi, -phi],
                [phi - 1, phi, -1],
                [phi - 1, phi, -1],
                [1, 1 - phi, -phi],
                [0, -2, 0],
                [-1, 1 - phi, phi],
                [1 - phi, phi, 1],
                [phi - 1, phi, -1],
                [1, 1 - phi, -phi],
                [0, -2, 0],
                [-1, 1 - phi, phi],
                [1 - phi, phi, 1],
            ]
        )

        self.e1[5:10] *= -1
        self.e1[15:20] *= -1

        self.b = self.a + self.e1

        self.e1 /= 2
        self.e2 = np.cross(self.k, self.e1)
        self.eB = [np.stack([self.e1[f], self.e2[f]]) for f in range(20)]

        self.c = 0.5 * (self.a + self.b) + np.sqrt(3) * self.e2

        self.a /= self.VtoC
        self.b /= self.VtoC
        self.c /= self.VtoC

        self.abc = [
            np.stack([self.a[f], self.b[f], self.c[f]]) for f in range(20)
        ]

        # Oriented edges of the face triangle in the face coordinate system
        self.Tr = np.array(
            [
                [-1 / 2, -1 / 2, 1],
                [np.sqrt(3) / 2, -np.sqrt(3) / 2, 0],
            ]
        )

        # Bisectors of the face triangle in the face coordinate system
        self.Bis = np.array(
            [
                [np.sqrt(3) / 2, -np.sqrt(3) / 2, 0],
                [1 / 2, 1 / 2, -1],
            ]
        )

        # Neighboring faces
        self.neighboring_face = np.array(
            [
                [1, 4, 5],
                [2, 0, 6],
                [3, 1, 7],
                [4, 2, 8],
                [0, 3, 9],
                [17, 16, 0],
                [16, 15, 1],
                [15, 19, 2],
                [19, 18, 3],
                [18, 17, 4],
                [11, 14, 15],
                [12, 10, 16],
                [13, 11, 17],
                [14, 12, 18],
                [10, 13, 19],
                [7, 6, 10],
                [6, 5, 11],
                [5, 9, 12],
                [9, 8, 13],
                [8, 7, 14],
            ]
        )

        self.projection = Projection()


class Projection:
    def __init__(self, base_poly: Icosahedron = None):
        self.base_poly = base_poly

    def project(self, X, face):
        raise NotImplementedError

    def inv_project(self, P, face):
        raise NotImplementedError
