
import numpy as np

# Radius of earth
R = 6371

# Golden Number
phi = (1 + np.sqrt(5)) / 2


def latlon_to_X(lat, lon):
    """
    Converts latlon cordinates into orthogonal coordinates
    """
    lat *= np.pi / 180
    lon *= np.pi / 180

    X = np.zeros(3)
    X[0] = np.cos(lat) * np.cos(lon)
    X[1] = np.cos(lat) * np.sin(lon)
    X[2] = np.sin(lat)

    return X


def X_to_latlon(X):
    """
    Converts orthogonal coordinates into latlon coordinates
    """
    PX = np.copy(X)
    PX[2] = 0

    lat = np.arctan2(X[2], np.linalg.norm(PX))
    if lat == 90 or lat == -90:
        lon = 0
    else:
        lon = np.arctan2(X[1], X[0])

    lat *= 180 / np.pi
    lon *= 180 / np.pi

    return [lat, lon]


def compute_dist(X1, X2, in_latlon=False):
    """
    Computes (spherical) distance between UNITARY vectors X1 and X2
    """
    if in_latlon:
        X1 = latlon_to_X(*X1)
        X2 = latlon_to_X(*X2)
    return np.arccos(X1.dot(X2)) * R


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
