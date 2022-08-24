
import numpy as np

from hexasphere.geometry import Projection, phi


class GnomonicProj(Projection):
    """
    A simple projection that preserves angles and shapes,
    but doesn't preserve area
    """

    def project(self, X, face):
        """
        Project X onto face (k,e1)
        """

        k = self.base_poly.k[face]
        eB = self.base_poly.eB[face]

        P = eB.dot(X)
        P = self.base_poly.FtoC * P / k.dot(X)

        return P

    def inv_project(self, P, face):
        """
        Project P from face to sphere
        """

        k = self.base_poly.k[face]
        eB = self.base_poly.eB[face]

        X = eB.T.dot(P) + self.base_poly.FtoC * k

        return X / np.linalg.norm(X)


class SnyderEAProj(Projection):
    """
    A projection that preserves areas, but (slightly) deforms shapes
    """

    def __init__(self, grid=None):
        super().__init__(grid)

        v0 = self.base_poly.a[0]
        v1 = self.base_poly.b[0]
        v1 = (v0 + v1) * self.base_poly.VtoC / (2 * phi)
        v2 = self.base_poly.k[0]

        self.V = np.linalg.det(np.stack([v0, v1, v2]))

    def project(self, X, face):
        """
        Project X onto face (k,e1)
        """

        abc = self.base_poly.abc[face]
        dist_to_V = np.argsort(abc.dot(X))
        _, v1, v0 = abc[dist_to_V]
        v1 = (v0 + v1) * self.base_poly.VtoC / (2 * phi)
        v2 = self.base_poly.k[face]

        if np.all(X == v0):
            X_P = X
        else:
            if (dist_to_V[1] - dist_to_V[2]) % 3 == 1:
                # np.cross(v0.T, v1.T).dot(v2) >= 0:
                K = self.find_EA_barycenter(X, v0, v1, v2)
                subface = np.stack(
                    [
                        self.base_poly.VtoC * v0,
                        phi * v1,
                        self.base_poly.FtoC * v2
                    ]
                )

            else:
                K = self.find_EA_barycenter(X, v0, v2, v1)
                subface = np.stack(
                    [
                        self.base_poly.VtoC * v0,
                        self.base_poly.FtoC * v2,
                        phi * v1
                    ]
                )
            X_P = subface.T.dot(K)

        eB = self.base_poly.eB[face]
        P = eB.dot(X_P)

        return P

    def find_EA_barycenter(self, X, v0, v1, v2):
        d = self.V * X - np.linalg.det(np.stack([X, v1, v2])) * v0
        d /= np.linalg.norm(d)
        h = np.sqrt((1 - v0.dot(X)) / (1 - v0.dot(d)))
        A = 2 * np.arctan(
            np.linalg.det(np.stack([v0, v1, d]))
            / (1 + v0.dot(v1) + v1.dot(d) + v0.dot(d))
        )
        A2 = np.pi / 30

        K = np.zeros(3)
        K[2] = h * A / A2
        K[1] = h - K[2]
        K[0] = 1 - h

        return K

    def inv_project(self, P, face):
        """
        Project P from face to sphere
        """

        abc = self.base_poly.abc[face]
        dist_to_V = np.argsort(-self.base_poly.Bis.T.dot(P))
        _, v1, v0 = abc[dist_to_V]
        v1 = (v0 + v1) * self.base_poly.VtoC / (2 * phi)
        v2 = self.base_poly.k[face]

        k = self.base_poly.k[face]
        eB = self.base_poly.eB[face]

        X_P = eB.T.dot(P) + self.base_poly.FtoC * k

        if np.all(X_P == self.base_poly.VtoC * v0):
            return X_P

        if (dist_to_V[1] - dist_to_V[2]) % 3 == 1:
            subface = np.stack(
                [self.base_poly.VtoC * v0, phi * v1, self.base_poly.FtoC * v2]
            )

        else:
            v1, v2 = np.copy(v2), np.copy(v1)
            subface = np.stack(
                [self.base_poly.VtoC * v0, self.base_poly.FtoC * v1, phi * v2]
            )

        K = np.linalg.solve(subface.T, X_P)

        c01 = v0.dot(v1)
        c12 = v1.dot(v2)
        c20 = v2.dot(v0)
        s = np.sqrt(1 - c12**2)

        h = 1 - K[0]
        A = (K[2] / h) * np.pi / 30
        S = np.sin(A)
        C = 1 - np.cos(A)
        f = S * self.V + C * (c01 * c12 - c20)
        g = C * s * (1 + c01)
        q = 2 * np.arctan2(g, f) / np.arccos(c12)
        d = self.slerp(v1, v2, q)
        t = np.arccos(1 + h**2 * (v0.dot(d) - 1)) / np.arccos(v0.dot(d))
        X = self.slerp(v0, d, t)

        return X

    def slerp(self, u, v, t):

        ang_dist = np.arccos(u.dot(v))
        return (
            np.sin((1 - t) * ang_dist) * u / np.sin(ang_dist)
            + np.sin(t * ang_dist) * v / np.sin(ang_dist)
        )
