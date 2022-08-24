
import numpy as np

from hexasphere.geometry import Icosahedron, R
from hexasphere.geometry import compute_dist, X_to_latlon, latlon_to_X


class HexGrid(Icosahedron):
    def __init__(self, face_A=None, overlap=0):
        """
        ### Parameters

        - face_A, optional : np.array, shape = (3, 3), dtype = float

        Orthogonal coordinates of the three vertices of face_A

        - overlap : float

        A positive overlap value (in km) will give a grid
        where hexes overlap over the given distance
        """
        super().__init__(face_A)

        self.overlap = overlap
        self.margin = (
            0.5 * overlap * np.sqrt(
                5 * np.sqrt(3) / (np.pi * R**2)
            )
        )

    def compute_n_for_radius(self, r):
        """
        Computes n (the number of hexes an edge of a face goes through),
        so that the average area of hexagon tiles is equal to the area
        of a disk of radius r (in kilometers)
        """
        h_area = np.pi * r**2
        nb_h = (4 * np.pi * R**2) / h_area

        return int(np.round(np.sqrt(nb_h / 10) - 1))

    def compute_radius_for_n(self, n):
        """
        Returns the approximative radius of hexes in a grid of resolution n+1
        The radius of an hex is the radius of the circle with an area equal to
        the area of the given hex
        """

        nb_h = 20 * (n + 1) ** 2 / 2
        h_area = (4 * np.pi * R**2) / nb_h
        eq_r = np.sqrt(h_area / np.pi)

        return eq_r

    def compute_n_for_height(self, h):
        """
        Computes n (the number of hexes an edge of a face goes through),
        so that the average height of hexagon tiles is h (in kilometers)
        """
        h_area = 2 * h**2 * np.sqrt(3)
        nb_h = (4 * np.pi * R**2) / h_area

        return int(np.round(np.sqrt(nb_h / 10) - 1))

    def compute_height_for_n(self, n):
        """
        Returns the approximative height of hexes in a grid of resolution n+1
        The heigh of an hex is the distance between its center and any of its
        edges
        """

        nb_h = 20 * (n + 1) ** 2 / 2
        h_area = (4 * np.pi * R**2) / nb_h
        height = np.sqrt(h_area * np.sqrt(3) / 6)

        return height

    def compute_n_for_side(self, s):
        """
        Computes n (the number of hexes an edge of a face goes through),
        so that the average side of hexagon tiles is s (in kilometers)
        """
        h_area = 3 * s**2 * np.sqrt(3) / 2
        nb_h = (4 * np.pi * R**2) / h_area

        return int(np.round(np.sqrt(nb_h / 10) - 1))

    def compute_side_for_n(self, n):
        """
        Returns the approximative side length of hexes in a grid of resolution
        n+1
        """

        nb_h = 20 * (n + 1) ** 2 / 2
        h_area = (4 * np.pi * R**2) / nb_h
        side = np.sqrt(2 * h_area * np.sqrt(3) / 9)

        return side

    def rectify_coordinates(self, face, pos, n):
        """
        Retrieves new face and new pos in face, if given pos is out of face
        """

        x, y, z = pos

        while True:

            if x > n + 1:
                if face % 10 < 5:
                    x, y, z = ((n + 1) - z, 2 * (n + 1) - x, (n + 1) - y)
                else:
                    x, y, z = (2 * (n + 1) - x, (n + 1) - y, (n + 1) - z)
                face = self.neighboring_face[face, 0]

            elif y > n + 1:
                if face % 10 < 5:
                    x, y, z = (2 * (n + 1) - y, (n + 1) - z, (n + 1) - x)
                else:
                    x, y, z = ((n + 1) - x, 2 * (n + 1) - y, (n + 1) - z)
                face = self.neighboring_face[face, 1]

            elif z > n + 1:
                x, y, z = ((n + 1) - x, (n + 1) - y, 2 * (n + 1) - z)
                face = self.neighboring_face[face, 2]

            else:
                break

        return face, (x, y, z)

    def latlon_to_hex(self, lat, lon, n, out_str=False):
        """
        Returns hex(es) to which the point (lat, lon) of the sphere belongs

        ## Parameters

        - lat, lon : float

        Latitude and longitude, in radians

        - n : int

        Number of hexagons on the edge of each face for the desired hexagonal
        resolution

        - overlap : float, optional

        Lenght (in km) of the overlapping section of two consecutive hexagons

        - out_str : bool, optional

        If True, then the function returns a string of 18 characters
        representing the hex
        """
        location = Location(self)
        # The point is projected on the icosahedron
        location.retrieve_by_projection(latlon=(lat, lon))
        # The hex to which the projected point belongs is retrieved
        return location.find_hex(n, out_str)

    def hex_to_latlon(self, hexagon, n=None, in_str=False):
        """
        Returns the (lat, lon) coordinates of the center of the hexagon

        ## Parameters

        - hexagon : Hexagon (or str)

        - in_str : bool, optional

        If True, hexagon parameter must be the string id of hexagon
        """
        if in_str:
            hexagon = Hexagon(self, str_id=hexagon, res=n + 1)
        X = self.projection.inv_project(hexagon.get_P(), hexagon.face)
        return X_to_latlon(X)


class Hexagon:

    face_char = [
        "A",
        "B",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "J",
        "K",
        "L",
        "M",
        "N",
        "O",
        "P",
        "Q",
        "R",
        "S",
        "T",
    ]

    def __init__(
        self,
        grid: HexGrid,
        face=None,
        pos=None,
        str_id=None,
        solve_conflicts=False,
        res=None,
    ):
        """
        ## Parameters

        - grid : HexGrid

        - face : int, optional

        - pos : (int, int, int), optional

        - str_id : str, optional

        If provided, then face and pos parameters aren't required and taken
        into account. Instead, they are deduced from str_id

        - solve_conflicts : bool, optional

        If True, a method is run to check which is the standard (face, pos)
        couple (aimed at hexagons on edges or vertices of the icosahedron)

        - res : int, optional
        """

        self.grid = grid

        if str_id is not None:
            face, pos = self.from_str_id(str_id)

        if res is None:
            res = sum(pos) // 2

        self.n = res - 1

        a, b, c = pos

        self.vertex_conflicts = (a == 0) or (b == 0) or (c == 0)
        self.edge_conflicts = (
            (a == self.n + 1) or (b == self.n + 1) or (c == self.n + 1)
        )
        self.edge_conflicts = self.edge_conflicts and not self.vertex_conflicts

        if solve_conflicts and (self.vertex_conflicts or self.edge_conflicts):
            face, pos = self.resolve_conflicts(face, pos)

        self.face = face
        self.pos = pos

        self.P = None

    def from_str_id(self, str_id):
        """
        Retrieve face and pos of an hexagon from a string identifier
        """
        face = self.face_char.index(str_id[0])
        raw_pos = str_id[1:].split("-")
        pos = tuple(int(raw_pos[i]) for i in range(3))
        return face, pos

    def to_str_id(self):
        """
        Returns the string identifier of hexagon
        """
        a, b, c = self.pos
        return self.face_char[self.face] + f"{a:05}-{b:05}-{c:05}"

    def __str__(self):
        return str(self.face) + " " + str(self.pos) + " / n = " + str(self.n)

    def get_P(self):
        """
        P represents the coordinates of the hex center
        in the face coordinate system
        """
        self.P = (
            2 * np.sqrt(3) * self.grid.Bis.dot(np.array(self.pos))
            / (3 * (self.n + 1))
        )
        return self.P

    def resolve_conflicts(self, face, pos):
        """
        If the hexagon falls on an edge or a vertex of the icosahedron,
        we must decide to which face it belongs
        """
        a, b, c = pos
        n = self.n

        if self.vertex_conflicts:

            if a == 0 and face % 10 < 5:
                return face, (a, b, c)

            elif a == 0:
                return (face - 4) % 5 + 10 * (face // 10), (a, b, c)

            elif b == 0 and face % 10 < 5:
                return (face + 1) % 5 + 10 * (face // 10), (b, a, c)

            elif b == 0:
                return face - 5, (b, a, c)

            elif c == 0 and face % 10 < 5:
                return 10 * (face // 10), (a, b, c)

            elif c == 0:
                return (2 - face) % 5 + 10 * (1 - face // 10), (c, a, b)

        if self.edge_conflicts:

            if c == n + 1 and face % 10 < 5:
                return face, (a, b, c)

            elif c == n + 1:
                return face - 5, (n + 1 - a, n + 1 - b, c)

            elif a == n + 1 and (face // 10 == 0 or face % 10 < 5):
                return face, (a, b, c)

            elif a == n + 1:
                return (2 - face) % 5 + 5, (a, n + 1 - b, n + 1 - c)

            elif b == n + 1 and face % 10 < 5:
                return (
                    (face - 1) % 5 + 10 * (face // 10),
                    (b, n + 1 - c, n + 1 - a)
                )

            elif b == n + 1 and face // 10 == 1:
                return (1 - face) % 5 + 5, (n + 1 - a, b, n + 1 - c)

            elif b == n + 1:
                return face, (a, b, c)

    def retrieve_polygon(
        self, overlap=0, out_latlon=False, out_lonlat=False, out_geojson=False
    ):
        """
        Returns the list of vertices of the hexagon, as 3D vectors
        """
        abc = np.array(self.pos)
        n = self.n

        P_V = []
        face_V = []

        if overlap == 0:
            v = (3 * np.eye(3) - np.ones((3, 3))) / 3
        else:
            v = (
                (1 + 0.5 * overlap / self.compute_height())
                * (3 * np.eye(3) - np.ones((3, 3)))
                / 3
            )

        for i in range(3):
            for s in [1, -1]:
                pos = abc + s * v[i]
                face_v, pos_v = self.grid.rectify_coordinates(
                    self.face,
                    pos,
                    n
                )
                P_v = (
                    2 * np.sqrt(3) * self.grid.Bis.dot(np.array(pos_v))
                    / (3 * (n + 1))
                )
                P_V.append(P_v)
                face_V.append(face_v)

        res = [
            self.grid.projection.inv_project(P_V[i], face_V[i])
            for i in [0, 3, 4, 1, 2, 5]
        ]

        if out_latlon:
            res = [X_to_latlon(X) for X in res]
            res.append(res[0])

        if out_lonlat or out_geojson:
            res = [X_to_latlon(X)[::-1] for X in res]
            res.append(res[0])

        if out_geojson:
            res.append(res[0])
            res = {
                "coordinates": [res],
                "type": "Polygon",
            }

        return res

    def compute_neighbor(self, dP=(0, 0, 0)):

        a, b, c = self.pos
        i, j, k = dP
        a, b, c = a + i, b + j, c + k

        face, pos = self.grid.rectify_coordinates(self.face, (a, b, c), self.n)

        return Hexagon(
            self.grid,
            face,
            pos,
            res=self.n + 1,
            solve_conflicts=True
        )

    def effective_radius(self):
        """
        Returns the average distance between the vertices of hex and its center
        """

        X_C = self.grid.projection.inv_project(self.get_P(), self.face)

        X_V = self.retrieve_polygon()

        avg_r = 0
        for X_v in X_V:
            avg_r += compute_dist(X_v, X_C) / len(X_V)

        return avg_r * np.sqrt((3 * np.sqrt(3)) / (2 * np.pi))

    def find_descendant_hexes(self):

        n = self.n
        face = self.face
        a, b, c = self.pos
        root_descendant_hex = Hexagon(
            self.grid, face, (4 * a, 4 * b, 4 * c), res=4 * (n + 1)
        )

        descendants = []

        for i in range(-2, 3):
            for j in range(-2, 3):
                for k in range(-2, 3):
                    if i + j + k == 0:

                        neighbor_hex = root_descendant_hex.compute_neighbor(
                            (i, j, k)
                        )

                        if neighbor_hex is not None:
                            descendants.append(neighbor_hex)

        return descendants

    def find_parent_hex(self, gen=1):

        if gen == 0:
            return [self]

        n = self.n
        face = self.face
        a, b, c = self.pos

        pos = np.stack([a, b, c])
        new_pos = np.zeros(3).astype(int)
        for i in range(3):
            if pos[i] % 4**gen <= 2 ** (2 * gen - 1):
                new_pos[i] = pos[i] // 4**gen
            else:
                new_pos[i] = pos[i] // 4**gen + 1

        delta = pos - (4**gen) * new_pos

        res = []

        if np.sum(delta) > 0:

            to_increment = np.argsort(delta)

            new_pos[to_increment[2]] += 1
            a, b, c = new_pos
            new_pos[to_increment[2]] -= 1
            res.append((a, b, c))

            if delta[to_increment[2]] == delta[to_increment[1]]:
                new_pos[to_increment[1]] += 1
                a, b, c = new_pos
                new_pos[to_increment[1]] -= 1
                res.append((a, b, c))

        elif np.sum(delta) < 0:

            to_decrement = np.argsort(delta)

            new_pos[to_decrement[0]] -= 1
            a, b, c = new_pos
            new_pos[to_decrement[0]] += 1
            res.append((a, b, c))

            if delta[to_decrement[0]] == delta[to_decrement[1]]:
                new_pos[to_decrement[1]] -= 1
                a, b, c = new_pos
                new_pos[to_decrement[1]] += 1
                res.append((a, b, c))

        else:
            res = [tuple(new_pos)]

        return [
            Hexagon(
                self.grid,
                face,
                pos,
                res=(n + 1) // 4**gen,
                solve_conflicts=True
            )
            for pos in res
        ]

    def compute_radius(self):
        return self.grid.compute_radius_for_n(self.n)

    def compute_height(self):
        return self.grid.compute_height_for_n(self.n)

    def compute_side(self):
        return self.grid.compute_side_for_n(self.n)


class Location:
    """
    A (geographic) location to be projected on a base icosahedron of an hexgrid

    ### Attributes

    - self.X : 3D orthogonal coordinates of the location
    - self.latlon : (lat,lon) coordinates of the location
    - self.face : face of the icosahedron the projected location belongs to
    - self.P : face coordinates of the projected location
    """

    def __init__(self, grid: HexGrid):

        self.X = None
        self.latlon = None

        self.grid = grid
        self.face = None
        self.P = None

    def retrieve_by_projection(self, X=None, latlon=None):

        if latlon is not None:
            self.latlon = latlon
            self.X = latlon_to_X(*latlon)
        else:
            self.X = X

        self.face = np.argmax(self.grid.k.dot(self.X))
        self.P = self.grid.projection.project(self.X, self.face)

    def find_pos_from_P_TrB(self, P_TrB, n):
        """
        Finds the hex pos (a, b, c) of a point P converted into triangular
        coordinates P_TrB (orthogonal projection on sides of the face triangle)
        """

        N = 2 * n + 1

        x, y, z = P_TrB

        u, v, w = x * (N + 1) / 2, y * (N + 1) / 2, z * (N + 1) / 2
        u, v, w = int(u), int(v), int(w)

        a = (2 + (N - v) + w) // 3
        b = (2 + (N - w) + u) // 3
        c = N + 1 - (a + b)

        if a < 0 or b < 0 or c < 0 or a > n + 1 or b > n + 1 or c > n + 1:
            raise ValueError("Projected point outside of face")

        return (a, b, c)

    def find_hex(self, n, out_str=False):
        """
        Finds which hexagon(s) the point P of a face belongs to
        """
        margin = self.grid.margin

        x, y, z = self.grid.Tr.T.dot(self.P)
        x, y, z = x + 1, y + 1, z + 1

        if margin > 0:
            s = set()

            for j in [-margin, margin]:
                s.add(self.find_pos_from_P_TrB((x + j, y, z), n))
                s.add(self.find_pos_from_P_TrB((x, y + j, z), n))
                s.add(self.find_pos_from_P_TrB((x, y, z + j), n))

            res = list(s)

        else:
            res = [self.find_pos_from_P_TrB((x, y, z), n)]

        if out_str:
            return [
                Hexagon(
                    self.grid, self.face, pos, res=n + 1, solve_conflicts=True
                ).to_str_id()
                for pos in res
            ]
        else:
            return [
                Hexagon(
                    self.grid,
                    self.face,
                    pos,
                    res=n + 1,
                    solve_conflicts=True
                )
                for pos in res
            ]


if __name__ == "__main__":

    # import geojson

    hexgrid = HexGrid()
    n = hexgrid.compute_n_for_radius(0.25)

    from hexasphere.projection import SnyderEAProj

    snyder_proj = SnyderEAProj(hexgrid)
    hexgrid.projection = snyder_proj

    print("Grid properties :")
    print("radius :", hexgrid.compute_radius_for_n(n))
    print("height :", hexgrid.compute_height_for_n(n))
    print("side :", hexgrid.compute_side_for_n(n))

    for _ in range(10000):
        input("")
        print("__________________")
        lat = np.random.random() * 180 - 90
        lon = np.random.random() * 360 - 180
        print(lat, lon)
        H = hexgrid.latlon_to_hex(lat, lon, n)[0]  # 31.717
        H_ll = hexgrid.hex_to_latlon(H)
        print(H)
        print(H_ll)
        print(
            "distance :",
            compute_dist(
                (lat, lon),
                H_ll,
                in_latlon=True
            )
        )

        """poly = H.retrieve_polygon(out_latlon=True)
        poly = [point[::-1] for point in poly]
        poly = geojson.Polygon([poly])
        poly_json = geojson.dumps(poly)
        with open(H.to_str_id() + ".geojson", mode="w") as f:
            f.write(poly_json)"""

        """descendants_H = H.find_descendant_hexes()
        for d in descendants_H:
            print("-", d)
            poly = d.retrieve_polygon(out_latlon=True)
            poly = [point[::-1] for point in poly]
            poly = geojson.Polygon([poly])
            poly_json = geojson.dumps(poly)
            with open(d.to_str_id() + ".geojson", mode="w") as f:
                f.write(poly_json)"""

        print("hex radius :", H.effective_radius())
