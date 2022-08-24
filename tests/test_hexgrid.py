import numpy as np

from src.hexasphere import hexgrid, projection
from src.hexasphere.geometry import X_to_latlon, phi

from unittest import TestCase


Xs = [
    np.array([1, 0, 0], dtype=float),
    np.array([phi, 0, 1], dtype=float),
    np.array([0, 0, 1], dtype=float),
    np.array([2 * phi + 1, phi, 0], dtype=float),
]

LATLONs = [
    X_to_latlon(X) for X in Xs
]

CODEs = [
    "E00036-00018-00018",
    "A00036-00036-00000",
    "C00018-00018-00036",
    "A00024-00024-00024",
]


class TestEncode(TestCase):

    def test_encode_snyder(self):
        grid = hexgrid.HexGrid()
        proj = projection.SnyderEAProj(grid)
        grid.projection = proj

        n = 35

        for X, hex_code in zip(Xs, CODEs):

            X /= np.linalg.norm(X)

            P = hexgrid.Location(grid)
            P.retrieve_by_projection(X=X)

            H = P.find_hex(n=n)[0]

            self.assertEqual(
                H.to_str_id(),
                hex_code
                )

    def test_encode_gnomonic(self):
        grid = hexgrid.HexGrid()
        proj = projection.GnomonicProj(grid)
        grid.projection = proj

        n = 35

        for X, hex_code in zip(Xs, CODEs):

            X /= np.linalg.norm(X)

            P = hexgrid.Location(grid)
            P.retrieve_by_projection(X=X)

            H = P.find_hex(n=n)[0]

            self.assertEqual(
                H.to_str_id(),
                hex_code
                )


class TestDecode(TestCase):

    def test_decode_snyder(self):
        grid = hexgrid.HexGrid()
        proj = projection.SnyderEAProj(grid)
        grid.projection = proj

        for latlon, hex_code in zip(LATLONs, CODEs):

            lat, lon = grid.hex_to_latlon(hex_code, in_str=True)

            self.assertAlmostEqual(
                lat,
                latlon[0]
                )

            if np.round(lat-90, 7) != 0:

                self.assertAlmostEqual(
                    lon,
                    latlon[1]
                    )

    def test_decode_gnomonic(self):
        grid = hexgrid.HexGrid()
        proj = projection.GnomonicProj(grid)
        grid.projection = proj

        for latlon, hex_code in zip(LATLONs, CODEs):

            lat, lon = grid.hex_to_latlon(hex_code, in_str=True)

            self.assertAlmostEqual(
                lat,
                latlon[0]
                )

            if np.round(lat-90, 7) != 0:

                self.assertAlmostEqual(
                    lon,
                    latlon[1]
                    )
