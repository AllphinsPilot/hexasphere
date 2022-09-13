from unittest import TestCase

from src.hexasphere import hexgrid, projection


class TestOverlapping(TestCase):

    def test_overlapping(self):

        grid = hexgrid.HexGrid()
        proj = projection.SnyderEAProj(grid)
        grid.projection = proj

        n = 35
        height = grid.compute_height_for_n(n)
        grid.set_overlap(2.000001*height)

        hexes = grid.latlon_to_hex(0, 0, n, out_str=True)

        self.assertCountEqual(
            hexes,
            [
                "E00036-00018-00018",
                "E00036-00017-00019",
                "E00036-00019-00017",
                "E00035-00019-00018",
                "E00035-00018-00019",
                "A00019-00035-00018",
                "A00018-00035-00019",
            ]
        )

    def test_no_overlapping(self):

        grid = hexgrid.HexGrid()
        proj = projection.SnyderEAProj(grid)
        grid.projection = proj

        n = 35
        height = grid.compute_height_for_n(n)
        grid.set_overlap(1.999999*height)

        hexes = grid.latlon_to_hex(0, 0, n, out_str=True)

        self.assertCountEqual(
            hexes,
            [
                "E00036-00018-00018"
            ]
        )