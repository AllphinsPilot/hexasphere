from unittest import TestCase

from src.hexasphere import hexgrid, projection


class TestPolygon(TestCase):

    def test_polygon(self):

        grid = hexgrid.HexGrid()
        proj = projection.SnyderEAProj(grid)
        grid.projection = proj

        H = hexgrid.Hexagon(grid, str_id="A00012-00012-00012")
        h1 = hexgrid.Hexagon(grid, str_id="A00013-00011-00012")
        h2 = hexgrid.Hexagon(grid, str_id="A00012-00011-00013")
        h3 = hexgrid.Hexagon(grid, str_id="A00011-00012-00013")
        h4 = hexgrid.Hexagon(grid, str_id="A00011-00013-00012")
        h5 = hexgrid.Hexagon(grid, str_id="A00012-00013-00011")
        h6 = hexgrid.Hexagon(grid, str_id="A00013-00012-00011")

        v1, v2, v3, v4, v5, v6, _ = H.retrieve_polygon(out_latlon=True)

        self.assertAlmostEqual(v1, h1.retrieve_polygon(out_latlon=True)[4])
        self.assertAlmostEqual(v2, h2.retrieve_polygon(out_latlon=True)[5])
        self.assertAlmostEqual(v3, h3.retrieve_polygon(out_latlon=True)[0])
        self.assertAlmostEqual(v4, h4.retrieve_polygon(out_latlon=True)[1])
        self.assertAlmostEqual(v5, h5.retrieve_polygon(out_latlon=True)[2])
        self.assertAlmostEqual(v6, h6.retrieve_polygon(out_latlon=True)[3])
