from unittest import TestCase

from src.hexasphere import hexgrid, projection


class TestNeighbors(TestCase):

    def test_compute_neighbor(self):

        grid = hexgrid.HexGrid()
        proj = projection.SnyderEAProj(grid)
        grid.projection = proj

        H = hexgrid.Hexagon(grid, str_id="A00006-00024-00018")

        self.assertEqual(H.compute_neighbor((1, -2, 1)).to_str_id(), "A00007-00022-00019")
        self.assertEqual(H.compute_neighbor((0, 1, -1)).to_str_id(), "E00023-00007-00018")
        self.assertEqual(H.compute_neighbor((18, 0, -18)).to_str_id(), "A00024-00024-00000")

    def test_k_ring(self):

        grid = hexgrid.HexGrid()
        proj = projection.SnyderEAProj(grid)
        grid.projection = proj

        H = hexgrid.Hexagon(grid, str_id="A00001-00024-00023")

        computed_ring = H.k_ring(2, out_str=True)

        true_ring = [
            'J00001-00024-00023',
            'J00002-00023-00023',
            'F00024-00001-00023',
            'A00000-00024-00024',
            'E00023-00001-00024',
            'E00022-00002-00024',
            'F00023-00002-00023',
            'A00001-00023-00024',
            'E00024-00001-00023',
            'E00023-00002-00023',
            'E00022-00003-00023',
            'A00002-00022-00024',
            'A00002-00023-00023',
            'E00024-00002-00022',
            'E00023-00003-00022',
            'A00003-00022-00023',
            'A00003-00023-00022',
            'E00024-00003-00021'
        ]

        for h_id in true_ring:
            self.assertIn(h_id, computed_ring)

        for h_id in computed_ring:
            self.assertIn(h_id, true_ring)
