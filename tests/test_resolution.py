from unittest import TestCase

from src.hexasphere import hexgrid


class TestResolution(TestCase):

    grid = hexgrid.HexGrid()

    n1 = grid.compute_n_for_height(0.25)
    n2 = grid.compute_n_for_height(2.5)
    n3 = grid.compute_n_for_height(50)

    def test_resolution(self):

        self.assertEqual(self.n1, 15348)
        self.assertEqual(self.n2, 1534)
        self.assertEqual(self.n3, 76)

    def test_compute_height(self):

        self.assertEqual(self.grid.compute_height_for_n(self.n1), 0.2499983483713541)
        self.assertEqual(self.grid.compute_height_for_n(self.n2), 2.4998206183400087)
        self.assertEqual(self.grid.compute_height_for_n(self.n3), 49.83408635262226)

    def test_compute_radius(self):

        self.assertEqual(self.grid.compute_radius_for_n(self.n1), 0.2625170496179913)
        self.assertEqual(self.grid.compute_radius_for_n(self.n2), 2.624999475300683)
        self.assertEqual(self.grid.compute_radius_for_n(self.n3), 52.32953499463051)
