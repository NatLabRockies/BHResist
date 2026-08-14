from unittest import TestCase

from bhr.utilities import coth


class TestUtilities(TestCase):
    def test_coth(self):
        self.assertAlmostEqual(coth(10), 1.0000000041223074)
        self.assertAlmostEqual(coth(1000), 1.0)
