from unittest import TestCase

from bhr.enums import BoundaryCondition
from bhr.utilities import coth, inch_to_m, set_boundary_condition_enum, smoothing_function


class TestUtilities(TestCase):
    def test_set_boundary_condition_enum(self):
        self.assertIs(
            set_boundary_condition_enum("uniform_heat_flux"),
            BoundaryCondition.UNIFORM_HEAT_FLUX,
        )
        with self.assertRaisesRegex(ValueError, "Invalid boundary condition"):
            set_boundary_condition_enum("unsupported")

    def test_inch_to_m(self):
        self.assertAlmostEqual(inch_to_m(1), 0.0254)

    def test_smoothing_function_clamps_outside_range(self):
        self.assertEqual(smoothing_function(-1, 0, 10, 20, 30), 20)
        self.assertEqual(smoothing_function(11, 0, 10, 20, 30), 30)

    def test_coth(self):
        self.assertAlmostEqual(coth(10), 1.0000000041223074)
        self.assertAlmostEqual(coth(1000), 1.0)
