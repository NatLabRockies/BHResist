import unittest
from math import log

from scp import get_fluid

from bhr.pipe import Pipe


class TestPipe(unittest.TestCase):
    def setUp(self):
        self.inputs = {
            "pipe_outer_diameter": 0.0334,
            "pipe_dimension_ratio": 11,
            "pipe_length": 100,
            "pipe_conductivity": 0.4,
            "fluid_type": "WATER",
        }

    def test_init_pipe(self):
        p = Pipe(**self.inputs)
        tol = 1e-4

        # props
        self.assertAlmostEqual(p.pipe_conductivity, 0.4, delta=tol)

        # geometry
        self.assertAlmostEqual(p.pipe_outer_diameter, 0.0334, delta=tol)
        self.assertAlmostEqual(p.pipe_inner_diameter, 0.0273, delta=tol)
        self.assertAlmostEqual(p.pipe_length, 100, delta=tol)
        self.assertAlmostEqual(p.thickness, 0.00303, delta=tol)

        # areas
        self.assertAlmostEqual(p.area_cr_outer, 8.761e-4, delta=tol)
        self.assertAlmostEqual(p.area_cr_inner, 5.628e-4, delta=tol)
        self.assertAlmostEqual(p.area_cr_pipe, 3.078e-4, delta=tol)
        self.assertAlmostEqual(p.area_s_outer, 10.4929, delta=tol)
        self.assertAlmostEqual(p.area_s_inner, 8.5851, delta=tol)

        # volumes
        self.assertAlmostEqual(p.total_vol, 0.0876, delta=tol)
        self.assertAlmostEqual(p.fluid_vol, 0.0586, delta=tol)
        self.assertAlmostEqual(p.pipe_wall_vol, 0.0289, delta=tol)

    def test_adopts_user_defined_fluid(self):
        custom_fluid = get_fluid(
            "user_defined",
            viscosity=0.002,
            specific_heat=3200.0,
            density=1050.0,
            conductivity=0.42,
        )
        inputs = self.inputs.copy()
        inputs.pop("fluid_type")
        pipe = Pipe(**inputs, fluid=custom_fluid)

        self.assertIs(pipe.fluid, custom_fluid)

        replacement = get_fluid("water")
        pipe.set_fluid(replacement)
        self.assertIs(pipe.fluid, replacement)

    def test_mdot_to_re(self):
        p = Pipe(**self.inputs)
        tol = 0.1
        self.assertAlmostEqual(p.mdot_to_re(0.1, 20), 4649.9, delta=tol)

    def test_calc_friction_factor(self):
        p = Pipe(**self.inputs)
        tol = 1e-4

        # laminar tests
        re = 100
        self.assertEqual(p.friction_factor(re), 64.0 / re)

        re = 1000
        self.assertEqual(p.friction_factor(re), 64.0 / re)

        re = 1400
        self.assertEqual(p.friction_factor(re), 64.0 / re)

        # transitional tests
        re = 2000
        self.assertAlmostEqual(p.friction_factor(re), 0.03213, delta=tol)

        re = 3000
        self.assertAlmostEqual(p.friction_factor(re), 0.03344, delta=tol)

        re = 4000
        self.assertAlmostEqual(p.friction_factor(re), 0.04127, delta=tol)

        # turbulent tests
        re = 5000
        self.assertEqual(p.friction_factor(re), (0.79 * log(re) - 1.64) ** (-2.0))

        re = 15000
        self.assertEqual(p.friction_factor(re), (0.79 * log(re) - 1.64) ** (-2.0))

        re = 25000
        self.assertEqual(p.friction_factor(re), (0.79 * log(re) - 1.64) ** (-2.0))

    def test_laminar_nusselt(self):
        p = Pipe(**self.inputs)
        tol = 0.01
        self.assertAlmostEqual(p.laminar_nusselt(), 4.01, delta=tol)

    def test_turbulent_nusselt(self):
        p = Pipe(**self.inputs)
        tol = 0.01
        self.assertAlmostEqual(p.turbulent_nusselt(3000, 20), 18.39, delta=tol)
        self.assertAlmostEqual(p.turbulent_nusselt(10000, 20), 79.50, delta=tol)

    def test_calc_conduction_resistance(self):
        p = Pipe(**self.inputs)
        tolerance = 0.00001
        self.assertAlmostEqual(p.calc_cond_resist(), 0.0798443, delta=tolerance)

    def test_calc_convection_resistance(self):
        p = Pipe(**self.inputs)
        temp = 20
        tol = 0.00001
        self.assertAlmostEqual(p.calc_conv_resist(0, temp), 0.13266, delta=tol)
        self.assertAlmostEqual(p.calc_conv_resist(0.07, temp), 0.020784, delta=tol)
        self.assertAlmostEqual(p.calc_conv_resist(2, temp), 0.00094, delta=tol)

    def test_calc_resist(self):
        pipe = Pipe(**self.inputs)
        tol = 0.00001
        self.assertAlmostEqual(pipe.calc_fluid_pipe_resist(0.5, 20), 0.082985, delta=tol)

    def test_pressure_loss(self):
        pipe = Pipe(**self.inputs)
        tol = 1
        self.assertAlmostEqual(pipe.pressure_loss(0.5, 20), 33533, delta=tol)
