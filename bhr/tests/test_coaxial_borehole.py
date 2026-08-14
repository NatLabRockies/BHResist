from unittest import TestCase

from bhr.coaxial_borehole import Coaxial


class TestCoaxialBorehole(TestCase):
    def setUp(self):
        self.inputs = {
            "borehole_diameter": 0.115,
            "outer_pipe_outer_diameter": 0.064,
            "outer_pipe_dimension_ratio": 11,
            "outer_pipe_conductivity": 0.389,
            "inner_pipe_outer_diameter": 0.032,
            "inner_pipe_dimension_ratio": 11,
            "inner_pipe_conductivity": 0.389,
            "length": 200,
            "grout_conductivity": 1.5,
            "soil_conductivity": 3,
            "fluid_type": "WATER",
            "fluid_concentration": 0,
        }

    def test_init(self):
        coax = Coaxial(**self.inputs)
        self.assertEqual(coax.borehole_diameter, 0.115)

    def test_re_annulus(self):
        coax = Coaxial(**self.inputs)

        # laminar flow
        self.assertAlmostEqual(coax.re_annulus(m_dot=0.1, temp=20), 1506.215, delta=1e-3)
        # transitional flow
        self.assertAlmostEqual(coax.re_annulus(m_dot=0.2, temp=20), 3012.430, delta=1e-3)
        # turbulent flow
        self.assertAlmostEqual(coax.re_annulus(m_dot=0.5, temp=20), 7531.077, delta=1e-3)

    def test_laminar_nusselt_annulus(self):
        coax = Coaxial(**self.inputs)

        tolerance = 1e-3

        self.assertAlmostEqual(coax.laminar_nusselt_annulus()[0], 5.4394, delta=tolerance)
        self.assertAlmostEqual(coax.laminar_nusselt_annulus()[1], 4.5980, delta=tolerance)

    def test_turbulent_nusselt_annulus(self):
        coax = Coaxial(**self.inputs)
        self.assertAlmostEqual(coax.turbulent_nusselt_annulus(re=10000, temp=20)[0], 72.0409, delta=1e-3)
        self.assertAlmostEqual(coax.turbulent_nusselt_annulus(re=20000, temp=20)[1], 125.4305, delta=1e-3)

    def test_convective_resist_annulus(self):
        coax = Coaxial(**self.inputs)

        # tests convective resistance of outside surface of inner pipe
        # laminar flow
        self.assertAlmostEqual(coax.calc_conv_resist_annulus(m_dot=0.1, temp=20)[0], 0.062236, delta=1e-3)
        # transitional flow
        self.assertAlmostEqual(coax.calc_conv_resist_annulus(m_dot=0.2, temp=20)[0], 0.051664, delta=1e-3)
        # turbulent flow
        self.assertAlmostEqual(coax.calc_conv_resist_annulus(m_dot=0.5, temp=20)[0], 0.005412, delta=1e-3)

        # tests convective resistance of inside surface of outer pipe
        # laminar flow
        self.assertAlmostEqual(coax.calc_conv_resist_annulus(m_dot=0.1, temp=20)[1], 0.04499, delta=1e-3)
        # transitional flow
        self.assertAlmostEqual(coax.calc_conv_resist_annulus(m_dot=0.2, temp=20)[1], 0.036135, delta=1e-3)
        # turbulent flow
        self.assertAlmostEqual(coax.calc_conv_resist_annulus(m_dot=0.5, temp=20)[1], 0.003314, delta=1e-3)

    def test_calc_local_bh_resistance(self):
        coax = Coaxial(**self.inputs)

        # laminar flow
        self.assertAlmostEqual(coax.calc_local_bh_resistance(m_dot=0.1, temp=20)[0], 0.3472, delta=1e-3)
        # transitional flow
        self.assertAlmostEqual(coax.calc_local_bh_resistance(m_dot=0.2, temp=20)[0], 0.3210, delta=1e-3)
        # turbulent flow
        self.assertAlmostEqual(coax.calc_local_bh_resistance(m_dot=0.5, temp=20)[0], 0.2381, delta=1e-3)

    def test_calc_fluid_pipe_resist_uses_outer_annulus_path(self):
        coax = Coaxial(**self.inputs)
        _, outer_pipe_conduction = coax.calc_cond_resist()
        _, outer_annulus_convection = coax.calc_conv_resist_annulus(m_dot=0.5, temp=20)

        expected = outer_pipe_conduction + outer_annulus_convection
        self.assertAlmostEqual(coax.calc_fluid_pipe_resist(m_dot=0.5, temp=20), expected)

    def test_calc_effective_bh_resistance_uhf(self):
        coax = Coaxial(**self.inputs)

        self.assertAlmostEqual(coax.calc_effective_bh_resistance_uhf(m_dot=0.02, temp=20), 7.07, delta=1e-3)

    def test_calc_effective_bh_resistance_ubwt(self):
        coax = Coaxial(**self.inputs)

        self.assertAlmostEqual(coax.calc_effective_bh_resistance_ubwt(m_dot=0.02, temp=20), 2.31, delta=1e-3)
