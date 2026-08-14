from unittest import TestCase
from unittest.mock import patch

from bhr.double_u_borehole import DoubleUTube


class TestDoubleUBorehole(TestCase):
    def setUp(self):
        # from example in Claesson & Javed 2019 paper
        self.inputs = {
            "borehole_diameter": 0.115,  # confirmed
            "pipe_outer_diameter": 0.032,  # confirmed
            "pipe_dimension_ratio": 20.164,  # tuned to set Rp = 0.05
            "length": 200,  # confirmed
            "shank_space": 0.02263,  # confirmed
            "pipe_conductivity": 0.389,  # reasonable value, but ultimately effects come out in Rp
            "pipe_inlet_arrangement": "DIAGONAL",  # correct
            "grout_conductivity": 1.5,  # confirmed
            "soil_conductivity": 3,  # confirmed
            "fluid_type": "WATER",  # effects come out in Rp
        }

        self.v_dot_bh = 1.5 / 3600  # m3/s
        self.v_dot_per_u = self.v_dot_bh / 2
        self.rho_test = 997
        self.m_dot_bh = self.v_dot_bh * self.rho_test
        self.m_dot_per_u = self.v_dot_per_u * self.rho_test

    def test_init(self):
        bh = DoubleUTube(**self.inputs)
        self.assertEqual(bh.bh_length, 200)

    def test_shank_space_assert(self):
        # checks if code is raising error when shank space is too large
        d = self.inputs.copy()
        with self.assertRaises(AssertionError):
            d.update({"shank_space": 0.115})
            DoubleUTube(**d)

        # checks if code is raising error when shank space is too small
        d = self.inputs.copy()
        with self.assertRaises(AssertionError):
            d.update({"shank_space": 0.002})
            DoubleUTube(**d)

    def test_update_beta_b1(self):
        bh = DoubleUTube(**self.inputs)

        tolerance = 1e-3

        self.assertAlmostEqual(bh.update_b1(m_dot_per_u_tube=self.m_dot_per_u, temperature=20), 0.359, delta=tolerance)
        self.assertAlmostEqual(bh.pipe_resist, 0.05, delta=tolerance)

    def test_calc_bh_resist_local(self):
        bh = DoubleUTube(**self.inputs)
        tolerance = 1e-3
        self.assertAlmostEqual(
            bh.calc_bh_resist_local(m_dot_per_u_tube=self.m_dot_per_u, temperature=20), 7.509e-02, delta=tolerance
        )

    def test_calc_internal_resist(self):
        bh = DoubleUTube(**self.inputs)
        tolerance = 1e-3
        self.assertAlmostEqual(
            bh.calc_internal_resist(m_dot_per_u_tube=self.m_dot_per_u, temperature=20), 0.1604, delta=tolerance
        )

    def test_calc_resistances_diagonal(self):
        d = self.inputs.copy()
        d.update({"pipe_inlet_arrangement": "DIAGONAL"})
        bh = DoubleUTube(**d)
        tolerance = 1e-3

        # values match test case data in Table 1
        self.assertAlmostEqual(
            bh.calc_effective_bh_resistance_uhf(m_dot=self.m_dot_bh, temp=20), 0.1302, delta=tolerance
        )
        self.assertAlmostEqual(
            bh.calc_effective_bh_resistance_ubwt(m_dot=self.m_dot_bh, temp=20), 0.1235, delta=tolerance
        )

    def test_calc_resistances_adjacent(self):
        d = self.inputs.copy()
        d.update({"pipe_inlet_arrangement": "ADJACENT"})
        bh = DoubleUTube(**d)
        tolerance = 1e-3

        # values match test case data in Table 1
        self.assertAlmostEqual(
            bh.calc_effective_bh_resistance_uhf(m_dot=self.m_dot_bh, temp=20), 0.1089, delta=tolerance
        )
        self.assertAlmostEqual(
            bh.calc_effective_bh_resistance_ubwt(m_dot=self.m_dot_bh, temp=20), 0.1062, delta=tolerance
        )

    def test_invalid_pipe_arrangement(self):
        d = self.inputs.copy()
        d.update({"pipe_inlet_arrangement": "UNSUPPORTED"})
        with self.assertRaises(AssertionError):
            DoubleUTube(**d)

    def test_resistance_methods_require_calculated_pipe_resistance(self):
        for method_name in ("calc_bh_resist_local", "calc_internal_resist"):
            with self.subTest(method=method_name):
                bh = DoubleUTube(**self.inputs)
                with (
                    patch.object(bh, "update_b1", return_value=0),
                    self.assertRaisesRegex(ValueError, "Pipe resistance has not been calculated"),
                ):
                    getattr(bh, method_name)(self.m_dot_per_u, 20)

    def test_internal_resistance_rejects_invalid_runtime_arrangement(self):
        bh = DoubleUTube(**self.inputs)
        bh.pipe_inlet_arrangement = None

        with self.assertRaisesRegex(AssertionError, "Invalid pipe inlet arrangement"):
            bh.calc_internal_resist(self.m_dot_per_u, 20)
