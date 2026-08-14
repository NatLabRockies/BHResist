from unittest import TestCase

from bhr.single_u_borehole import SingleUBorehole


class TestSingleUBorehole(TestCase):
    def setUp(self):
        self.inputs = {
            "borehole_diameter": 0.096,
            "pipe_outer_diameter": 0.032,
            "pipe_dimension_ratio": 18.53,
            # This is not a realistic value,
            # set to make pipe conduction + internal fluid resistance = 0.05 to match tests
            "length": 100,
            "shank_space": 0.016,
            "pipe_conductivity": 0.389,
            "grout_conductivity": 0.6,
            "soil_conductivity": 4.0,
            "fluid_type": "WATER",
        }

    def test_init(self):
        bh = SingleUBorehole(**self.inputs)
        self.assertEqual(bh.length, 100)

    def test_update_beta(self):
        bh = SingleUBorehole(**self.inputs)
        bh.update_beta(m_dot=0.5, temp=20)
        self.assertAlmostEqual(bh.pipe_resist, 0.05, delta=0.00001)

    def test_calc_internal_and_grout_resistance(self):
        bh = SingleUBorehole(**self.inputs)
        flow_rate = 0.5
        temperature = 20
        tolerance = 1e-3
        self.assertAlmostEqual(bh.theta_1, 0.33333, delta=tolerance)
        self.assertAlmostEqual(bh.theta_2, 3.0, delta=tolerance)
        self.assertAlmostEqual(bh.calc_total_internal_bh_resistance(flow_rate, temperature), 0.32365, delta=tolerance)
        self.assertAlmostEqual(bh.calc_grout_resistance(flow_rate, temperature), 0.17701, delta=tolerance)

        self.inputs.update({"soil_conductivity": 1.0, "grout_conductivity": 3.6})
        bh = SingleUBorehole(**self.inputs)
        bh.update_beta(m_dot=0.5, temp=20)
        self.assertAlmostEqual(bh.calc_total_internal_bh_resistance(flow_rate, temperature), 0.17456, delta=tolerance)
        self.assertAlmostEqual(bh.calc_grout_resistance(flow_rate, temperature), 0.03373, delta=tolerance)

    def test_grout_resistance_requires_calculated_pipe_resistance(self):
        bh = SingleUBorehole(**self.inputs)

        with self.assertRaisesRegex(ValueError, "Pipe resistance has not been calculated"):
            bh.calc_grout_resistance(m_dot=0.5, temp=20)

    def test_calc_effective_bh_resistance_uhf(self):
        bh = SingleUBorehole(**self.inputs)
        tolerance = 1e-3
        self.assertAlmostEqual(bh.calc_effective_bh_resistance_uhf(m_dot=0.5, temp=20), 0.20435, delta=tolerance)

    def test_calc_effective_bh_resistance_ubwt(self):
        bh = SingleUBorehole(**self.inputs)
        tolerance = 1e-3
        self.assertAlmostEqual(bh.calc_effective_bh_resistance_ubwt(m_dot=0.5, temp=20), 0.20435, delta=tolerance)

    def test_calc_direct_coupling_resistance_preserves_negative_value(self):
        inputs = self.inputs | {
            "borehole_diameter": 0.14,
            "pipe_dimension_ratio": 11,
            "shank_space": 0.04,
            "pipe_conductivity": 0.4,
            "grout_conductivity": 1.2,
            "soil_conductivity": 2.5,
        }
        bh = SingleUBorehole(**inputs)
        total_internal_resistance = bh.calc_total_internal_bh_resistance(m_dot=0.5, temp=20)
        local_borehole_resistance = bh.calc_local_bh_resistance(m_dot=0.5, temp=20)
        expected = (4 * total_internal_resistance * local_borehole_resistance) / (
            4 * local_borehole_resistance - total_internal_resistance
        )

        self.assertLess(expected, 0)
        actual, _ = bh.calc_direct_coupling_resistance(m_dot=0.5, temp=20)
        self.assertAlmostEqual(actual, expected)
