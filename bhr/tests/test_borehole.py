import unittest

from scp import get_fluid

from bhr.borehole import Borehole


class TestBorehole(unittest.TestCase):
    def test_adopts_user_defined_fluid(self):
        custom_fluid = get_fluid(
            "user_defined",
            name="BoreholeFluid",
            viscosity=0.002,
            specific_heat=3200.0,
            density=1050.0,
            conductivity=0.42,
        )
        bh = Borehole()
        bh.init_single_u_borehole(
            borehole_diameter=0.14,
            pipe_outer_diameter=0.042,
            pipe_dimension_ratio=11,
            length=100,
            shank_space=0.01,
            pipe_conductivity=0.4,
            grout_conductivity=1.2,
            soil_conductivity=2.5,
            fluid=custom_fluid,
        )

        self.assertIsNotNone(bh._bh)
        self.assertIs(bh._bh.fluid, custom_fluid)
        self.assertGreater(bh.calc_bh_resist(temperature=20, mass_flow_rate=0.5), 0)

        replacement = get_fluid("water")
        bh.set_fluid(replacement)
        self.assertIs(bh._bh.fluid, replacement)

    def test_set_fluid_requires_initialized_borehole(self):
        with self.assertRaisesRegex(RuntimeError, "Initialize the borehole"):
            Borehole().set_fluid(get_fluid("water"))

    def test_init_single_u_uhf(self):
        bh = Borehole()
        bh.init_single_u_borehole(
            borehole_diameter=0.14,
            pipe_outer_diameter=0.042,
            pipe_dimension_ratio=11,
            length=100,
            shank_space=0.01,
            pipe_conductivity=0.4,
            grout_conductivity=1.2,
            soil_conductivity=2.5,
            fluid_type="PROPYLENEGLYCOL",
            fluid_concentration=0.2,
        )

        self.assertAlmostEqual(bh.calc_bh_resist(temperature=20, mass_flow_rate=0.5), 0.20425, delta=1e-4)
        self.assertAlmostEqual(bh.calc_pipe_cond_resist(), 0.07984, delta=1e-4)
        self.assertAlmostEqual(bh.calc_fluid_resist(temperature=20, mass_flow_rate=0.5), 0.006449, delta=1e-4)
        self.assertAlmostEqual(bh.calc_fluid_pipe_resist(temperature=20, mass_flow_rate=0.5), 0.086336, delta=1e-4)

    def test_init_single_u_ubwt(self):
        bh = Borehole()
        bh.init_single_u_borehole(
            borehole_diameter=0.14,
            pipe_outer_diameter=0.042,
            pipe_dimension_ratio=11,
            length=100,
            shank_space=0.01,
            pipe_conductivity=0.4,
            grout_conductivity=1.2,
            soil_conductivity=2.5,
            fluid_type="PROPYLENEGLYCOL",
            fluid_concentration=0.2,
            boundary_condition="uniform_borehole_wall_temp",
        )

        self.assertAlmostEqual(bh.calc_bh_resist(temperature=20, mass_flow_rate=0.5), 0.20414, delta=1e-4)
        self.assertAlmostEqual(bh.calc_pipe_cond_resist(), 0.07984, delta=1e-4)
        self.assertAlmostEqual(bh.calc_fluid_resist(temperature=20, mass_flow_rate=0.5), 0.006449, delta=1e-4)
        self.assertAlmostEqual(bh.calc_fluid_pipe_resist(temperature=20, mass_flow_rate=0.5), 0.086336, delta=1e-4)

    def test_init_single_u_from_dict(self):
        inputs = {
            "fluid_type": "PROPYLENEGLYCOL",
            "fluid_concentration": 0.2,
            "boundary_condition": "uniform_heat_flux",
            "borehole_type": "single_u_tube",
            "single_u_tube": {
                "pipe_outer_diameter": 0.042,
                "pipe_dimension_ratio": 11,
                "pipe_conductivity": 0.4,
                "shank_space": 0.01,
            },
            "grout_conductivity": 1.2,
            "soil_conductivity": 2.5,
            "length": 100,
            "borehole_diameter": 0.14,
        }

        bh = Borehole()
        bh.init_from_dict(inputs)

        # only pass flow rate, so pipe resistance should be computed in the process of this call

        self.assertAlmostEqual(bh.calc_bh_resist(temperature=20, mass_flow_rate=0.5), 0.20425, delta=1e-4)
        self.assertAlmostEqual(bh.calc_pipe_cond_resist(), 0.07984, delta=1e-4)
        self.assertAlmostEqual(bh.calc_fluid_resist(temperature=20, mass_flow_rate=0.5), 0.006449, delta=1e-4)
        self.assertAlmostEqual(bh.calc_fluid_pipe_resist(temperature=20, mass_flow_rate=0.5), 0.086336, delta=1e-4)

        # try bh wall temp boundary condition
        inputs.update({"boundary_condition": "uniform_borehole_wall_temp"})

        bh_2 = Borehole()
        bh_2.init_from_dict(inputs)

        # only pass flow rate, so pipe resistance should be computed in the process of this call

        self.assertAlmostEqual(bh_2.calc_bh_resist(temperature=20, mass_flow_rate=0.5), 0.20414, delta=1e-4)
        self.assertAlmostEqual(bh.calc_pipe_cond_resist(), 0.07984, delta=1e-4)
        self.assertAlmostEqual(bh.calc_fluid_resist(temperature=20, mass_flow_rate=0.5), 0.006449, delta=1e-4)
        self.assertAlmostEqual(bh.calc_fluid_pipe_resist(temperature=20, mass_flow_rate=0.5), 0.086336, delta=1e-4)

    def test_init_double_u_uhf(self):
        bh = Borehole()
        bh.init_double_u_borehole(
            borehole_diameter=0.115,
            pipe_outer_diameter=0.032,
            pipe_dimension_ratio=18.9,
            length=200,
            shank_space=0.02263,
            pipe_conductivity=0.389,
            pipe_inlet_arrangement="ADJACENT",
            grout_conductivity=1.5,
            soil_conductivity=3,
            fluid_type="WATER",
            fluid_concentration=0,
        )

        # only pass flow rate, so pipe resistance should be computed in the process of this call
        self.assertAlmostEqual(bh.calc_bh_resist(temperature=20, mass_flow_rate=0.4154), 0.1090, delta=1e-4)
        self.assertAlmostEqual(bh.calc_pipe_cond_resist(), 0.045761, delta=1e-4)
        self.assertAlmostEqual(bh.calc_fluid_resist(temperature=20, mass_flow_rate=0.5), 0.006077, delta=1e-4)
        self.assertAlmostEqual(bh.calc_fluid_pipe_resist(temperature=20, mass_flow_rate=0.5), 0.051838, delta=1e-4)

    def test_init_double_u_ubwt(self):
        bh = Borehole()
        bh.init_double_u_borehole(
            borehole_diameter=0.115,
            pipe_outer_diameter=0.032,
            pipe_dimension_ratio=18.9,
            length=200,
            shank_space=0.02263,
            pipe_conductivity=0.389,
            pipe_inlet_arrangement="ADJACENT",
            grout_conductivity=1.5,
            soil_conductivity=3,
            fluid_type="WATER",
            fluid_concentration=0,
            boundary_condition="uniform_borehole_wall_temp",
        )

        # only pass flow rate, so pipe resistance should be computed in the process of this call
        self.assertAlmostEqual(bh.calc_bh_resist(temperature=20, mass_flow_rate=0.4154), 0.1065, delta=1e-4)
        self.assertAlmostEqual(bh.calc_pipe_cond_resist(), 0.045761, delta=1e-4)
        self.assertAlmostEqual(bh.calc_fluid_resist(temperature=20, mass_flow_rate=0.5), 0.006077, delta=1e-4)
        self.assertAlmostEqual(bh.calc_fluid_pipe_resist(temperature=20, mass_flow_rate=0.5), 0.051838, delta=1e-4)

    def test_double_u_component_resistances_use_per_circuit_flow(self):
        bh = Borehole()
        bh.init_double_u_borehole(
            borehole_diameter=0.115,
            pipe_outer_diameter=0.032,
            pipe_dimension_ratio=18.9,
            length=200,
            shank_space=0.02263,
            pipe_conductivity=0.389,
            pipe_inlet_arrangement="ADJACENT",
            grout_conductivity=1.5,
            soil_conductivity=3,
            fluid_type="WATER",
        )

        total_mass_flow_rate = 0.5
        per_circuit_mass_flow_rate = total_mass_flow_rate / 2
        expected_fluid_resistance = bh._bh.calc_conv_resist(per_circuit_mass_flow_rate, 20)
        expected_fluid_pipe_resistance = bh._bh.calc_fluid_pipe_resist(per_circuit_mass_flow_rate, 20)

        self.assertAlmostEqual(bh.calc_fluid_resist(total_mass_flow_rate, 20), expected_fluid_resistance)
        self.assertAlmostEqual(bh.calc_fluid_pipe_resist(total_mass_flow_rate, 20), expected_fluid_pipe_resistance)

    def test_init_double_u_from_dict(self):
        inputs = {
            "fluid_type": "WATER",
            "fluid_concentration": 0,
            "boundary_condition": "uniform_heat_flux",
            "borehole_type": "double_u_tube",
            "double_u_tube": {
                "pipe_outer_diameter": 0.032,
                "pipe_dimension_ratio": 18.9,
                "pipe_conductivity": 0.389,
                "shank_space": 0.02263,
                "pipe_inlet_arrangement": "ADJACENT",  # or DIAGONAL
            },
            "grout_conductivity": 1.5,
            "soil_conductivity": 3,
            "length": 200,
            "borehole_diameter": 0.115,
        }

        bh = Borehole()
        bh.init_from_dict(inputs)

        # only pass flow rate, so pipe resistance should be computed in the process of this call
        self.assertAlmostEqual(bh.calc_bh_resist(temperature=20, mass_flow_rate=0.4154), 0.1090, delta=1e-4)
        self.assertAlmostEqual(bh.calc_pipe_cond_resist(), 0.045761, delta=1e-4)
        self.assertAlmostEqual(bh.calc_fluid_resist(temperature=20, mass_flow_rate=0.5), 0.006077, delta=1e-4)
        self.assertAlmostEqual(bh.calc_fluid_pipe_resist(temperature=20, mass_flow_rate=0.5), 0.051838, delta=1e-4)

    def test_init_coaxial_uhf(self):
        bh = Borehole()
        bh.init_coaxial_borehole(
            borehole_diameter=0.115,
            outer_pipe_outer_diameter=0.064,
            outer_pipe_dimension_ratio=11,
            outer_pipe_conductivity=0.389,
            inner_pipe_outer_diameter=0.032,
            inner_pipe_dimension_ratio=11,
            inner_pipe_conductivity=0.389,
            length=200,
            grout_conductivity=1.5,
            soil_conductivity=3.0,
            fluid_type="WATER",
            fluid_concentration=0,
        )

        # only pass flow rate, so pipe resistance should be computed in the process of this call
        self.assertAlmostEqual(bh.calc_bh_resist(mass_flow_rate=0.5, temperature=20), 0.18128, delta=1e-4)
        self.assertAlmostEqual(bh.calc_pipe_cond_resist(), 0.082102, delta=1e-4)
        self.assertAlmostEqual(bh.calc_fluid_resist(temperature=20, mass_flow_rate=0.5), 0.008727, delta=1e-4)
        self.assertAlmostEqual(bh.calc_fluid_pipe_resist(temperature=20, mass_flow_rate=0.5), 0.085416, delta=1e-4)

    def test_init_coaxial_ubwt(self):
        bh = Borehole()
        bh.init_coaxial_borehole(
            borehole_diameter=0.115,
            outer_pipe_outer_diameter=0.064,
            outer_pipe_dimension_ratio=11,
            outer_pipe_conductivity=0.389,
            inner_pipe_outer_diameter=0.032,
            inner_pipe_dimension_ratio=11,
            inner_pipe_conductivity=0.389,
            length=200,
            grout_conductivity=1.5,
            soil_conductivity=3.0,
            fluid_type="WATER",
            fluid_concentration=0,
            boundary_condition="uniform_borehole_wall_temp",
        )

        # only pass flow rate, so pipe resistance should be computed in the process of this call
        self.assertAlmostEqual(bh.calc_bh_resist(mass_flow_rate=0.5, temperature=20), 0.18454, delta=1e-4)
        self.assertAlmostEqual(bh.calc_pipe_cond_resist(), 0.082102, delta=1e-4)
        self.assertAlmostEqual(bh.calc_fluid_resist(temperature=20, mass_flow_rate=0.5), 0.00872, delta=1e-4)
        self.assertAlmostEqual(bh.calc_fluid_pipe_resist(temperature=20, mass_flow_rate=0.5), 0.085416, delta=1e-4)

    def test_init_coaxial_from_dict(self):
        inputs = {
            "fluid_type": "WATER",
            "fluid_concentration": 0,
            "boundary_condition": "uniform_heat_flux",
            "borehole_type": "coaxial",
            "coaxial": {
                "outer_pipe_outer_diameter": 0.064,
                "outer_pipe_dimension_ratio": 11,
                "outer_pipe_conductivity": 0.389,
                "inner_pipe_outer_diameter": 0.032,
                "inner_pipe_dimension_ratio": 11,
                "inner_pipe_conductivity": 0.389,
            },
            "grout_conductivity": 1.5,
            "soil_conductivity": 3,
            "length": 200,
            "borehole_diameter": 0.115,
        }

        bh = Borehole()
        bh.init_from_dict(inputs)

        # only pass flow rate, so pipe resistance should be computed in the process of this call
        self.assertAlmostEqual(bh.calc_bh_resist(mass_flow_rate=0.5, temperature=20), 0.18128, delta=1e-4)
        self.assertAlmostEqual(bh.calc_pipe_cond_resist(), 0.082102, delta=1e-4)
        self.assertAlmostEqual(bh.calc_fluid_resist(temperature=20, mass_flow_rate=0.5), 0.00872, delta=1e-4)
        self.assertAlmostEqual(bh.calc_fluid_pipe_resist(temperature=20, mass_flow_rate=0.5), 0.085416, delta=1e-4)
