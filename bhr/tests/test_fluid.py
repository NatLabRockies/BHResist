import unittest

from scp import get_fluid as scp_get_fluid

from bhr.fluid import get_fluid, resolve_fluid


class TestFluid(unittest.TestCase):
    def test_reexports_secondary_coolant_props_factory(self):
        self.assertIs(get_fluid, scp_get_fluid)

    def test_secondary_coolant_props_factory(self):
        fluids = {
            "ethyl_alcohol": 968.9,
            "ethylene_glycol": 1024.1,
            "methyl_alcohol": 966.7,
            "propylene_glycol": 1014.7,
        }
        for fluid_type, expected_density in fluids.items():
            with self.subTest(fluid_type=fluid_type):
                fluid = get_fluid(fluid_type, concentration=0.2)
                self.assertAlmostEqual(fluid.density(20), expected_density, delta=0.1)

        water = get_fluid("water")
        self.assertAlmostEqual(water.density(20), 998.2, delta=0.1)

    def test_resolve_fluid_supports_legacy_keys(self):
        fluid = resolve_fluid("PROPYLENEGLYCOL", 0.2)
        self.assertAlmostEqual(fluid.density(20), 1014.7, delta=0.1)

    def test_resolve_fluid_adopts_user_defined_instance(self):
        custom_fluid = get_fluid(
            "user_defined",
            name="BoreholeFluid",
            viscosity=0.002,
            specific_heat=3200.0,
            density=1050.0,
            conductivity=0.42,
        )

        self.assertIs(resolve_fluid(fluid=custom_fluid), custom_fluid)

    def test_resolve_fluid_rejects_ambiguous_or_invalid_inputs(self):
        custom_fluid = get_fluid(
            "user_defined",
            viscosity=0.002,
            specific_heat=3200.0,
            density=1050.0,
            conductivity=0.42,
        )

        with self.assertRaisesRegex(ValueError, "either fluid_type or fluid"):
            resolve_fluid("water", fluid=custom_fluid)
        with self.assertRaisesRegex(ValueError, "fluid_concentration"):
            resolve_fluid(fluid_concentration=0.2, fluid=custom_fluid)
        with self.assertRaisesRegex(TypeError, "BaseFluid"):
            resolve_fluid(fluid=object())
        with self.assertRaisesRegex(ValueError, "either fluid_type or fluid"):
            resolve_fluid()
