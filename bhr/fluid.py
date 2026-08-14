from scp import get_fluid
from scp.base_fluid import BaseFluid

# Deprecated: legacy fluid keys will be removed in the v1.0.0 release.
_LEGACY_FLUID_KEYS = {
    "ETHYLALCOHOL": "ethyl_alcohol",
    "ETHYLENEGLYCOL": "ethylene_glycol",
    "METHYLALCOHOL": "methyl_alcohol",
    "PROPYLENEGLYCOL": "propylene_glycol",
    "WATER": "water",
}


def resolve_fluid(
    fluid_type: str | None = None,
    fluid_concentration: float = 0.0,
    *,
    fluid: BaseFluid | None = None,
) -> BaseFluid:
    """
    Return a SecondaryCoolantProps fluid for use by BHResist.

    Pass either a built-in fluid key and concentration or an existing
    SecondaryCoolantProps BaseFluid instance. Existing instances support
    user-defined fluids created by scp.get_fluid.

    The compact, uppercase keys accepted by earlier BHResist releases remain
    supported and are translated to the canonical SecondaryCoolantProps keys.

    :param fluid_type: built-in fluid key accepted by scp.get_fluid
    :param fluid_concentration: mixture concentration fraction for a built-in fluid
    :param fluid: existing SecondaryCoolantProps fluid instance
    :return: resolved fluid instance
    """
    if fluid is not None:
        if fluid_type is not None:
            raise ValueError("Specify either fluid_type or fluid, not both.")
        if fluid_concentration != 0.0:
            raise ValueError("fluid_concentration cannot be used with an existing fluid instance.")
        if not isinstance(fluid, BaseFluid):
            raise TypeError("fluid must be an instance of scp.base_fluid.BaseFluid.")
        return fluid

    if fluid_type is None:
        raise ValueError("Specify either fluid_type or fluid.")

    fluid_key = _LEGACY_FLUID_KEYS.get(fluid_type.upper(), fluid_type)
    return get_fluid(fluid_key, concentration=fluid_concentration)


__all__ = ["BaseFluid", "get_fluid", "resolve_fluid"]
