from math import exp, tanh

from bhr.enums import BoundaryCondition


def set_boundary_condition_enum(bh_str: str) -> BoundaryCondition:
    boundary_condition_map = {bc.name: bc for bc in BoundaryCondition}
    if bh_str.upper() in boundary_condition_map:
        return boundary_condition_map[bh_str.upper()]

    raise ValueError(f"Invalid boundary condition: '{bh_str}'")


def inch_to_m(x_inch: float) -> float:
    """
    Convert inches to meters

    :param x_inch: value to convert, inches
    :return: float value, meters
    """
    return x_inch * 0.0254


def smoothing_function(x, x_low_limit, x_high_limit, y_low_limit, y_high_limit) -> float:
    """
    Sigmoid smoothing function

    https://en.wikipedia.org/wiki/Sigmoid_function

    :param x: independent variable
    :param x_low_limit: lower limit on x range
    :param x_high_limit: upper limit on x range
    :param y_low_limit: lower limit on y range
    :param y_high_limit: upper limit on y range
    :return: smoothed value between x_low_limit and x_high_limit. returns y_low_limit and y_high_limit
    below and above the x_low_limit and x_high_limit, respectively.
    """

    # check lower and upper bounds
    if x < x_low_limit:
        return y_low_limit

    if x > x_high_limit:
        return y_high_limit

    # interp x range to value between -5 to 5 for sigmoid function, which is picked to be
    # within approx. 1% of the desired end points of the interpolation range, 0 to 1.
    s_x_max = 5
    s_x_min = -5
    s_x = (x - x_low_limit) / (x_high_limit - x_low_limit) * (s_x_max - s_x_min) + s_x_min

    # get sigmoid function for interpolation between real y values
    s_y = 1 / (1 + exp(-s_x))

    # finally interpolate to the final value
    return s_y * (y_high_limit - y_low_limit) + y_low_limit


def coth(x):
    return 1 / tanh(x)
