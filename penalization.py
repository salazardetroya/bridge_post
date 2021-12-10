from firedrake.constant import Constant

def ramp(rho, ramp_p=20.0, val_1=1.0, val_0=0.0):
    """RAMP penalization

    Args:
        rho (float): Volume fraction function
        ramp_p (float, optional): Penalization parameter. Defaults to 30.0. Pick value from > 0.0
        val_1 (float, optional): Function value for rho=1. Defaults to 1.0.
        val_0 (float, optional): Function value for rho=0. Defaults to 0.0.

    Returns:
        float: Penalized material property
    """
    if isinstance(ramp_p, float):
        assert ramp_p > 0, "ramp_p has to be positive"
    elif isinstance(ramp_p, Constant):
        assert ramp_p.dat.data[0] > 0, "ramp_p has to be positive"
    if val_1 > val_0:
        return (rho) / (Constant(1.0) + ramp_p * (Constant(1.0) - rho)) * (
            Constant(val_1 - val_0)
        ) + Constant(val_0)
    else:
        return (Constant(1.0) - rho) / (Constant(1.0) + ramp_p * rho) * (
            Constant(val_0 - val_1)
        ) + Constant(val_1)


def inv_ramp(rho, ramp_p=20.0, val_1=1.0, val_0=0.0):
    """Inverse RAMP penalization
        Intermediate values become greater

    Args:
        rho (float): Volume fraction function
        ramp_p (float, optional): Penalization parameter. Defaults to 30.0. Pick value from > 0.0
        val_1 (float, optional): Function value for rho=1. Defaults to 1.0.
        val_0 (float, optional): Function value for rho=0. Defaults to 0.0.

    Returns:
        float: Penalized material property
    """
    if isinstance(ramp_p, float):
        assert ramp_p > 0, "ramp_p has to be positive"
    elif isinstance(ramp_p, Constant):
        assert ramp_p.dat.data[0] > 0, "ramp_p has to be positive"
    if val_0 > val_1:
        return (rho) / (Constant(1.0) + ramp_p * (Constant(1.0) - rho)) * (
            Constant(val_1 - val_0)
        ) + Constant(val_0)
    else:
        return (Constant(1.0) - rho) / (Constant(1.0) + ramp_p * rho) * (
            Constant(val_0 - val_1)
        ) + Constant(val_1)


def simp(rho, simp_p=3.0, val_1=1.0, val_0=0.0):
    """SIMP penalization

    Args:
        rho (float): Volume fraction
        simp_p (float, optional): Penalization parameter. Defaults to 3.0.
        val_1 (float, optional): Function value for rho=1. Defaults to 1.0.
        val_0 (float, optional): Function value for rho=0. Defaults to 0.0.

    Returns:
        float: Penalized material properties
    """
    if val_1 > val_0:
        return (rho) ** simp_p * (val_1 - val_0) + val_0
    else:
        return (Constant(1.0) - rho) ** Constant(simp_p) * (val_1 - val_0) + val_0


def volume(rho, val_1=1.0, val_0=0.0):
    if val_1 > val_0:
        return rho * (val_1 - val_0) + val_0
    else:
        return (Constant(1.0) - rho) * (val_0 - val_1) + val_1

