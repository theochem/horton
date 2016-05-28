try:
    from horton.units import angstrom
except ImportError:
    print("Warning, using built-in units")
    angstrom = 1.0e-10/0.5291772083e-10