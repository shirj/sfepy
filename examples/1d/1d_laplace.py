r"""
Laplace equation in 1D.
"""
filename_mesh = '1d.mesh'

# Or...
import numpy as nm

# Coordinates of mesh vertices. Elements are consecutive values.
#mesh = nm.linspace(0, 1, 11)

materials = {
    'coef' : ({'val' : 1.0},),
}

regions = {
    'Omega' : 'all', # or 'cells of group 6'
    'Gamma_Left' : ('vertices in (x < 0.00001)', 'facet'),
    'Gamma_Right' : ('vertices in (x > 3.99999)', 'facet'),
}

fields = {
    'temperature' : ('real', 1, 'Omega', 1),
}

variables = {
    't' : ('unknown field', 'temperature', 0),
    's' : ('test field',    'temperature', 't'),
}

ebcs = {
    't1' : ('Gamma_Left', {'t.0' : 2.0}),
    't2' : ('Gamma_Right', {'t.0' : -2.0}),
}

integrals = {
    'i' : 2,
}

equations = {
    'Temperature' : """dw_laplace.i.Omega( coef.val, s, t ) = 0"""
}

solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton',
                {'i_max'      : 1,
                 'eps_a'      : 1e-10,
    }),
}

options = {
    'nls' : 'newton',
    'ls' : 'ls',
}
