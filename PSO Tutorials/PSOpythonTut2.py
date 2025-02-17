# Import modules
import matplotlib.pyplot as plt
from IPython.display import Image

# Import PySwarms
import pyswarms as ps
from pyswarms.utils.functions import single_obj as fx
from pyswarms.utils.plotters import (plot_cost_history, plot_contour, plot_surface)

# import modules
import numpy as np

# create a parameterized version of the classic Rosenbrock unconstrained optimzation function
def rosenbrock_with_args(x, a, b, c=0):
    f = (a - x[:, 0]) ** 2 + b * (x[:, 1] - x[:, 0] ** 2) ** 2 + c
    return (-f)

from pyswarms.single.global_best import GlobalBestPSO

# instatiate the optimizer
x_max = 10 * np.ones(2)
x_min = -1 * x_max
bounds = (x_min, x_max)
options = {'c1': 0.5, 'c2': 0.3, 'w': 0.9}
optimizer = GlobalBestPSO(n_particles=10, dimensions=2, options=options, bounds=bounds)

# # now run the optimization, pass a=1 and b=100 as a tuple assigned to args

# cost, pos = optimizer.optimize(rosenbrock_with_args, 1000, a=1, b=100, c=0)


# kwargs={"a": 1.0, "b": 100.0, 'c':0}
# cost, pos = optimizer.optimize(rosenbrock_with_args, 1000, **kwargs)

cost, pos = optimizer.optimize(rosenbrock_with_args, 1000, a=1, b=100)

plot_cost_history(cost_history=optimizer.cost_history)
plt.show()