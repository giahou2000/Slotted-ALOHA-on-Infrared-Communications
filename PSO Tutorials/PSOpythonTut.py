# Import modules
import matplotlib.pyplot as plt
from IPython.display import Image
# Import modules
import numpy as np

# Import PySwarms
import pyswarms as ps
from pyswarms.utils.functions import single_obj as fx
from pyswarms.utils.plotters import (plot_cost_history, plot_contour, plot_surface)

# # Set-up hyperparameters
# options = {'c1': 0.5, 'c2': 0.3, 'w':0.9}

# # Call instance of PSO
# optimizer = ps.single.GlobalBestPSO(n_particles=10, dimensions=2, options=options)

# # Perform optimization
# cost, pos = optimizer.optimize(fx.sphere, iters=1000)

# # Set-up hyperparameters
# options = {'c1': 0.5, 'c2': 0.3, 'w':0.9, 'k': 2, 'p': 2}

# # Call instance of PSO
# optimizer = ps.single.LocalBestPSO(n_particles=10, dimensions=2, options=options)

# # Perform optimization
# cost, pos = optimizer.optimize(fx.sphere, iters=1000)

# Create bounds
max_bound = 5.12 * np.ones(2)
min_bound = - max_bound
bounds = (min_bound, max_bound)

# Initialize swarm
options = {'c1': 0.5, 'c2': 0.3, 'w':0.9}

# Call instance of PSO with bounds argument
optimizer = ps.single.GlobalBestPSO(n_particles=10, dimensions=2, options=options, bounds=bounds)

# Perform optimization
cost, pos = optimizer.optimize(fx.rastrigin, iters=1000)

plot_cost_history(cost_history=optimizer.cost_history)
plt.show()