import numpy as np

from scipy.linalg import solve_banded
import tqdm

import matplotlib.pyplot


def left_value(E, L=5, N=10 ** 3):
    xs, h = np.linspace(0, L, N, retstep=True)
    
    ab = np.empty((3, N))
    ab[0] = np.ones(N)
    ab[1, -1] = 2 * (1 + h * (2 * L - 2 * E) ** 0.5)
    ab[1, 1:-1] = np.array([-2 * (1 - h ** 2 * (E - m * h)) for m in range(1, N - 1)], dtype=float)
    ab[2] = np.ones(N)
    ab[2, -2] = 2 * (1 - h ** 2 * (E - (N - 2) * h)) - 4
    
    print(ab)
    print()
    
    b = np.zeros(N)
    b[-1] = 1
    
    ys = solve_banded((0, 2), ab, b)
    
    print(ys)
    
    return ab


def check_for_eigenvalue(energy, left_energy_bound=5, lattice_size=10 ** 3):
    energy_lattice, lattice_step = np.linspace(0, left_energy_bound, lattice_size, retstep=True)
    
    # описание трехдиагональной матрицы системы
    banded_matrix = np.empty(shape=(3, lattice_size))
    # побочные диагонали заполнены единицами кроме граничных значений
    banded_matrix[0] = np.ones(lattice_size)
    # "сшивка"
    banded_matrix[1, -1] = 2 * (1 + lattice_step * (2 * left_energy_bound - 2 * energy) ** 0.5)
    # главная диагональ
    banded_matrix[1, 1:-1] = np.array([
        -2 * (1 - lattice_step ** 2 * (energy - index * lattice_size)) for index in range(1, lattice_size - 1)
    ], dtype=float)
    banded_matrix[2] = np.ones(lattice_size)
    # "сшивка"
    banded_matrix[2, -2] = 2 * (1 - lattice_step ** 2 * (energy - (lattice_size - 2) * lattice_step)) - 4
    
    print(banded_matrix)
    print()
    
    
    column = np.zeros(lattice_size)
    column[-1] = 1
    
    # решение трехдиагональной системы
    solution = solve_banded((0, 2), banded_matrix, column)
    
    print(solution)
    
    return banded_matrix  # psi in 0


print(check_for_eigenvalue(2, 5, 10) - left_value(2, 5, 10))

# energy_lattice = np.linspace(1.7, 5, 10)
# bound_condition = np.array([check_for_eigenvalue(e) for e in energy_lattice])
# print(bound_condition)

# fig = matplotlib.pyplot.figure(figsize=(7, 5))

# ax = fig.add_subplot(111)

# ax.plot(energy_lattice, bound_condition)
# ax.grid()

# fig.show()
