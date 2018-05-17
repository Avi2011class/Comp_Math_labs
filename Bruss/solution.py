import numpy as np
import sympy as sp
import tqdm
import matplotlib.pyplot as plt

A = 1
B = 3
u0 = [2, 3]
left_bound = 0
right_bound = 20
N = 150

# Compile the functions
x, y = sp.symbols('x y')
u = sp.Matrix([x, y])

G = sp.Matrix([A + x ** 2 * y - (B + 1) * x,
               B * x - x ** 2 * y])
J = G.jacobian(u)

J_function = sp.lambdify((x, y), J)
G_function = sp.lambdify((x, y), G)


def g(vector):
    return G_function(*vector).flatten()


def j(vector):
    return J_function(*vector)


def newton_method(start_vector, function, jacobian, eps=1e-10, max_iterations=20):
    current_solution = np.copy(start_vector)
    for current_iteration in range(max_iterations):
        dx = np.linalg.solve(jacobian(current_solution), -function(current_solution))
        current_solution += dx
        if np.linalg.norm(dx, ord=np.inf) < eps:
            break
    return current_solution


time_lattice, step_size = np.linspace(left_bound, right_bound, N, retstep=True)

# Решение с помощью многошаговой неявной схемы
solution = np.empty((N, 2))
solution[0] = u0
solution[1] = newton_method(solution[0],
                            lambda vector: vector - solution[0] - step_size * g(vector),
                            lambda vector: np.eye(2) - step_size * j(vector))

for current_step in tqdm.tqdm(range(1, N - 1)):
    solution[current_step + 1] = newton_method(
        solution[current_step],
        lambda vector: 3 * vector - 4 * solution[current_step] + solution[current_step - 1]
                       - 2 * step_size * g(vector),
        lambda vector: 3 * np.eye(2) - 2 * step_size * j(vector)
    )

# Решение методом Эйлера с тем же шагом
e_solution = np.empty((N, 2))
e_solution[0] = u0

for current_step in tqdm.tqdm(range(N - 1)):
    e_solution[current_step + 1] = e_solution[current_step] + step_size * g(e_solution[current_step])

fig = plt.figure(figsize=(9, 10))
ax1, ax2, ax3 = fig.subplots(3)

ax = ax1

ax.plot(time_lattice, solution[:, 0], label='$x(t)$ неявная схема')
ax.plot(time_lattice, solution[:, 1], label='$y(t)$ неявная схема')

ax.plot(time_lattice, e_solution[:, 0], linewidth=1, linestyle='-.', label='$x(t)$ явный метод Эйлера')
ax.plot(time_lattice, e_solution[:, 1], linewidth=1, linestyle='-.', label='$y(t)$ явный метод Эйлера')

ax.set_xticks(np.arange(0, 26, 5))
ax.set_xticks(np.arange(0, 26, 1), minor=True)

ax.set_yticks(np.arange(0, 5, 1))
ax.set_yticks(np.arange(0, 5, 0.2), minor=True)

ax.grid(which='major', alpha=1, c='black')
ax.grid(which='minor', alpha=0.5)
ax.legend(loc='upper right')
ax.set_title('Поведение методов на заданном промежутке')

# ----------------------------------------------------------------------------------
ax = ax2

for N in tqdm.tqdm([100, 200, 500, 10000]):  # при значениях, меньшах 100, неустойчивость все портит
    time_lattice, step_size = np.linspace(left_bound, right_bound, N, retstep=True)

    e_solution = np.empty((N, 2))
    e_solution[0] = u0

    for current_step in tqdm.tqdm(range(N - 1)):
        e_solution[current_step + 1] = e_solution[current_step] + step_size * g(e_solution[current_step])

    ax.plot(time_lattice, e_solution[:, 0], linewidth=1, linestyle='-.',
            label='$x(t)$ явный метод Эйлера, $N={}$'.format(N))
    # ax.plot(time_lattice, e_solution[:, 1], linewidth=1, linestyle='-.', label='$y(t)$ явный метод Эйлера')

    ax.set_xticks(np.arange(0, 26, 5))
    ax.set_xticks(np.arange(0, 26, 1), minor=True)

    ax.set_yticks(np.arange(0, 5, 1))
    ax.set_yticks(np.arange(0, 5, 0.2), minor=True)

    ax.grid(which='major', alpha=1, c='black')
    ax.grid(which='minor', alpha=0.5)

    ax.legend(loc='upper right')
    ax.set_title('Поведение метода Эйлера при различном значении шага')
# ----------------------------2-----------------------------------------------------
ax = ax3
left_bound = 0
right_bound = 100
for N in tqdm.tqdm(list([100, 250, 10000])):
    time_lattice, step_size = np.linspace(left_bound, right_bound, N, retstep=True)

    solution = np.empty((N, 2))
    solution[0] = u0
    solution[1] = newton_method(solution[0],
                            lambda vector: vector - solution[0] - step_size * g(vector),
                            lambda vector: np.eye(2) - step_size * j(vector))

    for current_step in tqdm.tqdm(range(1, N - 1)):
        solution[current_step + 1] = newton_method(
            solution[current_step],
            lambda vector: 3 * vector - 4 * solution[current_step] + solution[current_step - 1]
                           - 2 * step_size * g(vector),
            lambda vector: 3 * np.eye(2) - 2 * step_size * j(vector)
        )
    ax.plot(time_lattice, solution[:, 0], label='$x_{%d}(t)$' % N)
    ax.plot(time_lattice, solution[:, 1], label='$y_{%d}(t)$' % N)


    ax.set_xticks(np.arange(0, 100, 20))
    ax.set_xticks(np.arange(0, 100, 5), minor=True)

    ax.set_yticks(np.arange(0, 7, 1))
    ax.set_yticks(np.arange(0, 7, 0.2), minor=True)

    ax.grid(which='major', alpha=1, c='black')
    ax.grid(which='minor', alpha=0.5)
    ax.legend(loc='upper right')
    ax.set_title('Поведение неявного метода в зависимости от величины шага')


# fig.show()
fig.savefig("plot.png")
