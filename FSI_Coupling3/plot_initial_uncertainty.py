import numpy as np
import matplotlib.pyplot as plt

nominal_hot_watergap = 0.002  # Nominal hot water gap in meters

def plot_assembly_deformation(assemblies, filename='deformation0.png'):
    plt.rcParams['axes.linewidth'] = 2.0
    fig, ax = plt.subplots(1, 1, figsize=(15, 8))
    plt.subplots_adjust(left=0.05, right=0.95, top=0.93, bottom=0.1)
    y = np.arange(1, 11)  # Grid positions (1 to 10)
    colors = {
        1: 'darkblue',
        2: 'darkblue',
        3: 'darkblue',
        14: 'darkblue',
        15: 'darkblue',
        6: 'steelblue',
        7: 'steelblue',
        8: 'steelblue',
        9: 'steelblue',
        10: 'steelblue',
        11: 'darkred',
        4: 'darkred',
        5: 'darkred',
        12: 'darkred',
        13: 'darkblue',
    }
    for i, x in enumerate(assemblies):
        x_plot = x[::-1]
        plt.xlim(nominal_hot_watergap - nominal_hot_watergap * 1.0, 15 * nominal_hot_watergap + nominal_hot_watergap * 1.0)
        plt.ylim(0., 11.)
        xtick_positions = [nominal_hot_watergap + nominal_hot_watergap * k for k in range(15)]
        plt.xticks(xtick_positions, [str(k + 1) for k in range(15)])
        plt.yticks(y, [str(k) for k in range(1, 11)])  # Grid labels
        plt.tick_params(labelleft=True, left=True, labelsize=15)
        plt.tick_params(axis='x', labelsize=15, labelrotation=0)
        plt.ylabel('Grids', fontsize=25)
        plt.xlabel('Assemblies', fontsize=25)
        plt.grid(True, which='major', linestyle='-')
        plt.plot(x_plot, y, marker='o', linestyle='-', color=colors.get(i+1))
    plt.savefig(filename)
    plt.close('all')

def GramSchmidt():
    f_0 = np.array([0.0, 1.09090909e-04, 3.63636364e-04, 6.72727273e-04, 9.27272727e-04,
                    1.00000000e-03, 8.54545455e-04, 5.63636364e-04, 2.54545455e-04, 0.0])
    f_1 = np.array([0.0, 5.4e-04, 1.0e-03, 9.6e-04, 3.8e-04, -3.8e-04, -8.6e-04,
                    -9.4e-04, -6.2e-04, 0.0])
    f_2 = np.array([0.0, 7.6e-04, 1.0e-03, 2.0e-04, -9.4e-04, -8.6e-04, 1.6e-04,
                    9.2e-04, 9.0e-04, 0.0])
    M = np.zeros((3, 10))
    M[0] = f_0 / np.linalg.norm(f_0)
    v = f_1 - np.dot(f_0, f_1) / np.dot(f_0, f_0) * f_0
    M[1] = v / np.linalg.norm(v)
    v = f_2 - np.dot(f_0, f_2) / np.dot(f_0, f_0) * f_0 - np.dot(M[1], f_2) / np.dot(M[1], M[1]) * M[1]
    M[2] = v / np.linalg.norm(v)
    return M[0], M[1], M[2]

# Example: plot one sample of initial uncertainty for CSW0
C= [-1.80640018e-04,  3.41280644e-04, -2.29175954e-03, -9.60219314e-05,
            -8.49123937e-05,  4.15388725e-03,  3.26323131e-03,  3.59113552e-04,
            -1.74406149e-03, -1.54797464e-03,  1.04884639e-05,  6.70227483e-05,
            -7.77740556e-05,  9.72333656e-04,  8.50323877e-04]
S = [8.00304835e-04,  1.07567214e-04, -1.44638224e-03, -1.68008800e-05,
        -1.70367470e-05, -2.67415013e-03, -1.01068999e-04,  2.60543181e-04,
        -1.35827353e-03,  3.66635031e-03,  6.36241473e-05,  3.47234554e-06,
        1.15619573e-03, -1.39904720e-03, -7.93415783e-05]
W = [9.41636894e-05,  3.19043207e-04, -3.31703156e-04,  2.23100689e-05,
            2.42452488e-05,  7.73366890e-04,  3.42280097e-04, -1.34272389e-03,
            -7.74612902e-04,  1.24295279e-04,  4.27688330e-05, -9.52701364e-06,
            -9.23258449e-04, -4.28579229e-04,  1.12371381e-04]

Sigma_mesurmentsC = 0.0003
Sigma_mesurmentsS = 0.0003
sigma_mesurmentsW = 0.0003
cov_diagC = Sigma_mesurmentsC ** 2
cov_diagS = Sigma_mesurmentsS ** 2
cov_diagW = sigma_mesurmentsW ** 2
SigmaC = np.diag([cov_diagC] * 15)
SigmaS = np.diag([cov_diagS] * 15)
SigmaW = np.diag([cov_diagW] * 15)


# Compute the standard deviation for each coefficient (diagonal of covariance)
std_C = np.sqrt(np.diag(SigmaC))
std_S = np.sqrt(np.diag(SigmaS))
std_W = np.sqrt(np.diag(SigmaW))

c, s, w = GramSchmidt()
assemblies = []
assemblies_std = []
for i in range(15):
    x = np.zeros(10)
    std_x = np.zeros(10)
    for j in range(10):
        x[j] = nominal_hot_watergap + (i) * nominal_hot_watergap
        deformation_C = C[i] * c[j]
        deformation_S = S[i] * s[j]
        deformation_W = W[i] * w[j]
        x[j] += (deformation_C + deformation_S + deformation_W)
        # Uncertainty propagation (sum of variances)
        var_deformation = (std_C[i] * c[j])**2 + (std_S[i] * s[j])**2 + (std_W[i] * w[j])**2
        std_x[j] = np.sqrt(var_deformation)
    assemblies.append(x)
    assemblies_std.append(std_x)

# Plot all assemblies in the same figure with their uncertainty
plt.rcParams['axes.linewidth'] = 2.0
fig, ax = plt.subplots(1, 1, figsize=(15, 8))
plt.subplots_adjust(left=0.05, right=0.95, top=0.93, bottom=0.1)
y = np.arange(1, 11)
colors = [
    'darkblue', 'darkblue', 'darkblue', 'darkred', 'darkred',
    'steelblue', 'steelblue', 'steelblue', 'steelblue', 'steelblue',
    'darkred', 'darkred', 'darkblue', 'darkblue', 'darkblue'
]
for i in range(15):
    x = assemblies[i]
    std_x = assemblies_std[i]
    x_plot = x[::-1]
    std_plot = std_x[::-1]
    ax.plot(x_plot, y, marker='o', linestyle='-', color=colors[i], label=f'A{i+1}' if i < 10 else None)
    ax.fill_betweenx(
        y,
        x_plot - 2 * std_plot,
        x_plot + 2 * std_plot,
        color=colors[i],
        alpha=0.15
    )
xtick_positions = [nominal_hot_watergap + nominal_hot_watergap * k for k in range(15)]
ax.set_xlim(nominal_hot_watergap - nominal_hot_watergap * 1.0, 15 * nominal_hot_watergap + nominal_hot_watergap * 1.0)
ax.set_ylim(0., 11.)
ax.set_xticks(xtick_positions)
ax.set_xticklabels([str(k + 1) for k in range(15)])
ax.set_yticks(y)
ax.set_yticklabels([str(k) for k in range(1, 11)])
ax.tick_params(labelleft=True, left=True, labelsize=15)
ax.tick_params(axis='x', labelsize=15, labelrotation=0)
ax.set_ylabel('Grids', fontsize=25)
ax.set_xlabel('Assemblies', fontsize=25)
ax.grid(True, which='major', linestyle='-')
ax.set_title("Initial Deformation Uncertainty for All Assemblies due to Measurement Uncertainty", fontsize=20)
plt.savefig('initial_uncertainty_all_assemblies.png')
plt.close('all')