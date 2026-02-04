import numpy as np
import matplotlib.pyplot as plt

# 1. Définir les vecteurs de modes
f_0 = np.array([0.0, 1.09090909e-04, 3.63636364e-04, 6.72727273e-04,
                9.27272727e-04, 1.00000000e-03, 8.54545455e-04,
                5.63636364e-04, 2.54545455e-04, 0.0])
f_1 = np.array([0.0, 5.4e-04, 1.0e-03, 9.6e-04, 3.8e-04, -3.8e-04,
                -8.6e-04, -9.4e-04, -6.2e-04, 0.0])
f_2 = np.array([0.0, 7.6e-04, 1.0e-03, 2.0e-04, -9.4e-04, -8.6e-04,
                1.6e-04, 9.2e-04, 9.0e-04, 0.0])

# 2. Orthonormalisation par Gram-Schmidt
M = np.zeros((3, 10))
M[0] = f_0 / np.linalg.norm(f_0)
v = f_1 - np.dot(f_0, f_1) / np.dot(f_0, f_0) * f_0
M[1] = v / np.linalg.norm(v)
v = f_2 - np.dot(f_0, f_2) / np.dot(f_0, f_0) * f_0 - np.dot(M[1], f_2) / np.dot(M[1], M[1]) * M[1]
M[2] = v / np.linalg.norm(v)

C,S,W = M[0], M[1], M[2]  # Coefficients modaux

# Transposer M pour avoir les coefficients en colonnes
plt.figure(figsize=(10, 6))  
plt.plot(M[0]+M[1] + M[2], label='Mode (C+S+W)', color='blue', linewidth=2)
plt.title('C+S+w')
plt.legend()
plt.grid()
plt.savefig('principal_orthogonal_modes.png')
plt.show()


plt.xlabel('Grid Points')
plt.ylabel('Coefficient Value')
#verify orthonormality
assert np.allclose(np.dot(M, M.T), np.eye(3)), "M is not orthonormal"
print(np.max(M[0]+M[1]+M[2]))  # Affiche les valeurs maximales des modes

#M1,M2,M3 are without unit and when we multipluy by the modal coefficients, we get the displacement in metre
M = M.T  # (10, 3)

# 3. Paramètres de simulation
cycle_coeffs = [
    np.array([0.002, 0.001, 0.0005]),  # Cycle 1
    np.array([0.003, 0.002, 0.001]),   # Cycle 2
    np.array([0.006, 0.005, 0.004])    # Cycle 3
]
cycle_names = ['Cycle 1', 'Cycle 2', 'Cycle 3']
cycle_colors = ['darkred', 'navy', 'deepskyblue']
hist_colors = ['#8B0000', '#000080', '#87CEFA']  # darkred, navy, light blue
labels = ['C', 'S', 'W']
sigma = 0.0003  # Bruit de mesure sur chaque grille (en m)
N_sim = 10000   # Nombre d’échantillons Monte Carlo

font_title = 18
font_label = 16
font_ticks = 14
font_legend = 14

for cycle_idx, (true_coeffs, cycle_name, color, hist_color) in enumerate(zip(cycle_coeffs, cycle_names, cycle_colors, hist_colors)):
    # 4. Simulation Monte Carlo
    C_samples = np.zeros((N_sim, 3))
    for i in range(N_sim):
        U_clean = M @ true_coeffs
        U_noisy = U_clean + np.random.normal(0, sigma, size=U_clean.shape)
        C_hat = np.linalg.inv(M.T @ M) @ M.T @ U_noisy
        C_samples[i] = C_hat

    # 5. Moyenne et écart-type estimés des coefficients
    mean_C = np.mean(C_samples, axis=0)
    std_C = np.std(C_samples, axis=0)

    print(f"{cycle_name} - Estimated mean of (C, S, W):", mean_C)
    print(f"{cycle_name} - Estimated std of (C, S, W):", std_C)

    # 6. Affichage des histogrammes pour chaque modal coefficient dans une seule figure
    fig, axs = plt.subplots(1, 3, figsize=(18, 5))
    for j in range(3):
        axs[j].hist(C_samples[:, j], bins=50, density=True, alpha=0.8, color=hist_color)
        axs[j].axvline(mean_C[j], color=color, linestyle='--', linewidth=2, label=f'Mean: {mean_C[j]:.2e} m')
        axs[j].set_xlabel(f'{labels[j]} value (m)', fontsize=font_label)
        axs[j].set_ylabel('Density', fontsize=font_label)
        axs[j].legend(fontsize=font_legend)
        axs[j].grid(True)
        axs[j].tick_params(axis='both', labelsize=font_ticks)
    fig.suptitle(
        f'Distributions of modal coefficient with respect to the uncertainty of measurements (σ = {sigma} m)\n{cycle_name}',
        fontsize=font_title+2
    )
    plt.tight_layout(rect=[0, 0.03, 1, 0.90])
    plt.savefig(f'modal_coeff_hist_{cycle_name.replace(" ", "_").lower()}.png')
    plt.show()
