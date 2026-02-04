import numpy as np
import matplotlib.pyplot as plt

nominal_hot_watergap = 0.002  # Nominal hot water gap in meters

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

# === AJOUT : Définition des statistiques à partir des commentaires ===
mean_C = np.array([-0.0031455, -0.00738935, -0.01167397, -0.01566438, -0.01924809, -0.01013561,
                   -0.00298925, 0.00049625, 0.00430913, 0.01305934, 0.01892787, 0.0148245,
                   0.01060557, 0.00657104, 0.00224777])
std_C = np.array([1.06916392e-04, 7.16050291e-05, 5.56956866e-05, 5.79815488e-05,
                  6.70085437e-05, 1.57490140e-04, 2.14054423e-04, 1.97275405e-04,
                  2.26069801e-04, 2.39776181e-04, 3.79869953e-05, 4.17768207e-05,
                  8.18487961e-05, 9.44770643e-05, 6.82994602e-05])
mean_S = np.array([-0.00048305, -0.00094681, -0.00147281, -0.00196022, -0.00209924, -0.00416291,
                   -0.00067327, 0.0002572, -0.0006619, 0.00546541, 0.0023652, 0.00189862,
                   0.00139252, 0.00091244, 0.00041147])
std_S = np.array([8.83623230e-06, 1.30731126e-05, 2.04672366e-05, 2.42014719e-05,
                  1.05357864e-05, 1.26835074e-05, 3.53580537e-05, 3.89226920e-05,
                  3.58606155e-05, 3.30593728e-05, 5.28488636e-06, 1.94565629e-05,
                  1.65848527e-05, 1.33529380e-05, 1.35207312e-05])
mean_W = np.array([-0.00058913, -0.00211514, -0.00299628, -0.0033701, -0.00327792, -0.00176462,
                   -0.00082146, -0.00146316, 0.00019522, 0.00258906, 0.00352268, 0.00296832,
                   0.00229252, 0.00171676, 0.00048564])
std_W = np.array([1.84549639e-05, 2.00308644e-05, 4.76084734e-05, 2.86683021e-05,
                  4.29087087e-05, 2.66454696e-05, 3.64311214e-05, 3.60448786e-05,
                  3.45671561e-05, 4.46757051e-05, 1.32659046e-05, 2.92425646e-05,
                  3.40272654e-05, 3.10801085e-05, 1.88712895e-05])
q05_C = np.array([-0.00328632, -0.00749254, -0.01178487, -0.01575135, -0.01934283, -0.01059046,
                  -0.00354463, 0.00030245, 0.00369154, 0.01294093, 0.01888385, 0.01476707,
                  0.01051104, 0.00643446, 0.00208254])
q95_C = np.array([-0.00303453, -0.00728527, -0.011605, -0.01556849, -0.01915247, -0.01007957,
                  -0.00288352, 0.00100607, 0.00470888, 0.01349279, 0.01901787, 0.01492389,
                  0.01076517, 0.00667697, 0.00234224])
q05_S = np.array([-0.00049681, -0.00096673, -0.00150076, -0.00198854, -0.00211695, -0.0041981,
                  -0.00070178, 0.00017159, -0.00067992, 0.00537048, 0.00235245, 0.00185353,
                  0.00135949, 0.00088628, 0.00039632])
q95_S = np.array([-0.0004697, -0.00092286, -0.0014278, -0.00189899, -0.00207641, -0.00415938,
                  -0.00057891, 0.00027467, -0.00057469, 0.00551276, 0.0023715, 0.00191966,
                  0.0014176, 0.00093061, 0.00043469])
q05_W = np.array([-6.17482229e-04, -2.15838556e-03, -3.03528307e-03, -3.39310891e-03,
                  -3.30108798e-03, -1.84220845e-03, -9.15276059e-04, -1.47691603e-03,
                  9.51415186e-05, 2.56735313e-03, 3.51369682e-03, 2.88550299e-03,
                  2.20481973e-03, 1.63903302e-03, 4.38952377e-04])
q95_W = np.array([-0.00056123, -0.0020872, -0.00286421, -0.00327892, -0.00315587, -0.00175519,
                  -0.00080336, -0.00136729, 0.0002418, 0.00267823, 0.00356062, 0.0029929,
                  0.00233175, 0.00175162, 0.00050814])
# === FIN DU CHARGEMENT DES STATISTIQUES ===

c, s, w = GramSchmidt()
assemblies = []
assemblies_std = []
assemblies_q05 = []
assemblies_q95 = []
for i in range(15):
    x = np.zeros(10)
    std_x = np.zeros(10)
    x_q05 = np.zeros(10)
    x_q95 = np.zeros(10)
    for j in range(10):
        x[j] = nominal_hot_watergap + (i) * nominal_hot_watergap
        deformation_C = mean_C[i] * c[j]
        deformation_S = mean_S[i] * s[j]
        deformation_W = mean_W[i] * w[j]
        x[j] += (deformation_C + deformation_S + deformation_W)
        # Uncertainty propagation (sum of variances)
        var_deformation = (std_C[i] * c[j])**2 + (std_S[i] * s[j])**2 + (std_W[i] * w[j])**2
        std_x[j] = np.sqrt(var_deformation)*2
        # Quantile propagation (borne inf et sup)
        deformation_C_q05 = q05_C[i] * c[j]
        deformation_C_q95 = q95_C[i] * c[j]
        # Pour S et W, utilisez les quantiles si disponibles, sinon mean ± 2*std
        deformation_S_q05 = q05_S[i] * s[j]
        deformation_S_q95 = q95_S[i] * s[j]
        deformation_W_q05 = q05_W[i] * w[j]
        deformation_W_q95 = q95_W[i] * w[j]
        x_q05[j] = nominal_hot_watergap + (i) * nominal_hot_watergap + (deformation_C_q05 + deformation_S_q05 + deformation_W_q05)
        x_q95[j] = nominal_hot_watergap + (i) * nominal_hot_watergap + (deformation_C_q95 + deformation_S_q95 + deformation_W_q95)
    assemblies.append(x)
    assemblies_std.append(std_x)
    assemblies_q05.append(x_q05)
    assemblies_q95.append(x_q95)
    print(assemblies_std)
# --- PLOT AVEC MOYENNE ET ±2*STD ---
plt.rcParams['axes.linewidth'] = 2.0
fig, ax = plt.subplots(1, 1, figsize=(15, 8))
plt.subplots_adjust(left=0.05, right=0.95, top=0.93, bottom=0.1)
y = np.arange(1, 11)
# Correction : 15 couleurs pour 15 assemblées
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
# Ajout : calcul et affichage du std max
std_max = np.max([np.max(std) for std in assemblies_std])
ax.text(
    0.01, 0.98,
    f"Maximal Standard Deviation in [mm] = {std_max*1000:.2f}",
    transform=ax.transAxes,
    fontsize=16,
    verticalalignment='top',
    horizontalalignment='left',
    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none')
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
ax.set_title("Coupled Surrogate Models Final Uncertainty  (mean ± 2 std)", fontsize=20)
plt.savefig('final_uncertainty_all_assemblies_std.png')
plt.close()

# --- PLOT AVEC QUANTILES (5% - 95%) ---
fig, ax = plt.subplots(1, 1, figsize=(15, 8))
plt.subplots_adjust(left=0.05, right=0.95, top=0.93, bottom=0.1)
for i in range(15):
    x = assemblies[i]
    x_q05 = assemblies_q05[i]
    x_q95 = assemblies_q95[i]
    x_plot = x[::-1]
    q05_plot = x_q05[::-1]
    q95_plot = x_q95[::-1]
    ax.plot(x_plot, y, marker='o', linestyle='-', color=colors[i], label=f'A{i+1}' if i < 10 else None)
    ax.fill_betweenx(
        y,
        q05_plot,
        q95_plot,
        color=colors[i],
        alpha=0.15
    )
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
ax.set_title("Coupled Surrogate Models Final Uncertainty  (empirical 5%-95% quantiles)", fontsize=20)
plt.savefig('final_uncertainty_all_assemblies_quantiles.png')
plt.close()