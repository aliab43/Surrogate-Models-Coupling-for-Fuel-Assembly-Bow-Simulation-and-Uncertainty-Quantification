import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt
from rootlogon import ROOT, DataServer



nominal_hot_watergap = 0.002  # Nominal hot water gap in meters


 # Function to plot assembly deformations
def plot_assembly_deformation(assemblies):
        plt.rcParams['axes.linewidth'] = 2.0
        fig, ax = plt.subplots(1, 1, figsize=(15, 8))
        plt.subplots_adjust(left=0.05, right=0.95, top=0.93, bottom=0.1)
        y = np.arange(1, 11)  # Grid positions (1 to 10)
        # Define colors for assemblies
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
        
        # Plot each assembly
        for i, x in enumerate(assemblies):
            # Reverse the displacement vector so grid 1 is at the bottom and grid 10 at the top
            x_plot = x[::-1]
            plt.xlim(nominal_hot_watergap - nominal_hot_watergap * 1.0, 15 * nominal_hot_watergap + nominal_hot_watergap * 1.0)
            plt.ylim(0., 11.)
            # Set x-ticks to match nominal water gap positions
            xtick_positions = [nominal_hot_watergap + nominal_hot_watergap * k for k in range(15)]
            plt.xticks(xtick_positions, [str(k + 1) for k in range(15)])
            plt.yticks(y, [str(k) for k in range(1, 11)])  # Grid labels
            plt.tick_params(labelleft=True, left=True, labelsize=15)
            plt.tick_params(axis='x', labelsize=15, labelrotation=0)
            plt.ylabel('Grids', fontsize=25)
            plt.xlabel('Assemblies', fontsize=25)
            plt.grid(True, which='major', linestyle='-')
            # Use .get(i+1, 'black') to avoid KeyError
            plt.plot(x_plot, y, marker='o', linestyle='-', color=colors.get(i+1))

            
        # Add the title
        #plt.title("Bow of a line of 15 fuel assemblies", fontsize=22, pad=20)
        
        plt.savefig('deformation0.png')
        plt.close('all')


   

#displacement of the mean deformation after irradiation 
def GramSchmidt():
        f_0 = np.array([0.0, 1.09090909e-04, 3.63636364e-04, 6.72727273e-04, 9.27272727e-04,
                        1.00000000e-03, 8.54545455e-04, 5.63636364e-04, 2.54545455e-04, 0.0])
        f_1 = np.array([0.0, 5.4e-04, 1.0e-03, 9.6e-04, 3.8e-04, -3.8e-04, -8.6e-04,
                        -9.4e-04, -6.2e-04, 0.0])
        f_2 = np.array([0.0, 7.6e-04, 1.0e-03, 2.0e-04, -9.4e-04, -8.6e-04, 1.6e-04,
                        9.2e-04, 9.0e-04, 0.0])
        # Construction of the orthogonal base using the Gram-Schmidt method
        M = np.zeros((3, 10))
        M[0] = f_0 / np.linalg.norm(f_0)
        v = f_1 - np.dot(f_0, f_1) / np.dot(f_0, f_0) * f_0
        M[1] = v / np.linalg.norm(v)
        v = f_2 - np.dot(f_0, f_2) / np.dot(f_0, f_0) * f_0 - np.dot(M[1], f_2) / np.dot(M[1], M[1]) * M[1]
        M[2] = v / np.linalg.norm(v)
        
        return M[0],M[1], M[2]   
    


        

# THE MODAL COEFFICIENTS FOR EACH ASSEMBLY AFTER IRRADIATION PHASE
C_1_a = -2.913489e-03
S_1_a = -4.726097e-04
W_1_a = -5.382505e-04
C_2_a = -7.182974e-03
S_2_a = -9.389901e-04
W_2_a = -1.997472e-03
C_3_a = -1.166399e-02
S_3_a = -1.454465e-03
W_3_a = -2.941596e-03
C_4_a = -1.572742e-02
S_4_a = -1.969498e-03
W_4_a = -3.238798e-03
C_5_a = -1.921459e-02
S_5_a = -2.190367e-03
W_5_a = -2.992801e-03
C_6_a = -1.008510e-02
S_6_a = -4.159922e-03
W_6_a = -1.755898e-03
C_7_a = -2.881014e-03
S_7_a = -6.847987e-04
W_7_a = -8.030385e-04
C_8_a = 4.513329e-04
S_8_a = 2.746055e-04
W_8_a = -1.476861e-03
C_9_a = 4.338604e-03
S_9_a = -6.797427e-04
W_9_a = 2.035914e-04
C_10_a = 1.294234e-02
S_10_a = 5.473753e-03
W_10_a = 2.567530e-03
C_11_a = 1.883300e-02
S_11_a = 2.277920e-03
W_11_a = 3.419408e-03
C_12_a = 1.484919e-02
S_12_a = 1.835079e-03
W_12_a = 2.944957e-03
C_13_a = 1.063202e-02
S_13_a = 1.326535e-03
W_13_a = 2.318231e-03
C_14_a = 6.604120e-03
S_14_a = 8.758312e-04
W_14_a = 1.651382e-03
C_15_a = 2.173880e-03
S_15_a = 3.952360e-04
W_15_a = 4.783874e-04
C_a = []
S_a = []
W_a = []



C_i_j_0 = np.array(C_1_a)
S_i_j_0 = np.array(S_1_a)
W_i_j_0 = np.array(W_1_a)
c, s, w = GramSchmidt()
assemblies = []
for i in range(15):
    x = np.zeros(10)
    for j in range(10):
        x[j] = nominal_hot_watergap + (i) * nominal_hot_watergap
        deformation_C = C_i_j_0[i] * c[j]
        deformation_S = S_i_j_0[i] * s[j]
        deformation_W = W_i_j_0[i] * w[j]
        x[j] += (deformation_C + deformation_S + deformation_W)
    assemblies.append(x)
print(assemblies)
plot_assembly_deformation(assemblies)
variance = []
for i in range(15):
    var_grid_j = np.zeros(10)
    for j in range(10):
        
        var_deformation_C = var_c_i_a[i] * c[j]**2
        var_deformation_S = var_s_i_a[i] * s[j]**2
        var_deformation_W = var_w_i_a[i] * w[j]**2
        var_grid_j[j] += (var_deformation_C + var_deformation_S + var_deformation_W)
    variance.append(var_grid_j)

print("std of the deformation for each grid and each assembly:")
for i in range(15):
    print(f"Assembly {i+1}: {np.max(np.sqrt(variance[i]))}")


# Define assembly indices for each cycle type based on the provided color mapping (1-based to 0-based)
cycle2_indices = [0, 1, 2, 13, 14, 12]      # darkblue: 1,2,3,14,15,13
cycle1_indices = [10, 3, 4, 11]              # darkred: 11,4,5,12
cycle3_indices = [5, 6, 7, 8, 9]             # steelblue: 6,7,8,9,10

cycle_groups = [
    ("Cycle 2", cycle2_indices, 'darkblue'),
    ("Cycle 1", cycle1_indices, 'darkred'),
    ("Cycle 3", cycle3_indices, 'steelblue')
]

def plot_cycle_assemblies_with_uncertainty(cycle_name, indices, color, assemblies, variance):
    plt.rcParams['axes.linewidth'] = 2.0
    fig, ax = plt.subplots(1, 1, figsize=(15, 8))
    plt.subplots_adjust(left=0.05, right=0.95, top=0.93, bottom=0.1)
    y = np.arange(1, 11)
    for i in indices:
        x = assemblies[i]
        x_plot = x[::-1]
        std_disp = np.sqrt(variance[i][::-1]) / 2.65  # Adjust for nominal hot water gap
        ax.plot(x_plot, y, marker='o', linestyle='-', color=color)
        ax.fill_betweenx(
            y,
            x_plot - 2 * std_disp,
            x_plot + 2 * std_disp,
            color=color,
            alpha=0.15
        )
    # Highlight nominal water gap region only between A7 and A8 and grids 5 and 6
    x_assembly_7 = nominal_hot_watergap + nominal_hot_watergap * 6
    x_assembly_8 = nominal_hot_watergap + nominal_hot_watergap * 7
    y_grid_5 = 5
    y_grid_6 = 6
    # Fill between A7 and A8, from grid 5 to grid 6
    ax.fill_betweenx(
        y=[y_grid_5, y_grid_6],
        x1=x_assembly_7,
        x2=x_assembly_8,
        color='deepskyblue',
        alpha=0.4,
        label='Nominal water gap = 2 mm'
    )
    # Draw vertical dashed lines only between grid 5 and 6
    ax.plot([x_assembly_7, x_assembly_7], [y_grid_5, y_grid_6], color='black', linestyle='--', linewidth=2)
    ax.plot([x_assembly_8, x_assembly_8], [y_grid_5, y_grid_6], color='black', linestyle='--', linewidth=2)
    ax.set_xlim(nominal_hot_watergap - nominal_hot_watergap * 1.0, 15 * nominal_hot_watergap + nominal_hot_watergap * 1.0)
    ax.set_ylim(0., 11.)
    xtick_positions = [nominal_hot_watergap + nominal_hot_watergap * k for k in range(15)]
    ax.set_xticks(xtick_positions)
    ax.set_xticklabels([str(k + 1) for k in range(15)])
    ax.set_yticks(y)
    ax.set_yticklabels([str(k) for k in range(1, 11)])
    ax.tick_params(labelleft=True, left=True, labelsize=15)
    ax.tick_params(axis='x', labelsize=15, labelrotation=0)
    ax.set_ylabel('Grids', fontsize=25)
    ax.set_xlabel('Assemblies', fontsize=25)
    ax.grid(True, which='major', linestyle='-')
    ax.set_title(f'Deformation and 95% confidence interval for the Fuel assembly {cycle_name}', fontsize=20)
    # Only show the nominal water gap in the legend
    handles, labels = ax.get_legend_handles_labels()
    new_handles = []
    new_labels = []
    for h, l in zip(handles, labels):
        if l.startswith('Nominal water gap'):
            new_handles.append(h)
            new_labels.append(l)
    if new_handles:
        ax.legend(new_handles, new_labels, fontsize=14)
    plt.savefig(f'deformation_{cycle_name.replace(" ", "_").lower()}.png')
    plt.close('all')

# Plot for each cycle group
for cycle_name, indices, color in cycle_groups:
    plot_cycle_assemblies_with_uncertainty(cycle_name, indices, color, assemblies, variance)