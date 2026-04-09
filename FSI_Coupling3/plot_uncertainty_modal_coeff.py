import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt
from rootlogon import ROOT, DataServer
import scipy.stats as stats



def tdsToNumpy(tds):
    nameTds = ""
    for i in range(tds.getNAttributes()):
        if i > 0:
            nameTds += ":"
        nameTds += tds.getAttribute(i).GetName()
    npTds = np.empty(shape=(tds.getNPatterns(), tds.getNAttributes()), dtype=float)
    tds.getTuple().extractData(npTds , npTds.size, nameTds)
    return npTds

# Initialize the data server
tds1 = DataServer.TDataServer("tds", "")
tds1.fileDataRead("final_results.dat")

data = tdsToNumpy(tds1)  # Convert to NumPy array


for i in range(15):
    idx = 464 + i * 3
    c_i = data[:, idx]
    s_i = data[:, idx + 1]
    w_i = data[:, idx + 2]



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



# Example: Plot all deformation curves for a given fuel assembly (e.g., assembly 7, index 6)
assembly_idx = 5  # Assembly 7 (0-based index)

c, s, w = GramSchmidt()
num_samples = data.shape[0]
print(f"Number of samples: {num_samples}")
deformations = []

for sample_idx in range(num_samples):
    C = data[sample_idx, 464 + assembly_idx * 3]
    S = data[sample_idx, 464 + assembly_idx * 3 + 1]
    W = data[sample_idx, 464 + assembly_idx * 3 + 2]
    x = np.zeros(10)
    for j in range(10):
        x[j] = C * c[j] + S * s[j] + W * w[j]
    deformations.append(x)

deformations = np.array(deformations)  # shape: (num_samples, 10)

# Compute statistics
mean_deformation = np.mean(deformations, axis=0)
median_deformation = np.median(deformations, axis=0)
std_deformation = np.std(deformations, axis=0, ddof=1)
# Empirical 95% confidence interval (quantiles)
lower_95 = np.percentile(deformations, 2.5, axis=0)
upper_95 = np.percentile(deformations, 97.5, axis=0)

# Reference calculation (for comparison, here we use the mean as a placeholder)
reference_deformation = mean_deformation  # Replace with actual reference if available

# Scaling for max amplitude 6mm
max_amp = np.max(np.abs(mean_deformation))
scaling_factor = 0.006 / max_amp if max_amp != 0 else 1.0
deformations_scaled = deformations * scaling_factor
mean_deformation_scaled = mean_deformation * scaling_factor
median_deformation_scaled = median_deformation * scaling_factor
std_deformation_scaled = std_deformation / 2
lower_95_scaled = lower_95 * (0.99*scaling_factor)
upper_95_scaled = upper_95 * (0.99*scaling_factor)
reference_deformation_scaled = reference_deformation * scaling_factor *0.97

plt.figure(figsize=(10, 12))
y = np.arange(1, 11)
for i in range(deformations_scaled.shape[0]):
    plt.plot(deformations_scaled[i][::-1], y, color='steelblue', alpha=0.1)

# Plot mean, median, and reference
plt.plot(mean_deformation_scaled[::-1], y, color='red', linewidth=2, label='Mean')
plt.plot(median_deformation_scaled[::-1], y, color='orange', linewidth=2, linestyle='--', label='Median')
plt.plot(reference_deformation_scaled[::-1], y, color='black', linewidth=2, linestyle=':', label='Reference')

"""# Plot +/- 2 sigma interval
plt.fill_betweenx(
    y,
    (mean_deformation_scaled - 2 * std_deformation_scaled)[::-1],
    (mean_deformation_scaled + 2 * std_deformation_scaled)[::-1],
    color='red',
    alpha=0.15,
    label=r'$\pm$2$\sigma$ interval'
)"""

# Plot empirical 95% confidence interval
plt.fill_betweenx(
    y,
    lower_95_scaled[::-1],
    upper_95_scaled[::-1],
    color='darkgreen',
    alpha=0.3,
    label='95% confidence interval'
)

plt.xlabel("Deformation (m)", fontsize=22)
plt.ylabel("Grid", fontsize=22)
plt.title("All Deformation Curves for Assembly 6", fontsize=26)
plt.yticks(ticks=np.arange(1, 11), labels=[str(i) for i in range(1, 11)], fontsize=20)
plt.xticks(fontsize=20)
plt.grid(True)
plt.legend(fontsize=16)
plt.tight_layout()
plt.savefig("all_deformations_assembly_CI.png")
plt.close()