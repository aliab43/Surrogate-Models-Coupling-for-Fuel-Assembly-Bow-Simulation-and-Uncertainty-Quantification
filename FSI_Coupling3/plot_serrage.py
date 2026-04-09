import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt
from rootlogon import ROOT, DataServer



# Données fournies
flue1 = np.array([0, 0.1e25, 0.2e25, 0.3e25, 0.4e25, 0.5e25, 0.7e25, 2e25, 4e25, 1e26])
fser1 = np.array([1, 0.65, 0.49, 0.35, 0.30, 0.24, 0.19, 0.11, 0.04, 0.01]) # * 40

# Interpolation linéaire par morceaux
interp_func = interp1d(flue1, fser1, kind='linear', bounds_error=False, fill_value=(fser1[0], fser1[-1]))

# Génération des valeurs interpolées pour l'affichage
flue_fit = np.linspace(min(flue1), max(flue1), 50)
fser_fit_interp = interp_func(flue_fit)

# Fonction pour trouver la fluence pour une force de serrage donnée
def find_fluence_for_serrage(serrage):
    # Fonction de coût pour minimiser la différence entre l'interpolation et la valeur de serrage
    cost_function = lambda x: np.abs(interp_func(x) - serrage)
    
    # Effectuer la minimisation dans la plage de données
    result = minimize_scalar(cost_function, bounds=(min(flue1), max(flue1)), method='bounded')
    
    # Retourner la fluence trouvée
    return result.x


# Function to convert TDataServer data to a NumPy matrix
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
tds1 = DataServer.TDataServer("tdsObs", "observations")
tds1.fileDataRead("Doe_for_surrogate.dat")

data = tdsToNumpy(tds1)
# Convertir les données en tableau numpy
serrages_initiaux = data[:,10]

# Incrément mensuel de fluence
incremente_fluence = 0.277e25  # Exemple d'incrément mensuel

# Pour stocker les résultats
fluences_initiales = []
fluences_maj = []
serrages_majs = []

# Traitement pour chaque valeur de serrage initiale
for serrage_initial in serrages_initiaux:
    # Étape 1 : Trouver la fluence correspondante à la valeur de serrage initiale
    fluence_initiale = find_fluence_for_serrage(serrage_initial)
    fluences_initiales.append(fluence_initiale)
    
    # Étape 2 : Ajouter un incrément à la fluence
    fluence_maj = fluence_initiale + incremente_fluence
    fluences_maj.append(fluence_maj)
    
    # Étape 3 : Remettre à jour la valeur de serrage correspondante après un mois
    serrage_maj = interp_func(fluence_maj)
    serrages_majs.append(serrage_maj)

# Convertir les listes en tableaux NumPy pour des traitements ultérieurs
fluences_initiales = np.array(fluences_initiales)
fluences_maj = np.array(fluences_maj)
serrages_majs = np.array(serrages_majs)

# Displaying the results for the first 10 elements for verification
print("Example of the first 10 results:")
for i in range(3):
    print(f"Initial Grid Clamping: {serrages_initiaux[i]}, "
          f"Initial fluence: {fluences_initiales[i]:.2e}, "
          f"Fluence after one month: {fluences_maj[i]:.2e}, "
          f"Updated Grid Clamping: {serrages_majs[i]:.2f}")

# Graphical display of the first 10 results for verification
plt.figure(figsize=(10, 6))
plt.plot(flue1, fser1, 'o', label='Data')
plt.plot(flue_fit, fser_fit_interp, '-', label='Piecewise Linear Interpolation')

# Marking the initial and updated points for the first 10 Grid Clampings
for i in range(3):
    plt.plot(fluences_initiales[i], serrages_initiaux[i], 'xr', label='Initial Grid Clamping' if i == 0 else "")
    plt.plot(fluences_maj[i], serrages_majs[i], 'xb', label='Updated Grid Clamping' if i == 0 else "")

plt.xlabel('Fluence (x 10^25)')
plt.ylabel('Grid Clamping Force fraction')
plt.title('Linear Interpolation and Monthly Evolution')
plt.grid(True)
plt.legend()
plt.savefig("Grid Clamping_curve_interpolation_monthly_update.png")
plt.show()

data[:,10] = serrages_majs



print(data[:,10])

extra_column = np.arange(2000,dtype=int).reshape(-1, 1)



# Add the extra column to the data matrix
data_with_extra_column = np.hstack((data, extra_column))

# Print the shape of the updated data matrix
print(np.shape(data_with_extra_column))

# Print the new column to verify
print(data_with_extra_column[:, -1])

format_list = ['%.20f'] * data.shape[1] + ['%d']  # Last column is formatted as integer

np.savetxt('data_matrix.txt', data_with_extra_column, delimiter=' ', fmt=format_list)





























"""x_new = x_initial + increment

# Calcul des nouvelles valeurs de y après 1 mois d'irradiation
y_new = exp_decreasing(x_new, a, b, c)

# Mettre à jour la colonne y dans le tableau de données
data[:, 10] = y_new

extra_column = np.arange(2000,dtype=int).reshape(-1, 1)



# Add the extra column to the data matrix
data_with_extra_column = np.hstack((data, extra_column))

# Print the shape of the updated data matrix
print(np.shape(data_with_extra_column))

# Print the new column to verify
print(data_with_extra_column[:, -1])

format_list = ['%.20f'] * data.shape[1] + ['%d']  # Last column is formatted as integer

np.savetxt('data_matrix.txt', data_with_extra_column, delimiter=' ', fmt=format_list)


"""