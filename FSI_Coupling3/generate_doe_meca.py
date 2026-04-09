import numpy as np
import matplotlib.pyplot as plt


class GenerateDOEMeca:
    
    def __init__(self):
        
        self.Nb_params = 13

    def generateModal_coef(self):


        sigma_mesurments = 0.0003  # Measurement noise on each grid (in m)
        muC1, muS1, muW1 = 0.0003, 0.0003, 0.0003  # Cycle 1
        muC2, muS2, muW2 = 0.003, 0.002, 0.001
        muC3, muS3, muW3 = 0.005, 0.004, 0.003

        # Define mean vectors for C, S, W (length 15, as specified)
        mean_C = [
            muC2, muC2, muC2, muC1, muC1, muC3, muC3, muC3, muC3, muC3,
            muC3, muC1, muC1, muC2, muC2
        ]
        mean_S = [
            muS2, muS2, muS2, muS1, muS1, muS3, muS3, muS3, muS3, muS3,
            muS3, muS1, muS1, muS2, muS2
        ]
        mean_W = [
            muW2, muW2, muW2, muW1, muW1, muW3, muW3, muW3, muW3, muW3,
            muW3, muW1, muW1, muW2, muW2
        ]

        # Covariance matrix: diagonal with variance sigma_mesurments^2
        cov_diag = sigma_mesurments ** 2
        Sigma = np.diag([cov_diag] * 15)

        # Generate multivariate samples with specified mean and diagonal covariance
        C_i_j_C = np.random.multivariate_normal(mean=mean_C, cov=Sigma).tolist()
        C_i_j_S = np.random.multivariate_normal(mean=mean_S, cov=Sigma).tolist()
        C_i_j_W = np.random.multivariate_normal(mean=mean_W, cov=Sigma).tolist()
        return C_i_j_C, C_i_j_S, C_i_j_W

    def generate_one_sample(self, nb_assembly):

        
        C_vectors, S_vectors, W_vectors = self.generateModal_coef()
        print(C_vectors)
        #F_h_all_calculations, modal_coefficients = self.getFh(nb_assembly+1,Ns)  

        
        
         # Initialize grid_clamping with a default value
        grid_clamping = None  # Or some default value, e.g., 0.0

        # Define the distribution of Grid Clamping with respect to assembly age and first configuration
        if nb_assembly in [3, 4, 11, 12]:
            mean_serrage = 37.0  
            std_dev_serrage = 0.5  
            grid_clamping = np.random.normal(loc=mean_serrage, scale=std_dev_serrage)
        elif nb_assembly in [0, 1, 2, 13, 14, 15]:
            mean_serrage = 10.0  
            std_dev_serrage = 1  
            grid_clamping = np.random.normal(loc=mean_serrage, scale=std_dev_serrage)
        elif nb_assembly in [5, 6, 7, 8, 9, 10]:
            mean_serrage = 5  # N
            std_dev_serrage = 1 # N
            grid_clamping = np.random.normal(loc=mean_serrage, scale=std_dev_serrage)
        else:
            # Handle unexpected values of nb_assembly
            print(f"Warning: nb_assembly {nb_assembly} not recognized. Defaulting grid_clamping.")
            grid_clamping = 32 
        
        

        # Extract deformation and modal coefficients
        def_C = C_vectors[nb_assembly]
        def_S = S_vectors[nb_assembly]
        def_W = W_vectors[nb_assembly]
        f_0, f_C, f_S, f_W = 0,0,0,0 #modal_coefficients[Ns]

        return [def_C, def_S, def_W, f_0, f_C, f_S, f_W,grid_clamping]


   

    def generateDOE_meca(self):
        
        nb_assembly = 15
        nb_params = self.Nb_params
        matrice_resultats = np.zeros((nb_assembly, nb_params))
        # Parameters for the normal distribution of neutron flux
        mean_flux_neutronique = 0.8e+18
        std_dev_flux_neutronique = 0.03 * mean_flux_neutronique  # Variability of +-3%

    
        mean_growth_creep = 1.
        std_dev_growth_creep = 0.3

        # Sampling random values based on specified distributions
        flux_neutronique = np.random.normal(loc=mean_flux_neutronique, scale=std_dev_flux_neutronique)
        temperature_entree = np.random.normal(553.15, 3.33)  # The temperature ~ Normal distribution
        pre_charge_axial = np.random.uniform(0.01, 0.025)  # ~ Uniform distribution
 

        constant_growth = np.random.normal(loc=mean_growth_creep, scale=std_dev_growth_creep)
        constant_creep = np.random.normal(loc=mean_growth_creep, scale=std_dev_growth_creep)
        
        for s in range(nb_assembly):
            file_name = "DOE_for_mechanics/Doe_for_surrogate.txt"
            
            # Utilisez 'a' pour ajouter sans écraser le contenu existant
            with open(file_name, 'a') as file: 

                # Définissez le texte pour l'en-tête seulement si c'est la première itération
                if s == 0:
                    txtFormatDefC = "Def_C | "
                    txtFormatDefS = "Def_S | "
                    txtFormatDefW = "Def_W | "
                    txtFormat0 = "F_0 | "
                    txtFormatC = "F_C | "
                    txtFormatS = "F_S | "
                    txtFormatW = "F_W | "
                    txtFormatD = " D |"
                    
                    # Ecrire l'en-tête du fichier une seule fois
                    file.write("#NAME: tdsCastem \n")
                    file.write("#COLUMN_NAMES: ")

                    file.write(txtFormatDefC)
                    file.write(txtFormatDefS)
                    file.write(txtFormatDefW)
                    file.write(txtFormat0)
                    file.write(txtFormatC)
                    file.write(txtFormatS)
                    file.write(txtFormatW)
                    file.write(f"T_in | fast_flux | MSI | grid_clamping | C_growth | C_creep \n")
                    
                    file.write("#COLUMN_TYPES: ")
                    for z in range(self.Nb_params):
                        file.write(txtFormatD)
                    file.write("\n\n")

                # Ecrire les données générées pour chaque itération
                txtFormat1 = "{price:.8f} "
                txtFloat = "{} "

                # Générer des données pour un échantillon unique
                vect_one = self.generate_one_sample(s)
                
                # Ecrire l'échantillon généré dans le fichier
                file.write(txtFormat1.format(price=vect_one[0]))
                file.write(txtFormat1.format(price=vect_one[1]))
                file.write(txtFormat1.format(price=vect_one[2]))
                file.write(txtFormat1.format(price=vect_one[3]))
                file.write(txtFormat1.format(price=vect_one[4]))
                file.write(txtFormat1.format(price=vect_one[5]))
                file.write(txtFormat1.format(price=vect_one[6]))
                file.write(txtFloat.format(temperature_entree))
                file.write(txtFloat.format(flux_neutronique))
                file.write(txtFloat.format(pre_charge_axial))
                file.write(txtFloat.format(vect_one[7]))
                file.write(txtFloat.format(constant_growth))
                file.write(txtFloat.format(constant_creep))
                file.write("\n")
            matrice_resultats[s, :] = [vect_one[0], vect_one[1], vect_one[2], vect_one[3], vect_one[4], vect_one[5],
                                            vect_one[6], temperature_entree, flux_neutronique, pre_charge_axial, vect_one[7], constant_growth,
                                            constant_creep]

        # Sauvegarder la matrice dans un fichier .npy
        npy_file_name = "DOE_for_mechanics/Doe_for_mechanics.npy"
        np.save(npy_file_name, matrice_resultats)
        print(matrice_resultats)
        print(f"Matrice sauvegardée dans le fichier {npy_file_name}")

s = GenerateDOEMeca().generateDOE_meca()




















"""    def getFh(self, nb_assembly, Ns):
        
        if nb_assembly > 15 :
            nb_assembly = 15
        # Define the 4 main modes for modal decomposition (size 10 for each mode)
        Q = np.array([
            [0, -0.34182351, -0.33986812, -0.35552546, -0.36900858, -0.37084271, -0.36102097, -0.34508458, -0.3437104, 0],  # F_0
            [0, -0.55937995, -0.39424153, -0.06243252, 0.30706583, 0.46809913, 0.36218138, 0.07404396, -0.27875758, 0],   # F_C
            [0, -0.00479597, 0.31635535, 0.46610865, 0.36727189, 0.06233596, -0.26210788, -0.47490219, -0.49963184, 0],   # F_S
            [0, 0.19894744, -0.4367867, -0.15343436, 0.43733209, 0.20649579, -0.39132532, -0.36346723, 0.47639416, 0]     # F_W
        ]).T  # Transpose Q so each column is a mode

        # Initialize lists for each grid
        F_h_G_2, F_h_G_3, F_h_G_4, F_h_G_5, F_h_G_6, F_h_G_7, F_h_G_8, F_h_G_9 =  [], [], [], [], [], [], [], []

        # Loop through the grids (from 2 to 9) to retrieve all F_h values
        for index in range(2, 10):
            hydrau_results = f"results_day_7_hydrau/predicted_F_h_G_{index}_A_{nb_assembly}.dat"

            with open(hydrau_results, 'r') as file:
                for line in file:
                    # Ignore comment lines
                    if line.startswith("#"):
                        continue

                    # Extract data if the line is not empty
                    if line.strip():
                        row = list(map(float, line.split()))
                        F_h_value = row[56]  # We take column 56 (F_h)

                        # Store in the correct list based on the grid number
                        if index == 2:
                            F_h_G_2.append(F_h_value)
                        elif index == 3:
                            F_h_G_3.append(F_h_value)
                        elif index == 4:
                            F_h_G_4.append(F_h_value)
                        elif index == 5:
                            F_h_G_5.append(F_h_value)
                        elif index == 6:
                            F_h_G_6.append(F_h_value)
                        elif index == 7:
                            F_h_G_7.append(F_h_value)
                        elif index == 8:
                            F_h_G_8.append(F_h_value)
                        elif index == 9:
                            F_h_G_9.append(F_h_value)

        # Now that each list F_h_G_2 to F_h_G_9 is filled, we will construct a global vector
        F_h_all_calculations = []  # List to store the global vector for each calculation

        # Iterate through each calculation row (there are 1000 rows of calculations for each grid)
        for Ns in range(len(F_h_G_2)):  
            # Construct the global vector for each calculation row
            F_h_vector = [0,
                F_h_G_2[Ns], F_h_G_3[Ns], F_h_G_4[Ns], F_h_G_5[Ns], 
                F_h_G_6[Ns], F_h_G_7[Ns], F_h_G_8[Ns], F_h_G_9[Ns], 0
            ]
            
            # Add this global vector to the list of all calculations
            F_h_all_calculations.append(F_h_vector)

        # Now that we have F_h_all_calculations, we can calculate the modal coefficients
        modal_coefficients_list = []  # List to store modal coefficients for each calculation

        # Loop through each global F_h vector to calculate the modal coefficients
        for F_h in F_h_all_calculations:
            # Calculate the modal coefficients c_0, c_C, c_S, c_W via least squares decomposition
            coefficients, _, _, _ = np.linalg.lstsq(Q, F_h, rcond=None)  # Solve Q * c = F_h
            
            # Store the modal coefficients for this calculation
            modal_coefficients_list.append(coefficients)

        # Return both the F_h vectors for each calculation and the modal coefficients
        return F_h_all_calculations, modal_coefficients_list"""










"""def get_coef(self, file_name):
        C_vectors = []
        S_vectors = []
        W_vectors = []

        with open(file_name, 'r') as file:
            for line in file:
                # Ignore comment lines
                if line.startswith("#"):
                    # Extract column names
                    if line.startswith("#COLUMN_NAMES"):
                        columns = line.split(":")[1].strip().split(" | ")
                    continue  # Move to the next lines

                # Extract data if the line is not empty
                if line.strip():
                    # Convert the line to a list of floating-point numbers
                    row = list(map(float, line.split()))

                    # Extract the C, S, W parts and add them to the respective vectors
                    C_vectors.append(row[:16])  # Columns C_1 to C_16
                    S_vectors.append(row[16:32])  # Columns S_1 to S_16
                    W_vectors.append(row[32:48])  # Columns W_1 to W_16

        return C_vectors, S_vectors, W_vectors"""





"""sigma_C_0 = 0.0002  # C_0 is a new Fuel assembly
        sigma_S_0 = 0.0001
        sigma_W_0 = 0.00005
        sigma_C_1 = 0.002  # C_1 is cycle 1
        sigma_S_1 = 0.001
        sigma_W_1 = 0.0005
        sigma_C_23 = 0.005  # C_23 is a fuel assembly in cycle 2 or 3
        sigma_S_23 = 0.003
        sigma_W_23 = 0.001
        # Variance vectors for the components (we are storing variances here, not standard deviations)
        C_i = [sigma_C_1**2] * 3 + [sigma_C_0**2] * 2 + [sigma_C_23**2] * 6 + [sigma_C_0**2] * 2 + [sigma_C_1**2] * 3
        S_i = [sigma_S_1**2] * 3 + [sigma_S_0**2] * 2 + [sigma_S_23**2] * 6 + [sigma_S_0**2] * 2 + [sigma_S_1**2] * 3
        W_i = [sigma_W_1**2] * 3 + [sigma_W_0**2] * 2 + [sigma_W_23**2] * 6 + [sigma_W_0**2] * 2 + [sigma_W_1**2] * 3
        # Convert to numpy arrays
        Ci = np.array(C_i)
        Si = np.array(S_i)
        Wi = np.array(W_i)
        # Initial covariance matrix (diagonal)
        Sigma_C = np.diag(Ci)
        Sigma_S = np.diag(Si)
        Sigma_W = np.diag(Wi)
        # Matrix A (vector of 1s because the sum of components must be zero)
        A = np.ones((1, len(C_i)))
        # Inverse matrix calculations for conditioning
        inv_A_Sigma_C_A_T = np.linalg.inv(A @ Sigma_C @ A.T)
        inv_A_Sigma_S_A_T = np.linalg.inv(A @ Sigma_S @ A.T)
        inv_A_Sigma_W_A_T = np.linalg.inv(A @ Sigma_W @ A.T)
        # Conditional covariance matrix for each component
        M1 = Sigma_C - (Sigma_C @ A.T @ inv_A_Sigma_C_A_T @ A @ Sigma_C)
        M2 = Sigma_S - (Sigma_S @ A.T @ inv_A_Sigma_S_A_T @ A @ Sigma_S)
        M3 = Sigma_W - (Sigma_W @ A.T @ inv_A_Sigma_W_A_T @ A @ Sigma_W)
        # Generating conditionally multivariate normal distributed samples for C_i, S_i, and W_i
        C_i_j_C = np.random.multivariate_normal(mean=np.zeros(16), cov=M1).tolist()
        C_i_j_S = np.random.multivariate_normal(mean=np.zeros(16), cov=M2).tolist()
        C_i_j_W = np.random.multivariate_normal(mean=np.zeros(16), cov=M3).tolist()"""