import numpy as np
import matplotlib.pyplot as plt
import param_scalar 
from rootlogon import ROOT, DataServer

class GenerateDOEMeca:
    
    def __init__(self):
        
        self.Nb_params = 13

    

    def generate_one_sample(self, nb_assembly):

        
        indexC = param_scalar.index_C
        indexS = param_scalar.index_S
        indexW = param_scalar.index_W
        nb_assemblage = 15

        # Lecture de la matrice du sample d'origine 
        tds = DataServer.TDataServer()   ########"TDS_Sample", "DoE for Sample")
        file_name =  "Matrix_M_N.dat"
        tds.fileDataRead(file_name)
        nS = tds.getNPatterns()
        #print(nS)

        txtFormat1 =  tds.getIteratorName() + "=={price:d}"
        #print(txtFormat1)
        #print(txtFormat1.format(price=indexC))
        C_i_j_C_mat = tds.getMatrix("Coeff_C", txtFormat1.format(price=indexC))
        #print(C_i_j_C_mat)
        C_i_j_S_mat = tds.getMatrix("Coeff_S", txtFormat1.format(price=indexS))
        C_i_j_W_mat = tds.getMatrix("Coeff_W", txtFormat1.format(price=indexW))
        #print(C_i_j_S_mat)
        #print(C_i_j_W_mat)

        # Create the ndarray with the good shape
        C_i_j_C = np.frombuffer(C_i_j_C_mat.GetMatrixArray(), dtype=np.float64,count=nb_assemblage).reshape(1, nb_assemblage)
        
        C_i_j_S = np.frombuffer(C_i_j_S_mat.GetMatrixArray(), dtype=np.float64,count=nb_assemblage).reshape(1, nb_assemblage)
        
        C_i_j_W = np.frombuffer(C_i_j_W_mat.GetMatrixArray(), dtype=np.float64,count=nb_assemblage).reshape(1, nb_assemblage)
        

        tds_coeff = "coef_C_S_W.dat"
        

                
        with open(tds_coeff, 'w') as file:
                    file.write("#COLUMN_NAMES:  C |S | W | iter")
                    file.writelines("\n\n")
                    for i in range(nb_assemblage):
                        file.write( "{} {} {} {} \n".format(C_i_j_C[0][i],  C_i_j_S[0][i],  C_i_j_W[0][i],i))

        data = np.loadtxt(tds_coeff)
        C_i_j_C = data[:, 0]
        C_i_j_S = data[:, 1]
        C_i_j_W = data[:, 2]
        
        
        
        
        
        # Initialize grid_clamping with a default value
        grid_clamping = param_scalar.grid_clamping

        # Define the distribution of Grid Clamping with respect to assembly age and first configuration
        if nb_assembly in [3, 4, 11, 12]:
           
            grid_clamping = grid_clamping 
        elif nb_assembly in [0, 1, 2, 13, 14, 15]:
            
            grid_clamping = grid_clamping * 0.27 
        elif nb_assembly in [5, 6, 7, 8, 9, 10]:
       
            grid_clamping =  grid_clamping * 0.13

        else:
            # Handle unexpected values of nb_assembly
            print(f"Warning: nb_assembly {nb_assembly} not recognized. Defaulting grid_clamping.")
            grid_clamping = 37
        
        

        # Extract deformation and modal coefficients
        def_C = C_i_j_C[nb_assembly]
        def_S = C_i_j_S[nb_assembly]
        def_W = C_i_j_W[nb_assembly]


        f_0, f_C, f_S, f_W = 0,0,0,0 #modal_coefficients[Ns]

        return [def_C, def_S, def_W, f_0, f_C, f_S, f_W,grid_clamping]


   

    def generateDOE_meca(self):
        
        nb_assembly = 15
        nb_params = self.Nb_params
        matrice_resultats = np.zeros((nb_assembly, nb_params))
        flux_neutronique = 0.8e+18
        temperature_entree =  553 # to be changed
        pre_charge_axial = param_scalar.pre_charge_axial
 

        constant_growth = param_scalar.constant_growth 
        constant_creep = param_scalar.constant_creep 
        
        for s in range(nb_assembly):
            
            
            

            # Générer des données pour un échantillon unique
            vect_one = self.generate_one_sample(s)
                
               
            matrice_resultats[s, :] = [vect_one[0], vect_one[1], vect_one[2], vect_one[3], vect_one[4], vect_one[5],
                                            vect_one[6], temperature_entree, flux_neutronique, pre_charge_axial, vect_one[7], constant_growth,
                                            constant_creep]

        # Sauvegarder la matrice dans un fichier .npy
        npy_file_name = "Doe_for_mechanics.npy"
        np.save(npy_file_name, matrice_resultats)
        #print(matrice_resultats)
        #print(f"Matrice sauvegardée dans le fichier {npy_file_name}")

s = GenerateDOEMeca().generateDOE_meca()









