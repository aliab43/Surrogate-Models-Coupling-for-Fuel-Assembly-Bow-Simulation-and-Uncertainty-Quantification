import numpy as np
import matplotlib.pyplot as plt

class GenerateDoeHydrau:

    def __init__(self, Nb_params):
        self.Nb_params = Nb_params  # 57




    def generate_one_sample(self):


        temperature = np.random.normal(300., 3.33)  # Temperature ~ Normal distribution
        v_m = np.random.normal(5, 0.05)  # Mean velocity ~ Normal distribution
        max_deviation_in = np.random.normal(0.05, 0.005)  # ~ Normal distribution
        l_offset_in = np.random.normal(0.0, 0.02)  # ~ Normal distribution
        max_deviation_out = np.random.normal(0.04, 0.004)  # ~ Normal distribution
        l_offset_out = np.random.normal(0.0, 0.01)  # ~ Normal distribution
        Cg_sensi = np.random.uniform(1.0, 1.232)  # Sensitivity factor ~ Uniform distribution
        hl_sensi = np.random.uniform(0.01, 0.03)  # Heat load sensitivity ~ Uniform distribution
        return [temperature, v_m, max_deviation_in, l_offset_in, max_deviation_out, l_offset_out, Cg_sensi, hl_sensi]




    def generateDoeHydrau(self, Ns, nb_assembly):
        
        
        file_name =  "DOE_for_surrogate_hydraulic.dat"

        matrice_resultats = []

        with open(file_name, 'w') as file : 
            
            txtFormatC = "C_{price:d} | "
            txtFormatS = "S_{price:d} | "
            txtFormatW = "W_{price:d} | "
            txtFormatD = " D |"
            file.write("#NAME: tds_Phorcys \n")
            file.write("#COLUMN_NAMES: ")
            for s in range(nb_assembly):
                file.write(txtFormatC.format(price = s + 1))
            for s in range(nb_assembly):
                file.write(txtFormatS.format(price = s + 1))
            for s in range(nb_assembly):
                file.write(txtFormatW.format(price = s + 1))
            
            file.write("T | v_m | max_deviation_in | l_offset_in | max_deviation_out | l_offset_out | Cg_sensi | hl_sensi| tdsEstim_hydrau__n__iter__ \n")
            file.write("#COLUMN_TYPES: ")
            for z in range(self.Nb_params):
                file.write(txtFormatD)
            file.write(" \n")
            file.write(" \n")
            txtFormat1 = "{price:.8f} "
            txtFormat2 = "{price:.8f} "
            txtFloat = "{} "
            for i in range(Ns):
                vect_one = self.generate_one_sample()
                
                """for j in range(nb_assembly):
                    file.write(txtFormat1.format(price = vect_one[0][j]))
            
                
                for j in range(nb_assembly):
                    file.write(txtFormat1.format(price = vect_one[1][j]))
            
                
                for j in range(nb_assembly):
                    file.write(txtFormat1.format(price = vect_one[2][j]))"""
            

                file.write(txtFloat.format(vect_one[0]))
                file.write(txtFloat.format(vect_one[1]))
                file.write(txtFloat.format(vect_one[2]))
                file.write(txtFloat.format(vect_one[3]))
                file.write(txtFloat.format(vect_one[4]))
                file.write(txtFloat.format(vect_one[5]))
                file.write(txtFloat.format(vect_one[6]))
                file.write(txtFloat.format(vect_one[7]))
                file.write(txtFloat.format(1.0*i))

            

                file.write("\n")
                matrice_resultats.append([vect_one[0], vect_one[1], vect_one[2], vect_one[3], vect_one[4],vect_one[5], vect_one[6],vect_one[7]])

        # Sauvegarder la matrice dans un fichier .npy
        npy_file_name = "DOE_for_mechanics/Doe_for_hydraulic.npy"
        print(matrice_resultats)
        np.save(npy_file_name, matrice_resultats)
        print(f"Matrice sauvegardée dans le fichier {npy_file_name}")


s = GenerateDoeHydrau(Nb_params=57).generateDoeHydrau(Ns=1, nb_assembly=16)


             
