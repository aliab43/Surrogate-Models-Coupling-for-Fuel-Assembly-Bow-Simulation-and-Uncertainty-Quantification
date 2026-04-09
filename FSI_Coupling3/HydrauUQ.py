import numpy as np
import matplotlib.pyplot as plt
from rootlogon import ROOT, DataServer
import param_scalar
class GenerateDoeHydrau:

    def generate_one_sample(self):


        temperature = 300 
        
        hl_sensi = param_scalar.hl_sensi
        indexHBC = param_scalar.indexHBC

        tds = DataServer.TDataServer() 
        file_name =  "Matrix_M_N.dat"
        tds.fileDataRead(file_name)
        nS = tds.getNPatterns()
        #print(nS)

        txtFormat1 =  tds.getIteratorName() + "=={price:d}"
        #print(txtFormat1)
        #print(txtFormat1.format(price=indexHBC))
        Hydraulic_BC = tds.getMatrix("Hydraulic_BC", txtFormat1.format(price=indexHBC))
       

        # Create the ndarray with the good shape
        Hydraulic_BC = np.frombuffer(Hydraulic_BC.GetMatrixArray(), dtype=np.float64,count=5).reshape(1, 5)
       

        tds_hydrau = "Hydrau_BC.dat"
        

                
        with open(tds_hydrau, 'w') as file:
                    file.write("#COLUMN_NAMES:  Hydraulic_BC | iter")
                    file.writelines("\n\n")
                    for i in range(5):
                        file.write( "{} {}\n".format(Hydraulic_BC[0][i], i))

        data = np.loadtxt(tds_hydrau)
        
        v_m = data[0,0]
        max_deviation_in = data[1,0]
        l_offset_in = data[2,0]
        max_deviation_out = data[3,0]
        l_offset_out = data[4,0]
        Cg_sensi = 1.12
        return [temperature, v_m, max_deviation_in, l_offset_in, max_deviation_out, l_offset_out, Cg_sensi, hl_sensi]




    def generateDoeHydrau(self):

        matrice_resultats = []

        
        
        vect_one = self.generate_one_sample()
                
            
        matrice_resultats.append([vect_one[0], vect_one[1], vect_one[2], vect_one[3], vect_one[4],vect_one[5], vect_one[6],vect_one[7]])

        # Sauvegarder la matrice dans un fichier .npy
        npy_file_name = "Doe_for_hydraulic.npy"
        #print(matrice_resultats)
        np.save(npy_file_name, matrice_resultats)
        #print(f"Matrice sauvegardée dans le fichier {npy_file_name}")


s = GenerateDoeHydrau().generateDoeHydrau()


             
