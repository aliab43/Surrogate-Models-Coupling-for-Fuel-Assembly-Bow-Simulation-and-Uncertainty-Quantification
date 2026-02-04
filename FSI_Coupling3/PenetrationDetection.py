import numpy as np
from pathlib import Path
from ClassMechanicalSurrogates import MechanicalSurrogate
import matplotlib.pyplot as plt
import time

BASE_DIR = Path(__file__).resolve().parent



class PenetrationDetection:

    def __init__(self,X, step):
        
        self.X = X
        self.step = step
        self.surrogate_model = MechanicalSurrogate()
        self.nominal_hot_watergap = 0.002

        

    # Function to plot assembly deformations
    def plot_assembly_deformation(self,assemblies):
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
            # Do not reverse the displacement vector
            x_plot =  x[::-1]
            plt.xlim(self.nominal_hot_watergap - self.nominal_hot_watergap * 1.0, 15 * self.nominal_hot_watergap + self.nominal_hot_watergap * 1.0)
            plt.ylim(0., 11.)
            # Set x-ticks to match nominal water gap positions
            xtick_positions = [self.nominal_hot_watergap + self.nominal_hot_watergap * k for k in range(15)]
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
        
        output_dir = BASE_DIR / "iteration_deformation_init"
        output_dir.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_dir / "deformation.png")
        plt.close('all')



    def getCoefAfterSimulation(self, X_news):
        # Initial deformation values
        c_0, s_0, w_0 = [], [], []
        for s in range(15):
            c0 = X_news[s][0]
            s0 = X_news[s][1]
            w0 = X_news[s][2]
            c_0.append(c0)
            s_0.append(s0)
            w_0.append(w0) 
        
        # Deformation from surrogate model
        C_i_C, C_i_S, C_i_W =  [], [], []
        for i in range(15):
            if self.step == 1 :
                c, s, w = self.surrogate_model.callSurrogateModelsday7(X_news[i])           
            elif self.step == 2:
                c, s, w = self.surrogate_model.callSurrogateModelCreep(X_news[i])
            else :
                c,s,w = self.surrogate_model.callSurrogateModelday41(X_news[i])
            # Forcer le scalaire même si c, s, w sont des tableaux ou listes
            C_i_C.append(float(np.atleast_1d(c).item()))
            C_i_S.append(float(np.atleast_1d(s).item()))
            C_i_W.append(float(np.atleast_1d(w).item()))
                
        # Add initial deformation to the new deformation element-wise
        C_i_f = [c_0[i] + C_i_C[i] for i in range(15)]
        S_i_f = [s_0[i] + C_i_S[i] for i in range(15)]
        W_i_f = [w_0[i] + C_i_W[i] for i in range(15)]
        
        return C_i_f, S_i_f, W_i_f

    

    def GramSchmidt(self):
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
    

    def GenerateFiltredDeformation(self): 
        # Get initial coefficients
        C_i_0, S_i_0, W_i_0 = self.getCoefAfterSimulation(self.X)
        C_i_j_0 = np.array(C_i_0)
        S_i_j_0 = np.array(S_i_0)
        W_i_j_0 = np.array(W_i_0)

        c, s, w = self.GramSchmidt()

        # Convergence parameters
        converged = False
        max_iterations = 100
        iteration = 0
        lambda_value = 2e+5
        damping_factor = -0.4 if self.step == 3 else 1  # Damping for step 3
        # Start timing
        start_time = time.time()
        while not converged and iteration < max_iterations:
            assemblies = []

            # Calculate deformation for each assembly
            for i in range(15):
                x = np.zeros(10)
                for j in range(10):
                    x[j] = self.nominal_hot_watergap + (i) * self.nominal_hot_watergap
                    deformation_C = C_i_j_0[i] * c[j]
                    deformation_S = S_i_j_0[i] * s[j]
                    deformation_W = W_i_j_0[i] * w[j]
                    x[j] += (deformation_C + deformation_S + deformation_W)
                assemblies.append(x)

            # Initialize penalty forces matrix
            penalty_forces = np.zeros((15, 10))

            # Reset convergence flag
            convergence_flag = True

            # Check for penetration and apply penalty forces
            for i in range(len(assemblies)):
                for j in range(10):
                    # Check penetration between neighboring assemblies
                    if i < len(assemblies) - 1:
                        penetration_depth = assemblies[i + 1][j] - assemblies[i][j]
                        if penetration_depth < 0:
                            penalty_force_magnitude = 500 * (1 - np.exp(-lambda_value * np.abs(penetration_depth))) / (1 + 10 * np.abs(penetration_depth))
                            penalty_forces[i][j] -= penalty_force_magnitude
                            penalty_forces[i + 1][j] += penalty_force_magnitude
                            convergence_flag = False

                        # Check penetration with boundaries (reactor vessel walls)


                    # Modify penalty logic for boundaries
                    if assemblies[0][j] < 0.0005:  # Boundary for assembly 1
                        penetration_depth =  assemblies[0][j] - 0.0005
                        penalty_force_magnitude = 500 * (1 - np.exp(-lambda_value * np.abs(penetration_depth))) / (1 + 10 * np.abs(penetration_depth))
                        penalty_forces[0][j] += penalty_force_magnitude*0.05
                        convergence_flag = False
                    if assemblies[14][j] > 0.0315 :  # Boundary for assembly 15
                        penetration_depth = assemblies[14][j] - 0.0315
                        penalty_force_magnitude = 500 * (1 - np.exp(-lambda_value * np.abs(penetration_depth))) / (1 + 10 * np.abs(penetration_depth))
                        penalty_forces[14][j] -= penalty_force_magnitude*0.05 
                        convergence_flag = False


            # Modal basis (10x4, 10 rows for degrees of freedom, 4 columns for modes)
            q = np.array([
                [0, -0.34182351, -0.33986812, -0.35552546, -0.36900858, -0.37084271, -0.36102097, -0.34508458, -0.3437104, 0],  # F_0
                [0, -0.55937995, -0.39424153, -0.06243252, 0.30706583, 0.46809913, 0.36218138, 0.07404396, -0.27875758, 0],   # F_C
                [0, -0.00479597, 0.31635535, 0.46610865, 0.36727189, 0.06233596, -0.26210788, -0.47490219, -0.49963184, 0],   # F_S
                [0, 0.19894744, -0.4367867, -0.15343436, 0.43733209, 0.20649579, -0.39132532, -0.36346723, 0.47639416, 0]     # F_W
            ]).T 

            # Number of modes
            num_modes = q.shape[1]

            # Initialize modal coefficients
            modal_coefficients = np.zeros((penalty_forces.shape[0], num_modes))
            #print(penalty_forces)
            # Compute modal coefficients for each assembly
            for i in range(penalty_forces.shape[0]):  # Loop over assemblies
                for j in range(num_modes):  # Loop over modes
                   
                    dot_product = np.dot(penalty_forces[i, :], q[:, j])  # Dot product
                    modal_coefficients[i, j] = dot_product / 64
            # Update X with new modal coefficients
            #print(modal_coefficients)
            for i in range(15):
                self.X[i][3:7] += damping_factor * modal_coefficients[i]

            # Get new coefficients for the next iteration
            C_i_1, S_i_1, W_i_1 = self.getCoefAfterSimulation(self.X)
            C_i_j_0 = np.array(C_i_1)
            S_i_j_0 = np.array(S_i_1)
            W_i_j_0 = np.array(W_i_1)


            if convergence_flag :
                converged = True
                break

            iteration += 1
            print(f"Mechanical loop Iteration {iteration}")
            #self.plot_assembly_deformation(assemblies)
        # End timing
        end_time = time.time()

        # Print convergence time
        convergence_time = end_time - start_time
        print(f"Mechanical loop converged in {convergence_time:.2f} seconds" if converged else f"Max iterations reached in {convergence_time:.2f} seconds")
        # Final deformation plot and convergence status
        #self.plot_assembly_deformation(assemblies)
        print("Mechanical loop Converged" if converged else "Max iterations reached without full convergence")
        return C_i_j_0, S_i_j_0, W_i_j_0














