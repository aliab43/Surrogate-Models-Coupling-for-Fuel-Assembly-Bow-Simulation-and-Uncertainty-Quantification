import numpy as np
from pathlib import Path
from ClassHydraulicSurrogate import HydraulicSurrogateModel
from PenetrationDetection import PenetrationDetection
from GridClampingUpdate import GridClampingUpdater
import copy
import concurrent.futures

BASE_DIR = Path(__file__).resolve().parent

def load_doe(filename):
    candidates = [
        BASE_DIR / "DOE_for_mechanics" / filename,
        BASE_DIR / filename,
    ]
    for path in candidates:
        if path.exists():
            return np.load(path)
    raise FileNotFoundError(
        f"DOE file not found. Tried: {', '.join(str(p) for p in candidates)}"
    )




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





def GetFhCoeff(fh):
    """
    Projects a force field onto the provided modal basis to obtain modal coefficients.

    Parameters:
        fh : ndarray (assemblies x degrees of freedom)
            Force field for each assembly (rows: assemblies, columns: DOFs).

    Returns:
        modal_coefficients : ndarray (assemblies x modes)
            Calculated modal coefficients.
    """
    # Modal basis (provided)
    q = np.array([
        [0, -0.34182351, -0.33986812, -0.35552546, -0.36900858, -0.37084271, -0.36102097, -0.34508458, -0.3437104, 0],  # F_0
        [0, -0.55937995, -0.39424153, -0.06243252,  0.30706583,  0.46809913,  0.36218138,  0.07404396, -0.27875758, 0],  # F_C
        [0, -0.00479597,  0.31635535,  0.46610865,  0.36727189,  0.06233596, -0.26210788, -0.47490219, -0.49963184, 0],  # F_S
        [0,  0.19894744, -0.4367867,  -0.15343436,  0.43733209,  0.20649579, -0.39132532, -0.36346723,  0.47639416, 0]   # F_W
    ]).T  # Transposed to match the (DOFs x modes) format

    # Use the full force field and modal basis
    fh_full = fh
    q_full = q

    # Calculate modal coefficients
    num_modes = q_full.shape[1]
    modal_coefficients = np.zeros((fh_full.shape[0], num_modes))

    for i in range(fh_full.shape[0]):  # Loop over assemblies
        for j in range(num_modes):  # Loop over modes
            modal_coefficients[i, j] = np.dot(fh_full[i, :], q_full[:, j])

    return modal_coefficients

def compute_displacements_and_save(C_i_j_0, S_i_j_0, W_i_j_0, filename="", phase=""):
    """
    Computes displacements, saves them to a text file with appropriate phase prefixes,
    and also saves the modal coefficients (C, S, W).

    Parameters:
        C_i_j_0, S_i_j_0, W_i_j_0: Modal coefficients for displacement calculation.
        filename: Name of the output text file.
        phase: Prefix indicating the simulation phase ('b' for before, 'i' for during, 'a' for after).
    """
    assemblies = []
    c_basis, s_basis, w_basis = GramSchmidt()
    
    # Compute displacements
    for i in range(15):  # 15 fuel assemblies
        x = np.zeros(10)  # 10 grids per assembly
        for j in range(10):
            deformation_C = C_i_j_0[i] * c_basis[j]
            deformation_S = S_i_j_0[i] * s_basis[j]
            deformation_W = W_i_j_0[i] * w_basis[j]
            x[j] = deformation_C + deformation_S + deformation_W
        assemblies.append(x)
    
    # Save results to a text file
    with open(filename, "w") as f:
        # Save modal coefficients
        for i in range(15):
            f.write(f"C_{i+1}_{phase} = {C_i_j_0[i]:.6e}\n")
            f.write(f"S_{i+1}_{phase} = {S_i_j_0[i]:.6e}\n")
            f.write(f"W_{i+1}_{phase} = {W_i_j_0[i]:.6e}\n")

        f.write("\n")  # Separate modal coefficients from displacements

        # Save displacements
        for i in range(15):
            for j in range(10):
                f.write(f"U_G_{phase}_{j+1}_A_{i+1} = {assemblies[i][j]:.6e}\n")


def RunBeforeIrradiationCoupledSimulation(max_iterations=100, relaxation_factor= 0.2, tolerance=2.5e-3, step = 1):
    """
    Runs the coupled simulation between mechanical and hydraulic surrogate models
    with relaxation to ensure convergence.
    
    Returns:
        c, s, w : final converged modal coefficients
    """
    # Load initial mechanical DOE
    x_0 = load_doe("Doe_for_mechanics.npy")
    
    # Load hydraulic DOE
    X_hydraulic_0 = load_doe("Doe_for_hydraulic.npy")
    
    # Initialize hydraulic surrogate model
    hydraulic_model = HydraulicSurrogateModel()
    
    # Generate initial deformation
    print("=============================================\n    Fuel Assemblies insertion \n=============================================")
    print("               |")
    print("               v")
    #generate_deformation = PenetrationDetectionforTest(x_0,step).GenerateFiltredDeformationforTest()
    generate_deformation = PenetrationDetection(x_0,step).GenerateFiltredDeformation()
    c, s, w = generate_deformation
    # Store initial modal coefficients
    c_temp, s_temp, w_temp = c.copy(), s.copy(), w.copy()
    print("=============================================\n    Pump Activation \n=============================================")
    for iteration in range(max_iterations):
        
        print("               |")
        print("               v")
        print(f"--- FSI Coupling iteration {iteration + 1} ---")
        print("               |")
        print("               v")
        
        
        # Update hydraulic input and predict force field
        X_hydraulic = np.concatenate((c, s, w, X_hydraulic_0[0]))
      
        force_matrix = hydraulic_model.predict_force_field(X_hydraulic)
        #print(force_matrix)
        
        
        relaxed_force = force_matrix * relaxation_factor
        #print(force_matrix)
        # Compute modal coefficients for the force field
        fh_modal_coeff = GetFhCoeff(relaxed_force.T) 
        #print(fh_modal_coeff)
        # Update mechanical DOE with relaxed coefficients
        for i in range(15):
            x_0[i][3:7] += fh_modal_coeff[i]
        
        # Recompute filtered deformation
        generate_deformation = PenetrationDetection(x_0, step).GenerateFiltredDeformation()
        c, s, w = generate_deformation
        
        # Check convergence
        delta_c = np.abs(c - c_temp)
        delta_s = np.abs(s - s_temp)
        delta_w = np.abs(w - w_temp)
        
        if np.max(delta_c) < tolerance and np.max(delta_s) < tolerance and np.max(delta_w) < tolerance:
            print("Convergence achieved!")
            return c, s, w
        
        # Update previous modal coefficients for next iteration
        c_temp, s_temp, w_temp = c.copy(), s.copy(), w.copy()
    
    print("Maximum iterations reached without convergence.")
    return c , s , w 



def RunIrradiationPhase(initial_c, initial_s, initial_w, time=18, max_iterations=100, relaxation_factor=0.05, tolerance=2.5e-03):
    """
    Runs the irradiation phase simulation for 18 steps, including grid clamping updates.
    
    Args:
        initial_c, initial_s, initial_w : Modal coefficients from the pre-irradiation phase.
        time (int): Number of steps (e.g., 18 for 18 months).
        max_iterations (int): Max iterations per FSI loop for each step.
        relaxation_factor (float): Relaxation factor for modal coefficient updates.
        tolerance (float): Convergence tolerance for the FSI loop.
    
    Returns:
        final_c, final_s, final_w : Modal coefficients at the end of irradiation phase.
        final_grid_clamping : List of grid clamping values at the end of irradiation.
    """
    # Load initial mechanical DOE and hydraulic DOE
    x_0 = load_doe("Doe_for_mechanics.npy")
    X_hydraulic_0 = load_doe("Doe_for_hydraulic.npy")
    
    # Initialize hydraulic surrogate model
    hydraulic_model = HydraulicSurrogateModel()
    
    # Set initial modal coefficients
    c, s, w = initial_c, initial_s, initial_w
    final_grid_clamping = [0] * 15  # To store the final grid clamping values for each grid
    
    for month in range(1, time + 1):
        print(f"=============================================\n    Irradiation Step : {month} month(s)   \n=============================================")
        print("               |")
        print("               v")
        
        # Update mechanical DOE with current modal coefficients
        for i in range(15):  
            x_0[i][0] = c[i]
            x_0[i][1] = s[i]
            x_0[i][2] = w[i]
        
        for iteration in range(max_iterations):
            print(f"--- FSI Coupling iteration {iteration + 1} ---")
            
            # Update hydraulic input and predict force field
            X_hydraulic = np.concatenate((c, s, w, X_hydraulic_0[0]))
            force_matrix = hydraulic_model.predict_force_field(X_hydraulic)
            #print(force_matrix)

        
            relaxed_force = force_matrix * relaxation_factor
            fh_modal_coeff = GetFhCoeff(relaxed_force.T)

            
            # Update mechanical DOE with relaxed coefficients
            for i in range(15):
                x_0[i][3:7] += fh_modal_coeff[i]
            
            # Update grid clamping values
            for i in range(15):
                initial_grid_clamping = x_0[i][10]  # Get current grid clamping value
                updater = GridClampingUpdater(initial_grid_clamping)  # Initialize updater
                updated_grid_clamping = updater.update_Grid_Clamping()  # Update grid clamping
                
                # Update x_0 with the new grid clamping value
                x_0[i][10] = updated_grid_clamping
                final_grid_clamping[i] = updated_grid_clamping  # Save the latest clamping value
            
            # Recompute filtered deformation
            generate_deformation = PenetrationDetection(x_0, step=2).GenerateFiltredDeformation()
            c_new, s_new, w_new = generate_deformation
            
            # Check convergence
            delta_c = np.abs(c_new - c)
            delta_s = np.abs(s_new - s)
            delta_w = np.abs(w_new - w)
            
            if np.max(delta_c) < tolerance and np.max(delta_s) < tolerance and np.max(delta_w) < tolerance:
                print("Coupling convergence achieved!")
                c, s, w = c_new, s_new, w_new
                break  # Exit the loop for this step
            
            # Update for next iteration
            c, s, w = c_new, s_new, w_new
        else:
            print(f"Step {month} did not converge within {max_iterations} iterations.")
    
    return c, s, w, final_grid_clamping




def RunAfterIrradiationPhase(irradiation_c, irradiation_s, irradiation_w, final_grid_clamping, max_iterations = 1, relaxation_factor=0.1, tolerance=2.5e-03):
    """
    Runs the post-irradiation phase simulation.

    Args:
        irradiation_c, irradiation_s, irradiation_w : Modal coefficients from the irradiation phase.
        final_grid_clamping : List of grid clamping values from the end of irradiation.
        max_iterations (int): Max iterations for FSI loop.
        relaxation_factor (float): Relaxation factor for modal coefficient updates.
        tolerance (float): Convergence tolerance for the FSI loop.

    Returns:
        final_c, final_s, final_w : Modal coefficients at the end of the post-irradiation phase.
    """
    # Load mechanical DOE and hydraulic DOE
    x_0 = load_doe("Doe_for_mechanics.npy")
    X_hydraulic_0 = load_doe("Doe_for_hydraulic.npy")
    
    # Initialize hydraulic surrogate model
    hydraulic_model = HydraulicSurrogateModel()
    
    # Set initial modal coefficients
    c, s, w = irradiation_c, irradiation_s, irradiation_w
    
    print(f"=============================================\n   Reactor Power Off: Post-Irradiation Phase   \n=============================================")
    print("               |")
    print("               v")
    
    # Update mechanical DOE with current modal coefficients and grid clamping values
    for i in range(15):
        x_0[i][0] = c[i]
        x_0[i][1] = s[i]
        x_0[i][2] = w[i]
        x_0[i][10] = final_grid_clamping[i]  # Use the final grid clamping from irradiation phase
    c_temp, s_temp, w_temp = c.copy(), s.copy(), w.copy()
    for iteration in range(max_iterations):
        
        print(f"--- FSI Coupling iteration {iteration + 1} ---")
        print("               |")
        print("               v")
        
        
        # Update hydraulic input and predict force field
        X_hydraulic = np.concatenate((c, s, w, X_hydraulic_0[0]))
        force_matrix = hydraulic_model.predict_force_field(X_hydraulic)
        relaxed_force = force_matrix * relaxation_factor
        fh_modal_coeff = GetFhCoeff(relaxed_force.T)
        # Update mechanical DOE with relaxed coefficients
        for i in range(15):
            x_0[i][3:7] += fh_modal_coeff[i]
        
        
        # Recompute filtered deformation
        generate_deformation = PenetrationDetection(x_0, step = 3).GenerateFiltredDeformation()
        c, s, w = generate_deformation
        
        # Check convergence
        delta_c = np.abs(c - c_temp)
        delta_s = np.abs(s - s_temp)
        delta_w = np.abs(w - w_temp)
        
        if np.max(delta_c) < tolerance and np.max(delta_s) < tolerance and np.max(delta_w) < tolerance:
            print("Convergence achieved!")
            return c, s, w
        
        # Update previous modal coefficients for next iteration
        c_temp, s_temp, w_temp = c.copy(), s.copy(), w.copy()
    
    print("End of cycle")
    return c , s , w 



# Step 1: Run Pre-Irradiation Phase
initial_c, initial_s, initial_w = RunBeforeIrradiationCoupledSimulation(max_iterations=2, relaxation_factor=0.1, tolerance=2.5e-3, step=1)


compute_displacements_and_save(np.array(initial_c), np.array(initial_s), np.array(initial_w), filename="displacement_before_irradiation.txt", phase="b")

# Step 2: Run Irradiation Phase
final_c, final_s, final_w, final_grid_clamping = RunIrradiationPhase(initial_c, initial_s, initial_w, time=18, max_iterations=100, relaxation_factor=0.1, tolerance=2.5e-03)

compute_displacements_and_save(np.array(final_c), np.array(final_s), np.array(final_w), filename="displacement_irradiation.txt", phase="i")


# Step 3: Run Post-Irradiation Phase
post_irradiation_c, post_irradiation_s, post_irradiation_w = RunAfterIrradiationPhase(final_c, final_s, final_w, final_grid_clamping, max_iterations=2, relaxation_factor=0.2, tolerance=2.5e-03)

compute_displacements_and_save(np.array(post_irradiation_c), np.array(post_irradiation_s), np.array(post_irradiation_w), filename="displacement_after_irradiation.txt", phase="a")


def run_central_cycle_and_collect_inputs():
    """
    Exécute le cycle central (moyenne des PG) et collecte tous les points d'entrée
    utilisés pour chaque métamodèle lors des boucles de couplage.
    Retourne :
        - central_results: résultats centraux
        - inputs_mechanical: liste des inputs pour chaque appel au métamodèle mécanique
        - inputs_hydraulic: liste des inputs pour chaque appel au métamodèle hydraulique
    """
    # Initialisation des listes pour stocker les inputs
    inputs_mechanical = []
    inputs_hydraulic = []

    # Fonctions wrappers pour collecter les inputs
    from ClassMechanicalSurrogates import MechanicalSurrogate
    from ClassHydraulicSurrogate import HydraulicSurrogateModel

    class MechanicalSurrogateWithLog(MechanicalSurrogate):
        def callSurrogateModelsday7(self, xnew):
            inputs_mechanical.append(('day7', np.copy(xnew)))
            return super().callSurrogateModelsday7(xnew)
        def callSurrogateModelCreep(self, xnew):
            inputs_mechanical.append(('creep', np.copy(xnew)))
            return super().callSurrogateModelCreep(xnew)
        def callSurrogateModelday41(self, xnew):
            inputs_mechanical.append(('day41', np.copy(xnew)))
            return super().callSurrogateModelday41(xnew)

    class HydraulicSurrogateModelWithLog(HydraulicSurrogateModel):
        def predict_force_field(self, x):
            inputs_hydraulic.append(np.copy(x))
            return super().predict_force_field(x)

    # Patch les modules utilisés dans PenetrationDetection et le reste
    import PenetrationDetection
    PenetrationDetection.MechanicalSurrogate = MechanicalSurrogateWithLog
    import ClassHydraulicSurrogate
    ClassHydraulicSurrogate.HydraulicSurrogateModel = HydraulicSurrogateModelWithLog

    # Exécuter le cycle complet comme d'habitude
    initial_c, initial_s, initial_w = RunBeforeIrradiationCoupledSimulation(max_iterations=2, relaxation_factor=0.1, tolerance=2.5e-3, step=1)
    compute_displacements_and_save(np.array(initial_c), np.array(initial_s), np.array(initial_w), filename="displacement_before_irradiation.txt", phase="b")
    final_c, final_s, final_w, final_grid_clamping = RunIrradiationPhase(initial_c, initial_s, initial_w, time=18, max_iterations=100, relaxation_factor=0.1, tolerance=2.5e-03)
    compute_displacements_and_save(np.array(final_c), np.array(final_s), np.array(final_w), filename="displacement_irradiation.txt", phase="i")
    post_irradiation_c, post_irradiation_s, post_irradiation_w = RunAfterIrradiationPhase(final_c, final_s, final_w, final_grid_clamping, max_iterations=2, relaxation_factor=0.2, tolerance=2.5e-03)
    compute_displacements_and_save(np.array(post_irradiation_c), np.array(post_irradiation_s), np.array(post_irradiation_w), filename="displacement_after_irradiation.txt", phase="a")

    central_results = {
        "before": (initial_c, initial_s, initial_w),
        "irradiation": (final_c, final_s, final_w),
        "after": (post_irradiation_c, post_irradiation_s, post_irradiation_w)
    }
    return central_results, inputs_mechanical, inputs_hydraulic

def generate_pg_realization(mean, cov):
    """
    Génère une réalisation d'un processus gaussien multivarié corrélé.
    """
    return np.random.multivariate_normal(mean, cov)

def monte_carlo_cycle(n_mc=1000):
    """
    Boucle Monte Carlo pour la quantification d'incertitude (parallélisée).
    """
    central_results, inputs_mechanical, inputs_hydraulic = run_central_cycle_and_collect_inputs()

    # Préparer les arguments pour chaque réalisation
    args = (inputs_mechanical, inputs_hydraulic)
    args_list = [args for _ in range(n_mc)]

    all_solutions = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=20) as executor:
        for sol in executor.map(run_one_mc_realization, args_list):
            all_solutions.append(sol)

    # 4. Analyse statistique
    
    after_c = np.array([sol["after"][0] for sol in all_solutions])
    after_s = np.array([sol["after"][1] for sol in all_solutions])
    after_w = np.array([sol["after"][2] for sol in all_solutions])
    mean_c = np.mean(after_c, axis=0)
    mean_s = np.mean(after_s, axis=0)
    mean_w = np.mean(after_w, axis=0)
    std_c = np.std(after_c, axis=0)
    print("Empirical mean (C):", mean_c)
    print("Empirical std (C):", std_c)
    print("Empirical mean (S):", mean_s)
    print("Empirical std (S):", np.std(after_s, axis=0))        
    print("Empirical mean (W):", mean_w)
    print("Empirical std (W):", np.std(after_w, axis=0))
    # Pour les quantiles :
    q05_c = np.quantile(after_c, 0.05, axis=0)
    q95_c = np.quantile(after_c, 0.95, axis=0)
    q95_s = np.quantile(after_s, 0.95, axis=0)
    q95_w = np.quantile(after_w, 0.95, axis=0)
    q_05_s = np.quantile(after_s, 0.05, axis=0)
    q_05_w = np.quantile(after_w, 0.05, axis=0)
           
    print("5% quantile (C):", q05_c)
    print("95% quantile (C):", q95_c)
    print("5% quantile (S):", q_05_s)
    print("95% quantile (S):", q95_s)   
    print("5% quantile (W):", q_05_w)
    print("95% quantile (W):", q95_w)


def run_one_mc_realization(args):
    # args = (inputs_mechanical, inputs_hydraulic)
    inputs_mechanical, inputs_hydraulic = args

    from ClassMechanicalSurrogates import MechanicalSurrogate
    from ClassHydraulicSurrogate import HydraulicSurrogateModel

    meca = MechanicalSurrogate()
    meca_realizations = {}
    for mode, x in inputs_mechanical:
        # Générer une perturbation pour chaque sortie (c, s, w)
        if mode == 'day7':
            c_mean, s_mean, w_mean = meca.kriging_modelC7.kriging_function(x), meca.kriging_modelS7.kriging_function(x), meca.kriging_modelW7.kriging_function(x)
            c_std, s_std, w_std = np.sqrt(2.5e-07), np.sqrt(2.77779236530208e-09), np.sqrt(1.61e-08)
        elif mode == 'creep':
            c_mean, s_mean, w_mean = meca.kriging_modelVC.kriging_function(x), meca.kriging_modelVS.kriging_function(x), meca.kriging_modelVW.kriging_function(x)
            c_std, s_std, w_std = np.sqrt(5.8e-13), np.sqrt(4.8e-13), np.sqrt(4.7e-13)
        elif mode == 'day41':
            c_mean, s_mean, w_mean = meca.kriging_modelC41.kriging_function(x), meca.kriging_modelS41.kriging_function(x), meca.kriging_modelW41.kriging_function(x)
            c_std, s_std, w_std = np.sqrt(3.9e-07), np.sqrt(2.6e-09), np.sqrt(1.1e-08)
        else:
            continue
        # Générer une réalisation gaussienne pour chaque sortie
        c_real = c_mean + c_std * np.random.randn(*np.atleast_1d(c_mean).shape)
        s_real = s_mean + s_std * np.random.randn(*np.atleast_1d(s_mean).shape)
        w_real = w_mean + w_std * np.random.randn(*np.atleast_1d(w_mean).shape)
        c_diff = c_real - c_mean
        s_diff = s_real - s_mean
        w_diff = w_real - w_mean
        meca_realizations[(mode, tuple(x))] = (c_diff, s_diff, w_diff)

    # Pour chaque métamodèle hydraulique, générer une réalisation
    hydrau = HydraulicSurrogateModel()
    hydrau_realizations = {}
    for x in inputs_hydraulic:
        # Idem : moyenne et std pour modèle hydraulique)
        mean = np.array(hydrau.predict_force_field(x))
        std = 15 # Newton écat_type de prédiction de force   
        realization = mean + std * np.random.randn(*mean.shape)
        diff = realization - mean
        hydrau_realizations[tuple(x)] = diff

    # Patch pour cette réalisation
    class MechanicalSurrogateWithPerturb(MechanicalSurrogate):
        def callSurrogateModelsday7(self, xnew):
            key = ('day7', tuple(xnew))
            c_delta, s_delta, w_delta = meca_realizations.get(key, (0, 0, 0))
            c, s, w = super().callSurrogateModelsday7(xnew)
            return np.array(c) + c_delta, np.array(s) + s_delta, np.array(w) + w_delta
        def callSurrogateModelCreep(self, xnew):
            key = ('creep', tuple(xnew))
            c_delta, s_delta, w_delta = meca_realizations.get(key, (0, 0, 0))
            c, s, w = super().callSurrogateModelCreep(xnew)
            return np.array(c) + c_delta, np.array(s) + s_delta, np.array(w) + w_delta
        def callSurrogateModelday41(self, xnew):
            key = ('day41', tuple(xnew))
            c_delta, s_delta, w_delta = meca_realizations.get(key, (0, 0, 0))
            c, s, w = super().callSurrogateModelday41(xnew)
            return np.array(c) + c_delta, np.array(s) + s_delta, np.array(w) + w_delta

    class HydraulicSurrogateModelWithPerturb(HydraulicSurrogateModel):
        def predict_force_field(self, x):
            delta = hydrau_realizations.get(tuple(x), 0)
            return super().predict_force_field(x) + delta

    import PenetrationDetection
    PenetrationDetection.MechanicalSurrogate = MechanicalSurrogateWithPerturb
    import ClassHydraulicSurrogate
    ClassHydraulicSurrogate.HydraulicSurrogateModel = HydraulicSurrogateModelWithPerturb

    # Exécuter le cycle complet
    initial_c, initial_s, initial_w = RunBeforeIrradiationCoupledSimulation(max_iterations=2, relaxation_factor=0.1, tolerance=2.5e-3, step=1)
    final_c, final_s, final_w, final_grid_clamping = RunIrradiationPhase(initial_c, initial_s, initial_w, time=18, max_iterations=100, relaxation_factor=0.1, tolerance=2.5e-03)
    post_irradiation_c, post_irradiation_s, post_irradiation_w = RunAfterIrradiationPhase(final_c, final_s, final_w, final_grid_clamping, max_iterations=2, relaxation_factor=0.2, tolerance=2.5e-03)

    return {
        "before": (initial_c, initial_s, initial_w),
        "irradiation": (final_c, final_s, final_w),
        "after": (post_irradiation_c, post_irradiation_s, post_irradiation_w)
    }






monte_carlo_cycle(n_mc=100)





