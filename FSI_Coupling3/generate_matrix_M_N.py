import numpy as np


nb_assembly = 15

def generate_one_sample():
    
    mean_C= [-1.80640018e-04,  3.41280644e-04, -2.29175954e-03, -9.60219314e-05,
            -8.49123937e-05,  4.15388725e-03,  3.26323131e-03,  3.59113552e-04,
            -1.74406149e-03, -1.54797464e-03,  1.04884639e-05,  6.70227483e-05,
            -7.77740556e-05,  9.72333656e-04,  8.50323877e-04]
    mean_S = [8.00304835e-04,  1.07567214e-04, -1.44638224e-03, -1.68008800e-05,
        -1.70367470e-05, -2.67415013e-03, -1.01068999e-04,  2.60543181e-04,
        -1.35827353e-03,  3.66635031e-03,  6.36241473e-05,  3.47234554e-06,
        1.15619573e-03, -1.39904720e-03, -7.93415783e-05]
    mean_W = [9.41636894e-05,  3.19043207e-04, -3.31703156e-04,  2.23100689e-05,
            2.42452488e-05,  7.73366890e-04,  3.42280097e-04, -1.34272389e-03,
            -7.74612902e-04,  1.24295279e-04,  4.27688330e-05, -9.52701364e-06,
            -9.23258449e-04, -4.28579229e-04,  1.12371381e-04]
    
    Sigma_mesurmentsC = 0.0003
    Sigma_mesurementsS = 0.0003
    sigma_mesurmentsW = 0.0003
    cov_diagC = Sigma_mesurmentsC ** 2
    cov_diagS = Sigma_mesurementsS ** 2
    cov_diagW = sigma_mesurmentsW ** 2
    SigmaC = np.diag([cov_diagC] * 15)
    SigmaS = np.diag([cov_diagS] * 15)
    SigmaW = np.diag([cov_diagW] * 15)

    # Generate multivariate samples with specified mean and diagonal covariance
    C_i_j_C = np.random.multivariate_normal(mean=mean_C, cov=SigmaC).tolist()
    C_i_j_S = np.random.multivariate_normal(mean=mean_S, cov=SigmaS).tolist()
    C_i_j_W = np.random.multivariate_normal(mean=mean_W, cov=SigmaW).tolist()

    v_m = np.random.normal(5, 0.0005)  # the mean velocity ~ Normal distribution
    max_deviation_in = np.random.normal(0.05, 0.00005) # ~ Normal distribution
    l_offset_in = np.random.normal(0.0, 0.01)  # ~ Normal distribution
    max_deviation_out = np.random.normal(0.04, 0.00004)  # ~ Normal distributionn
    l_offset_out = np.random.normal(0.0, 0.005)  # ~ Normal distribution
    hl_sensi =  np.random.uniform(0.01,0.03) # ~ Uniform distribution
    pre_charge_axial = np.random.uniform(0.01, 0.025)  # ~ Uniform distribution
    grid_clamping = np.random.normal(37, 0.5)

    constant_growth = np.random.normal(loc = 1, scale = 0.3)
    constant_creep = np.random.normal(loc = 1, scale = 0.3)
    Hydraulic_BC = [v_m, max_deviation_in, l_offset_in, max_deviation_out, l_offset_out]
    return [C_i_j_C, C_i_j_S, C_i_j_W, Hydraulic_BC, hl_sensi, pre_charge_axial,  
            grid_clamping, constant_growth, constant_creep]
    

Ns = 200
file_name =  "Matrix_M_N.dat"

mean_C= [-1.80640018e-04,  3.41280644e-04, -2.29175954e-03, -9.60219314e-05,
            -8.49123937e-05,  4.15388725e-03,  3.26323131e-03,  3.59113552e-04,
            -1.74406149e-03, -1.54797464e-03,  1.04884639e-05,  6.70227483e-05,
            -7.77740556e-05,  9.72333656e-04,  8.50323877e-04]
mean_S = [8.00304835e-04,  1.07567214e-04, -1.44638224e-03, -1.68008800e-05,
        -1.70367470e-05, -2.67415013e-03, -1.01068999e-04,  2.60543181e-04,
        -1.35827353e-03,  3.66635031e-03,  6.36241473e-05,  3.47234554e-06,
        1.15619573e-03, -1.39904720e-03, -7.93415783e-05]
mean_W = [9.41636894e-05,  3.19043207e-04, -3.31703156e-04,  2.23100689e-05,
            2.42452488e-05,  7.73366890e-04,  3.42280097e-04, -1.34272389e-03,
            -7.74612902e-04,  1.24295279e-04,  4.27688330e-05, -9.52701364e-06,
            -9.23258449e-04, -4.28579229e-04,  1.12371381e-04]


with open(file_name, 'w') as file : 
    file.write("#NAME: Matrix_M_N\n")
    file.write("#COLUMN_NAMES: Coeff_C | Coeff_S | Coeff_W | Hydraulic_BC | hl_sensi | pre_charge_axial | grid_clamping | constant_creep | constant_growth | Matrix_M_N__n__iter__ \n")
    file.write("#COLUMN_TYPES: V | V | V | V | D | D | D | D | D | D  \n")
    file.write("\n")

    txtFormat1 = "{price:.8f},"
    txtFormat2 = "{price:.8f}] "
    txtFloat = "{} "
    for i in range(Ns):
        if i == 0:
            # First sample: use mean values only (no variance)
            vect_one = [
                mean_C.copy(),
                mean_S.copy(),
                mean_W.copy(),
                [5.0, 0.05, 0.0, 0.04, 0.0],  # mean values for Hydraulic_BC
                0.025,                        # mean for hl_sensi (or pick your preferred mean)
                0.0175,                       # mean for pre_charge_axial
                37.0,                         # mean for grid_clamping
                1.0,                          # mean for constant_growth
                1.0                           # mean for constant_creep
            ]
        else:
            vect_one = generate_one_sample()
        file.write("[")
        for j in range(nb_assembly-1):
            file.write(txtFormat1.format(price = vect_one[0][j]))
        file.write(txtFormat2.format(price = vect_one[0][nb_assembly-1]))
        file.write("[")
        for j in range(nb_assembly-1):
            file.write(txtFormat1.format(price = vect_one[1][j]))
        file.write(txtFormat2.format(price = vect_one[1][nb_assembly-1]))
        file.write("[")
        for j in range(nb_assembly-1):
            file.write(txtFormat1.format(price = vect_one[2][j]))
        file.write(txtFormat2.format(price = vect_one[2][nb_assembly-1]))
        file.write("[")
        for j in range(4):
            file.write(txtFormat1.format(price = vect_one[3][j]))
        file.write(txtFormat2.format(price = vect_one[3][4]))

        file.write(txtFloat.format(vect_one[4]))
        file.write(txtFloat.format(vect_one[5]))
        file.write(txtFloat.format(vect_one[6]))
        file.write(txtFloat.format(vect_one[7]))
        file.write(txtFloat.format(vect_one[8]))
        file.write(txtFloat.format(1.0*i))
        file.write("\n")




