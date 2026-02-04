import numpy as np
import matplotlib.pyplot as plt

v_m = np.random.normal(5, 0.005)  # the mean velocity ~ Normal distribution
max_deviation_in = np.random.normal(0.05, 0.00005) # ~ Normal distribution
l_offset_in = np.random.normal(0.0, 0.01)  # ~ Normal distribution
max_deviation_out = np.random.normal(0.04, 0.00004)  # ~ Normal distributionn
l_offset_out = np.random.normal(0.0, 0.005)  # ~ Normal distribution

def parabolic_profile(mean_axial_velocity, maximal_deviation, lateral_offset):
    x_min = 0.
    x_max = 1.
    denom = (x_max**2 - 2*x_max*x_min + x_min**2 + 12*lateral_offset**2)
    a = -12 * maximal_deviation * mean_axial_velocity / denom
    b = (12 * (x_max + x_min + 2 * lateral_offset)) * maximal_deviation * mean_axial_velocity / denom
    c = -2 * mean_axial_velocity * (
        (maximal_deviation - 1/2) * x_max**2 +
        ((4*maximal_deviation + 1)*x_min + 6*lateral_offset*maximal_deviation)*x_max +
        (maximal_deviation - 1/2)*x_min**2 +
        6*x_min*lateral_offset*maximal_deviation - 6*lateral_offset**2
    ) / denom
    def p(x):
        return a * x**2 + b * x + c
    return p

# Reference (mean) parameters for injection (entrant)
mean_v_m_in = 5.0
mean_max_deviation_in = 0.05
mean_l_offset_in = 0.0
std_v_m_in = 0.005
std_max_deviation_in = 0.0005
std_l_offset_in = 0.02

# Reference (mean) parameters for withdrawal (sortant)
mean_v_m_out = 5.0
mean_max_deviation_out = -0.04
mean_l_offset_out = 0.0
std_v_m_out = 0.0005
std_max_deviation_out = 0.0004
std_l_offset_out = 0.01

N = 2000
x_vals = np.linspace(0, 1, 2000)
profiles_in = []
profiles_out = []

for _ in range(N):
    # Injection (entrant)
    v_m_sample_in = np.random.normal(mean_v_m_in, std_v_m_in)
    max_dev_sample_in = np.random.normal(mean_max_deviation_in, std_max_deviation_in)
    l_offset_sample_in = np.random.normal(mean_l_offset_in, std_l_offset_in)
    p_in = parabolic_profile(v_m_sample_in, max_dev_sample_in, l_offset_sample_in)
    profiles_in.append(p_in(x_vals))
    # Withdrawal (sortant)
    v_m_sample_out = np.random.normal(mean_v_m_out, std_v_m_out)
    max_dev_sample_out = np.random.normal(mean_max_deviation_out, std_max_deviation_out)
    l_offset_sample_out = np.random.normal(mean_l_offset_out, std_l_offset_out)
    p_out = parabolic_profile(v_m_sample_out, max_dev_sample_out, l_offset_sample_out)
    profiles_out.append(p_out(x_vals))

profiles_in = np.array(profiles_in)
profiles_out = np.array(profiles_out)

mean_profile_in = np.mean(profiles_in, axis=0)
std_profile_in = np.std(profiles_in, axis=0)
mean_profile_out = np.mean(profiles_out, axis=0)
std_profile_out = np.std(profiles_out, axis=0)

plt.figure(figsize=(10, 6))
# Plot a few random profiles for illustration
for i in range(0, N, max(1, N // 10)):
    plt.plot(x_vals, profiles_in[i], color='green', alpha=0.2, linewidth=1)
    plt.plot(x_vals, profiles_out[i], color='brown', alpha=0.2, linewidth=1)
# Plot mean profiles
plt.plot(x_vals, mean_profile_in, color='green', label='Mean injection profile', linewidth=2)
plt.plot(x_vals, mean_profile_out, color='brown', label='Mean withdrawal profile', linewidth=2)
# Plot confidence intervals
plt.fill_between(x_vals, mean_profile_in - 2*std_profile_in, mean_profile_in + 2*std_profile_in, color='lime', alpha=0.3, label='Injection 95% CI')
plt.fill_between(x_vals, mean_profile_out - 2*std_profile_out, mean_profile_out + 2*std_profile_out, color='orange', alpha=0.3, label='Withdrawal 95% CI')
plt.xlabel('Adimensional length', fontsize=14)
plt.ylabel('Speed (m/s)', fontsize=14)
plt.title('Velocity profiles with uncertainty (injection and withdrawal)', fontsize=16)
plt.legend(fontsize=12)
plt.grid(True)
plt.tight_layout()
plt.savefig('velocity_profiles_injection_withdrawal_uncertainty.png')
plt.show()

