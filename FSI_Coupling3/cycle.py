

import numpy as np
import matplotlib.pyplot as plt

# Adjusted time and displacement values for the phases
time_before_irradiation = [0, 7]  # Time points before irradiation
displacement_before_irradiation = [0.002, 0.003]  # Non-zero starting displacement for clarity

# Time for irradiation phase (18 months, 30-day steps)
num_steps = 18
step_duration = 30
time_irradiation = [7 + i * step_duration for i in range(num_steps)]

# Nonlinear trajectory: Overlaying a sine wave on a linear increase
base_displacement = 0.003  # Starting displacement
linear_increment_per_step = 0.00005  # Linear increment for reference
amplitude = 0.0002  # Amplitude of sine wave (controls non-linearity)
frequency = 2 * np.pi / len(time_irradiation)  # Frequency of sine wave

displacement_irradiation = [
    base_displacement + linear_increment_per_step * i + amplitude * np.sin(frequency * i)
    for i in range(len(time_irradiation))
]

# Time and displacement for post-irradiation phase
time_after_irradiation = [time_irradiation[-1], time_irradiation[-1] + 30]
displacement_after_irradiation = [displacement_irradiation[-1], 0.0028]  # Smooth transition to post-irradiation

# Combine all phases
time_combined = time_before_irradiation + time_irradiation + time_after_irradiation
displacement_combined = displacement_before_irradiation + displacement_irradiation + displacement_after_irradiation

# Plotting
plt.figure(figsize=(16, 10))
plt.plot(time_combined, displacement_combined, color='b', marker='o', linestyle='-', linewidth=2)

# Adding annotations
plt.axvline(x=time_before_irradiation[1], color='crimson', linestyle='--', linewidth=3, label='Power on')
plt.axvline(x=time_irradiation[-1], color='black', linestyle='--', linewidth=3, label='Power off')
plt.text(35, 0.0026, "Before Irradiation \n(Including Pump Activation)", color='green', fontsize=20, ha='center')
plt.text(time_irradiation[num_steps // 2], 0.0032, "Irradiation Phase", color='Red', fontsize=20, ha='center')
plt.text(time_after_irradiation[1], 0.0031, "After Irradiation \n(Upper plate openning)", color='Blue', fontsize=20, ha='center')
plt.text(180, 0.0035, "Fh(t)", color='green', fontsize=20, ha='center')
plt.text(250, 0.0035, "Fh(t+dt)", color='green', fontsize=20, ha='center')

# Labels and title
plt.xlabel("Time (days)", fontsize=20)
plt.ylabel("Lateral displacement [m]", fontsize=20)
plt.title("Evolution of the Lateral Displacement of a Fuel Assembly Grid Over Reactor Cycle", fontsize=20)

# Legend
plt.legend(fontsize=20)

# Display the graph
plt.grid(True)
plt.xticks(fontsize=18)  # Adjust the font size for x-axis ticks
plt.yticks(fontsize=18)  # Adjust the font size for y-axis ticks
plt.tight_layout()
plt.savefig("irradiation_example2.png")
plt.show()
plt.close()

