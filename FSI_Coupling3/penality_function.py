import numpy as np
import matplotlib.pyplot as plt

# Define the constant lambda
lambda_value = 2 * 10**4

# Define the function f(x)
def f(x):
    return 500 * (1 - np.exp(-lambda_value * x))

# Generate the x-values
x_values = np.linspace(0, 1e-3, 100)  # from 0 to 10^-3
y_values = f(x_values)

# Create the plot
plt.figure(figsize=(12, 6))
plt.plot(x_values, y_values, label=r"$f(x) = 500(1 - e^{-\lambda x})$", color='b')

# Add a horizontal line for the asymptotic value
asymptotic_value = 500
plt.axhline(y=asymptotic_value, color='red', linestyle='--', label='Asymptotic Value = 500')

# Annotate the plot with larger text
plt.text(0.00018, 300, 
         "Threshold of 500 [N] to avoid computational instability\n\n"
         r"The slope at 0 should be equal to $10^7 N/m$"
         "\n\nfor contact precision",
         fontsize=19, color='black', ha='left')

# Set labels and title with larger font sizes
plt.xlabel("Penetration Depth [m]", fontsize=22)
plt.ylabel("Penalty Force [N]", fontsize=22)
plt.title("Non-linear Penalty Function", fontsize=26)

# Increase tick label font size
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

# Increase legend font size
plt.legend(fontsize=18)

# Grid and show
plt.grid(True)
plt.savefig("Non_Linear_Penality.png")
plt.show()
