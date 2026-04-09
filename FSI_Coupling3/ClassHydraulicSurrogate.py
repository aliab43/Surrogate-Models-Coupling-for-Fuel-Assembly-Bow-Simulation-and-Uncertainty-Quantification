import sys
import importlib
import numpy as np

class HydraulicSurrogateModel:
    def __init__(self):
        """
        Initialize the HydraulicSurrogateModel class.
        """
        # Path to the directory containing the surrogate model modules
        self.model_path = "TGaussianProcess_Hydraulic/"
        sys.path.append(self.model_path)

        # Grid indices (only 2 to 9 are considered for force calculation)
        self.grid_indices = range(2, 10)

        # Assembly indices (1 to 15)
        self.assembly_indices = range(1, 16)

    def predict_force_field(self, x):
        """
        Predict hydraulic forces for all grid and assembly indices and organize them into a matrix.

        Args:
            x (list or np.ndarray): Input to the surrogate model.

        Returns:
            np.ndarray: A matrix where rows represent grid indices and 
                        columns represent assembly indices.
        """
        # Initialize the result matrix (10 rows for grids, 15 columns for assemblies)
        results = np.zeros((10, 15))

        # Fill the matrix for grids 2 to 9
        for i in self.grid_indices:  # Loop through grid indices 2 to 9
            for j in self.assembly_indices:  # Loop through assembly indices 1 to 15
                # Dynamically construct the class and function names
                class_name = f"Krig_F_h_G_{i}_A_{j}"
                function_name = f"kriging_function_F_h_G_{i}_A_{j}"

                try:
                    # Dynamically import the module
                    module = importlib.import_module(class_name)

                    # Retrieve the class from the module
                    KrigingClass = getattr(module, class_name)

                    # Instantiate the class
                    kriging_model = KrigingClass()

                    # Call the prediction function
                    if hasattr(kriging_model, function_name):
                        prediction = getattr(kriging_model, function_name)(x)
                        # Assign the prediction to the correct position in the matrix
                        results[i - 1, j - 1] = prediction  # Adjust indices for Python (0-based)
                    else:
                        print(f"Function {function_name} not found in {class_name}")
                except ModuleNotFoundError:
                    print(f"Module {class_name} not found.")
                except AttributeError as e:
                    print(f"Error accessing {class_name} or {function_name}: {e}")
                except Exception as e:
                    print(f"An error occurred while processing {class_name}: {e}")

        # Return the transpose of the result matrix for proper alignment
        return results
