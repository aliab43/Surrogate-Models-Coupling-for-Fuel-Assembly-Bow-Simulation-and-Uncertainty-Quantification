import importlib
import sys
from pathlib import Path

import numpy as np

BASE_DIR = Path(__file__).resolve().parent


class HydraulicSurrogateModel:
    """Load hydraulic kriging exports and assemble the force matrix."""

    def __init__(self):
        self.model_path = BASE_DIR / "TGaussianProcess_Hydraulic"
        model_path_str = str(self.model_path)
        if model_path_str not in sys.path:
            sys.path.append(model_path_str)

        self.grid_indices = range(2, 10)
        self.assembly_indices = range(1, 16)

    def predict_force_field(self, x):
        """Return the hydraulic force field as a `(10, 15)` array."""
        results = np.zeros((10, 15))

        for grid_index in self.grid_indices:
            for assembly_index in self.assembly_indices:
                class_name = f"Krig_F_h_G_{grid_index}_A_{assembly_index}"
                function_name = f"kriging_function_F_h_G_{grid_index}_A_{assembly_index}"

                try:
                    module = importlib.import_module(class_name)
                    kriging_class = getattr(module, class_name)
                    kriging_model = kriging_class()

                    if hasattr(kriging_model, function_name):
                        prediction = getattr(kriging_model, function_name)(x)
                        results[grid_index - 1, assembly_index - 1] = prediction
                    else:
                        print(f"Function {function_name} not found in {class_name}")
                except ModuleNotFoundError:
                    print(f"Module {class_name} not found.")
                except AttributeError as error:
                    print(f"Error accessing {class_name} or {function_name}: {error}")
                except Exception as error:
                    print(f"An error occurred while processing {class_name}: {error}")

        return results
