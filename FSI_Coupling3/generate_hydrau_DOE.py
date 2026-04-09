from pathlib import Path

import numpy as np

BASE_DIR = Path(__file__).resolve().parent


class HydraulicDoeGenerator:
    """Generate hydraulic DOE samples in both text and NumPy formats."""

    def __init__(self, n_params=None, Nb_params=None):
        self.n_params = n_params if n_params is not None else Nb_params
        if self.n_params is None:
            raise ValueError("Provide `n_params` or the legacy `Nb_params` value.")

    def generate_one_sample(self):
        temperature = np.random.normal(300.0, 3.33)
        mean_velocity = np.random.normal(5.0, 0.05)
        max_deviation_in = np.random.normal(0.05, 0.005)
        offset_in = np.random.normal(0.0, 0.02)
        max_deviation_out = np.random.normal(0.04, 0.004)
        offset_out = np.random.normal(0.0, 0.01)
        cg_sensitivity = np.random.uniform(1.0, 1.232)
        hl_sensitivity = np.random.uniform(0.01, 0.03)
        return [
            temperature,
            mean_velocity,
            max_deviation_in,
            offset_in,
            max_deviation_out,
            offset_out,
            cg_sensitivity,
            hl_sensitivity,
        ]

    def generate_doe(
        self,
        n_samples,
        n_assemblies,
        dat_filename="DOE_for_surrogate_hydraulic.dat",
        npy_filename="DOE_for_mechanics/Doe_for_hydraulic.npy",
    ):
        dat_path = BASE_DIR / dat_filename
        npy_path = BASE_DIR / npy_filename
        npy_path.parent.mkdir(parents=True, exist_ok=True)

        results_matrix = []

        with dat_path.open("w") as file:
            text_format_c = "C_{price:d} | "
            text_format_s = "S_{price:d} | "
            text_format_w = "W_{price:d} | "
            text_format_d = " D |"

            file.write("#NAME: tds_Phorcys \n")
            file.write("#COLUMN_NAMES: ")
            for assembly_index in range(n_assemblies):
                file.write(text_format_c.format(price=assembly_index + 1))
            for assembly_index in range(n_assemblies):
                file.write(text_format_s.format(price=assembly_index + 1))
            for assembly_index in range(n_assemblies):
                file.write(text_format_w.format(price=assembly_index + 1))

            file.write(
                "T | v_m | max_deviation_in | l_offset_in | max_deviation_out | "
                "l_offset_out | Cg_sensi | hl_sensi | tdsEstim_hydrau__n__iter__ \n"
            )
            file.write("#COLUMN_TYPES: ")
            for _ in range(self.n_params):
                file.write(text_format_d)
            file.write(" \n\n")

            for sample_index in range(n_samples):
                sample = self.generate_one_sample()
                for value in sample:
                    file.write(f"{value} ")
                file.write(f"{float(sample_index)} \n")
                results_matrix.append(sample)

        np.save(npy_path, results_matrix)
        print(results_matrix)
        print(f"Matrix saved to {npy_path}")

    # Backward-compatible wrapper for the original research scripts.
    def generateDoeHydrau(self, Ns, nb_assembly):
        self.generate_doe(Ns, nb_assembly)


GenerateDoeHydrau = HydraulicDoeGenerator


if __name__ == "__main__":
    HydraulicDoeGenerator(Nb_params=57).generate_doe(n_samples=1, n_assemblies=16)
