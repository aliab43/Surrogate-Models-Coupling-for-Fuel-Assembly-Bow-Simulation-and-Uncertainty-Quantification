from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

BASE_DIR = Path(__file__).resolve().parent
DEFAULT_OUTPUT = BASE_DIR / "irradiation_example2.png"


def build_cycle_profile():
    """Build a simple illustrative displacement profile over the reactor cycle."""
    time_before_irradiation = [0, 7]
    displacement_before_irradiation = [0.002, 0.003]

    num_steps = 18
    step_duration = 30
    time_irradiation = [7 + step_duration * index for index in range(num_steps)]

    base_displacement = 0.003
    linear_increment_per_step = 0.00005
    amplitude = 0.0002
    frequency = 2 * np.pi / len(time_irradiation)

    displacement_irradiation = [
        base_displacement + linear_increment_per_step * index + amplitude * np.sin(frequency * index)
        for index in range(len(time_irradiation))
    ]

    time_after_irradiation = [time_irradiation[-1], time_irradiation[-1] + 30]
    displacement_after_irradiation = [displacement_irradiation[-1], 0.0028]

    time_combined = time_before_irradiation + time_irradiation + time_after_irradiation
    displacement_combined = (
        displacement_before_irradiation + displacement_irradiation + displacement_after_irradiation
    )
    return time_combined, displacement_combined, time_before_irradiation, time_irradiation, time_after_irradiation


def plot_cycle_profile(output_path=DEFAULT_OUTPUT, show=True):
    """Plot the illustrative cycle profile and save it next to this script by default."""
    (
        time_combined,
        displacement_combined,
        time_before_irradiation,
        time_irradiation,
        time_after_irradiation,
    ) = build_cycle_profile()

    fig, ax = plt.subplots(figsize=(16, 10))
    ax.plot(time_combined, displacement_combined, color="b", marker="o", linestyle="-", linewidth=2)

    ax.axvline(
        x=time_before_irradiation[1],
        color="crimson",
        linestyle="--",
        linewidth=3,
        label="Power on",
    )
    ax.axvline(
        x=time_irradiation[-1],
        color="black",
        linestyle="--",
        linewidth=3,
        label="Power off",
    )
    ax.text(
        35,
        0.0026,
        "Before Irradiation\n(Including Pump Activation)",
        color="green",
        fontsize=20,
        ha="center",
    )
    ax.text(
        time_irradiation[len(time_irradiation) // 2],
        0.0032,
        "Irradiation Phase",
        color="red",
        fontsize=20,
        ha="center",
    )
    ax.text(
        time_after_irradiation[1],
        0.0031,
        "After Irradiation\n(Upper plate opening)",
        color="blue",
        fontsize=20,
        ha="center",
    )
    ax.text(180, 0.0035, "Fh(t)", color="green", fontsize=20, ha="center")
    ax.text(250, 0.0035, "Fh(t+dt)", color="green", fontsize=20, ha="center")

    ax.set_xlabel("Time (days)", fontsize=20)
    ax.set_ylabel("Lateral displacement [m]", fontsize=20)
    ax.set_title(
        "Evolution of the Lateral Displacement of a Fuel Assembly Grid Over a Reactor Cycle",
        fontsize=20,
    )
    ax.legend(fontsize=20)
    ax.grid(True)
    ax.tick_params(axis="x", labelsize=18)
    ax.tick_params(axis="y", labelsize=18)

    output_path = Path(output_path)
    fig.tight_layout()
    fig.savefig(output_path)

    if show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    plot_cycle_profile()
