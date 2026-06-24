# NACA Optimization Tool & Aerodynamic Suite

## 📖 The Project

This open-source repository was born from the merger of two parallel projects created for the analysis and study of airfoils. Originally developed across different environments (a MATLAB Panel Method solver and a Python XFOIL optimization tool), we have now successfully merged and consolidated the entire workflow into a **single, unified Python suite**.

The result is a powerful and versatile aerodynamic tool that seamlessly bridges algorithmic inverse design (viscous flow via XFOIL) with rapid, detailed inviscid flow field analysis (vortex/source panel method) in one continuous execution.

## ✨ Key Features

* **NACA Airfoil Generation:** Parametric creation of 4-digit NACA airfoils using a cosine point distribution, essential for ensuring stability in boundary layer analysis.
* **Inverse Design (SLSQP Optimization):** Automatic search for optimal geometric parameters (maximum camber *m*, position *p*, thickness *t*) to reach a target $C_L$, respecting a maximum $C_D$ limit set by the user.
* **Smart XFOIL Wrapper:** Complete automation of viscous analysis with autonomous handling of crashes, alpha ramp-up, and data extraction. Automatically detects your OS environment to run on native Windows (`xfoil.exe`) or Linux/WSL (`xfoil`).
* **Integrated Panel Method (Sources & Vortices):** Resolves the potential flow directly in Python. It calculates tangential velocities, pressure coefficient ($C_p$) distribution, and utilizes heavily optimized vectorization (`matplotlib.path`) for instant high-resolution streamline visualization.
* **Export and Logging:** Automatic saving of the final geometry in `.dat` format (ready for CAD or CFD meshing), high-resolution plots, and optimization history in `.csv`.

## 🛠️ Prerequisites & Installation

Ensure you have **Python 3.8+** installed on your system.

Install the required dependencies via pip:

```bash
pip install numpy scipy matplotlib

```

* **Crucial note on XFOIL:** The project requires the XFOIL executable to work.
* **Windows users:** Ensure that the `xfoil.exe` file is present in the same folder as the main script.
* **Linux/WSL users:** Ensure `xfoil` is installed and accessible via your system's PATH, or placed in the working directory.



## 🚀 Usage

Run the main script from your terminal:

```bash
python naca_aero_suite.py

```

The program operates in two phases, guided by a command-line interface:

**Phase 1: Optimization**

1. **Fluid selection:** Air (Standard SL) or Water.
2. **Design parameters:** Speed and chord (for automatic calculation of Mach and Reynolds numbers).
3. **Targets:** Angle of attack, transition parameter $N_{crit}$, desired $C_L$, and maximum tolerated $C_D$.

**Phase 2: Inviscid Analysis**
Once XFOIL finds the optimal geometry, the suite will prompt you to run the advanced Panel Method simulation on the newly found airfoil to visualize the $C_p$ distribution and streamlines.

Upon completion, the tool creates a dedicated results folder (e.g., `Results_Re1500000_Alpha2.0_Cl0.6`) containing all output data, plots, and `.dat` coordinates.

## 🤝 Contributing

We are aerospace engineering students and enthusiasts. Now that the core physics and optimization engines are unified, any pull requests are absolutely welcome! Future implementation ideas include:

* Adding support and parametrization for 5-digit NACA airfoils.
* Implementing a Graphical User Interface (GUI) (e.g., PyQt or Tkinter).
* Extending the aerodynamic calculations to 3D finite wings (Lifting Line Theory).
* Enhancing code robustness and error handling.