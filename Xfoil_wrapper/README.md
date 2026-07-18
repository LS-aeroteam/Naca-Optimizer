# XFOIL Wrapper

This module hosts the project's legacy code. Unlike the In-House solver (which works in potential flow), this tool leverages a Python wrapper to launch and automate the well-known **XFOIL** executable in the background.

## When to use it?
You should use this module whenever you need high-fidelity viscous estimates. Since the code calls XFOIL, the optimization will take into account not only Lift (Cl) but also Drag (Cd) and boundary layer effects, offering a more complete simulation at the expense of execution speed (opening and closing XFOIL for each airfoil generates significant overhead).

## Usage
Make sure you are in this folder and run the main script:
```bash
python run_suite.py
```
If the program does not find the XFOIL executable on your system, it will attempt to download it automatically and prepare it for use (currently supported on Windows). By following the on-screen instructions, you can set your parameters and start the viscous optimization.
