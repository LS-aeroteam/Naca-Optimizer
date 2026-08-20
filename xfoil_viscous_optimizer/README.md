# XFOIL Viscous Optimizer

This module hosts the project's legacy code. Unlike the In-House solver (which works in potential flow), this tool leverages a Python wrapper to launch and automate the well-known **XFOIL** executable in the background.

## When to use it?
You should use this module whenever you need high-fidelity viscous estimates. Since the code calls XFOIL, the optimization will take into account not only Lift (Cl) but also Drag (Cd) and boundary layer effects, offering a more complete simulation at the expense of execution speed (opening and closing XFOIL for each airfoil generates significant overhead).

## Usage
Make sure you are in this folder and run the main script:
```bash
python run.py
```

### Pre-Run Automation
If the program does not find the XFOIL executable on your system, it will attempt to download it automatically and prepare it for use (currently supported on Windows). Additionally, it will **automatically install** any missing Python libraries without prompting.

### Workflow
By following the on-screen instructions, you can set your parameters and start the viscous optimization. The input prompts have been unified with the In-House solver for a seamless experience. 
Once launched, the process runs from start to finish without pausing for confirmations. All final aerodynamic plots, XFOIL polars, and best-candidate geometries are neatly exported into a categorized subfolder inside the `Results/` directory.
