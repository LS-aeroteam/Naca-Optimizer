import os
import shutil
import sys
import urllib.request
import zipfile
import tempfile
import subprocess

# --- Python Library Checks ---

def check_python_libraries():
    """
    Checks if essential Python libraries are installed.
    Exits with an informative message if a library is missing.
    """
    print("[+] Checking required Python libraries...")
    required_libraries = {
        "numpy": "for numerical operations",
        "scipy": "for optimization functions",
        "matplotlib": "for plotting results"
    }
    
    missing_count = 0
    for lib, reason in required_libraries.items():
        try:
            __import__(lib)
            print(f"    - {lib}: Found.")
        except ImportError:
            print(f"    - {lib}: NOT FOUND. This library is required {reason}.")
            missing_count += 1
            
    if missing_count > 0:
        print("\n[!] Please install the missing libraries to proceed.")
        print(f"    You can install them by running: pip install numpy scipy matplotlib")
        sys.exit(1) # Exit the program
    print("    All Python libraries are installed.\n")

def perform_all_checks():
    """Runs all pre-run checks for dependencies and environment."""
    print("======================================================================")
    print("               RUNNING PRE-FLIGHT ENVIRONMENT CHECKS")
    print("======================================================================")
    check_python_libraries()
    print("    - XFOIL check skipped (Using In-House Panel Method Solver).")
    print("======================================================================")
    print("                  ALL CHECKS PASSED. SYSTEM IS READY.")
    print("======================================================================")

