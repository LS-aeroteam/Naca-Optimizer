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
    
    missing_libs = []
    for lib, reason in required_libraries.items():
        try:
            __import__(lib)
            print(f"    - {lib}: Found.")
        except ImportError:
            print(f"    - {lib}: NOT FOUND. This library is required {reason}.")
            missing_libs.append(lib)
            
    if missing_libs:
        print("\n[i] Missing libraries detected. Attempting to install them automatically...")
        try:
            subprocess.check_call([sys.executable, "-m", "pip", "install"] + missing_libs)
            print("    [+] Successfully installed missing libraries.\n")
        except subprocess.CalledProcessError as e:
            print(f"\n[!] Failed to install libraries automatically. Error: {e}")
            print(f"    Please install manually: pip install {' '.join(missing_libs)}")
            sys.exit(1)
    else:
        print("    All Python libraries are installed.\n")

def perform_all_checks():
    """Runs all pre-run checks for dependencies and environment."""
    print("======================================================================")
    print("               RUNNING PRE-FLIGHT ENVIRONMENT CHECKS")
    print("======================================================================")
    check_python_libraries()
    print("======================================================================")
    print("                  ALL CHECKS PASSED. SYSTEM IS READY.")
    print("======================================================================")

