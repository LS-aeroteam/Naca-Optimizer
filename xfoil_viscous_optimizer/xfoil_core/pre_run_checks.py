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

# --- XFOIL Executable Check ---

def find_xfoil_executable():
    """
    Searches for the XFOIL executable in the project root and system PATH.

    Returns:
        str: The full path to the XFOIL executable if found, otherwise None.
    """
    xfoil_name = "xfoil.exe" if os.name == 'nt' else "xfoil"

    # 1. Check in the project's root directory (one level up from this file's package)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    local_path = os.path.join(project_root, xfoil_name)
    if os.path.isfile(local_path):
        print(f"    - Found XFOIL executable in project root: {local_path}")
        return local_path
    
    # 2. Check in system's PATH
    system_path = shutil.which(xfoil_name)
    if system_path:
        print(f"    - Found XFOIL executable in system PATH: {system_path}")
        return system_path
        
    return None

def check_xfoil():
    """
    Checks for the XFOIL executable and provides instructions if not found.
    If on Windows and missing, attempts to download it automatically.
    Exits the program if XFOIL cannot be found or downloaded.
    """
    print("[+] Checking for XFOIL executable...")
    xfoil_path = find_xfoil_executable()
    
    if xfoil_path is None:
        print("\n[!] XFOIL executable not found.")
        
        if os.name == 'nt':
            print("    [i] Attempting to download XFOIL automatically for Windows...")
            try:
                url = "https://web.mit.edu/drela/Public/web/xfoil/XFOIL6.99.zip"
                script_dir = os.path.dirname(os.path.abspath(__file__))
                project_root = os.path.dirname(script_dir)
                
                with tempfile.TemporaryDirectory() as temp_dir:
                    zip_path = os.path.join(temp_dir, "xfoil.zip")
                    urllib.request.urlretrieve(url, zip_path)
                    
                    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                        zip_ref.extractall(temp_dir)
                    
                    # Find the exe in the extracted contents and move it to project root
                    exe_found = False
                    for root, dirs, files in os.walk(temp_dir):
                        for file in files:
                            if file.lower() == 'xfoil.exe':
                                src_exe = os.path.join(root, file)
                                dest_exe = os.path.join(project_root, 'xfoil.exe')
                                shutil.copy(src_exe, dest_exe)
                                exe_found = True
                                break
                        if exe_found:
                            break
                            
                if exe_found:
                    print("    [+] Successfully downloaded and installed XFOIL!")
                    xfoil_path = find_xfoil_executable()
                else:
                    raise Exception("Downloaded zip did not contain xfoil.exe.")
            except Exception as e:
                print(f"    [!] Failed to download XFOIL automatically: {e}")
                sys.exit(1)
        else:
            print("    This program requires the XFOIL command-line tool for aerodynamic calculations.")
            print("\n    To fix this, please:")
            print("    1. Install XFOIL for your operating system (e.g., sudo apt-get install xfoil).")
            print("    2. Alternatively, place the executable in the main project directory.")
            sys.exit(1) # Exit the program
    
    # Optional: Test if the executable can run
    try:
        proc = subprocess.run([xfoil_path], input='\n\nQUIT\n', capture_output=True, text=True, timeout=10)
        if "XFOIL" not in proc.stdout and "XFOIL" not in proc.stderr:
             raise IOError("XFOIL did not respond as expected.")
        print("    - XFOIL executable is valid and responsive.")
    except (subprocess.TimeoutExpired, IOError, OSError) as e:
        print(f"\n[!] Found XFOIL at '{xfoil_path}', but it failed to run correctly.")
        print(f"    Error: {e}")
        print("    Please ensure it is a valid, executable file and that you have permission to run it.")
        sys.exit(1)

    print("    XFOIL dependency check passed.\n")


def perform_all_checks():
    """Runs all pre-run checks for dependencies and environment."""
    print("======================================================================")
    print("               RUNNING PRE-FLIGHT ENVIRONMENT CHECKS")
    print("======================================================================")
    check_python_libraries()
    check_xfoil()
    print("======================================================================")
    print("                  ALL CHECKS PASSED. SYSTEM IS READY.")
    print("======================================================================")

