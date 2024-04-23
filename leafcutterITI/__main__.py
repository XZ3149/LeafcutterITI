from pathlib import Path
from runpy import run_path

pkg_dir = Path(__file__).resolve().parent

def leafcutterITI_clustering():
    script_pth = pkg_dir / "clustering" / "leafcutterITI_clustering.py"
    run_path(str(script_pth), run_name="__main__")
    
def leafcutterITI_map_gen():
    script_pth = pkg_dir / "map_gen" / "leafcutterITI_map_gen.py"
    run_path(str(script_pth), run_name="__main__")
    
    
def leafcutterITI_scITI():
    script_pth = pkg_dir / "scITI" / "leafcutterITI_scITI.py"
    run_path(str(script_pth), run_name="__main__")