from pathlib import Path
from runpy import run_path

pkg_dir = Path(__file__).resolve().parent

def tealeaf_clustering():
    script_pth = pkg_dir / "clustering" / "tealeaf_clustering.py"
    run_path(str(script_pth), run_name="__main__")
    
def tealeaf_map_gen():
    script_pth = pkg_dir / "map_gen" / "tealeaf_map_gen.py"
    run_path(str(script_pth), run_name="__main__")
    
    
def tealeaf_sc():
    script_pth = pkg_dir / "sc" / "tealeaf_sc.py"
    run_path(str(script_pth), run_name="__main__")
    

    
def tealeaf_ggsashimi():
    script_pth = pkg_dir / "ggsashimi" / "tealeaf_ggsashimi.py"
    run_path(str(script_pth), run_name="__main__")