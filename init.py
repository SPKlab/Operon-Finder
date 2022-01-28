from pathlib import Path
import subprocess
from urllib.request import urlretrieve

file_name = "data.7z"

if not Path(file_name).exists:
    subprocess.run(["curl", "https://github.com/tejasvi/operon/releases/download/data/data.7z", "-o", file_name)
    subprocess.run(["atool", "x", file_name])