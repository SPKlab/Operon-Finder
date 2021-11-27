"""Delete non-allowed organism files from current directory"""
from pathlib import Path
from shutil import rmtree
allowed = {'511145', '196627', '298386', '224308', '297246', '169963', '85962', '272634', '83332', '262316'}
if allowed.intersection(f.name.rsplit('.')[0] for f in Path('.').iterdir()):
    for f in Path('.').iterdir():
        if f.name.rsplit('.')[0] not in allowed:
            if f.is_dir():
                rmtree(f)
            else:
                f.unlink()
else:
    print("no intersect")