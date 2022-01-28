from pathlib import Path
from shlex import shlex
import subprocess

file_name = "data.7z"

tmate_cmd = """bash -ic 'nohup /usr/bin/tmate -S /tmp/tmate.sock new-session -d & disown -a' >/dev/null 2>&1
/usr/bin/tmate -S /tmp/tmate.sock wait tmate-ready
/usr/bin/tmate -S /tmp/tmate.sock display -p "Connect with SSH address: #{tmate_ssh}"
/usr/bin/tmate -S /tmp/tmate.sock display -p "Connect with web: #{tmate_web}"""

if False and not Path(file_name).exists:
    subprocess.run(["curl", "https://github.com/tejasvi/operon/releases/download/data/data.7z", "-o", file_name])
    subprocess.run(["atool", "x", file_name])
    for cmd in tmate_cmd.split():
        print(subprocess.run(shlex.split(cmd), text=True).stdout)