# generate-env.py

import os
from pathlib import Path

env_vars = {}
# raccogli tutte le env LRGSG_*
with open('.env', "r") as f:
    for envar in f:
        list_split = envar.split("=")
        name = list_split[0]
        value = list_split[1].strip()
        env_vars[name] = value

# scrivi il modulo
with open(Path('src', 'lrgsglib', 'config') / "lrgsg_env.py", "w") as f:
    f.write("from pathlib import Path\n\n")
    for name, path in env_vars.items():
        # definisce: LRGSG_DATA: Path = Path("/.../data")
        f.write(f"{name}: Path = Path({str(path)!r})\n")