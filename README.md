# LRGSGlib

LRGSGlib is a set of Python modules and C/C++ extensions implementing the theoretical tools of the **Laplacian Renormalisation Group for Signed Graphs**.  It provides utilities for building signed networks, running renormalisation flows and simulating statistical physics models such as the Ising, contact process and voter dynamics.  Additional helpers for plotting, logging and networking patches are also included.

The C sources for the performance critical parts live in `src/lrgsglib/Ccore` and can be built along with the Python modules.

## Installation

1. **Clone the repository and initialise submodules**
   ```bash
   git clone https://github.com/giulioiannelli/lrgsglib
   cd lrgsglib
   git submodule init
   git submodule update
   ```
2. **Create the conda environment** (default name: `lrgsgenv`)
   ```bash
   conda env create -f lrgsgenv.yml
   conda activate lrgsgenv
   ```
3. **Build the project**
   ```bash
   make all
   ```
4. **Install in editable mode**
   ```bash
   pip install -e .
   ```

## Git helper

To keep the `main` branch in sync with `with-ipynb` but exclude the `ipynb/` directory, add the following alias to your `.gitconfig`:

```ini
[alias]
    sync-with-ipynb-to-main = "!sh -c '\
        trap \"git checkout with-ipynb\" EXIT; \
        git checkout main && \
        git pull origin main && \
        git checkout with-ipynb -- . :\(exclude\)ipynb/** && \
        git add -A && \
        if ! git diff --cached --quiet; then \
            git commit -m \"Sync from with-ipynb (exclude ipynb/)\" && \
            git push origin main; \
        else \
            echo \"No edits\"; \
        fi\
    '"
```

## More information

See the READMEs under `src/lrgsglib/Ccore` for details about the available C programs and statistical physics models.

