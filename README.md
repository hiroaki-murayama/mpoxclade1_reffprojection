# mpoxclade1_eigen
Code accompanying [Murayama & Asakura et al. "Roles of community and sexual contacts as drivers of clade I mpox outbreaks"](https://www.medrxiv.org/content/10.1101/2024.10.15.24315554)

## Instructions
* Codes are reproducible via [Codespaces](https://github.com/features/codespaces) or [Docker Desktop](https://www.docker.com/products/docker-desktop/) through Compose plugin. Docker build may take up to about 30 mins.
* Use [Jupytext](https://jupytext.readthedocs.io/en/latest/) commands to generate Jupyter Notebooks from the source files in /notebooks. Main analysis may take up to about 30 mins.
* If using VSCode-based Codespaces, install [Julia](https://marketplace.visualstudio.com/items?itemName=julialang.language-julia) and [Jupyter](https://marketplace.visualstudio.com/items?itemName=ms-toolsai.jupyter) extensions before opening the notebook.
* If Jupyter produces a kernel error, launch Julia from the terminal and try `]activate .` and `]instantiate` before reopening the notebook.
* See ./devcontainer/Dockerfile and Manifest.toml for environments and software dependencies
