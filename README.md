# mpoxclade1_reffprojection
Source code accompanying Endo et al. "Time to discuss reinstating an orthopoxvirus immunisation programme". <br>
## Instructions
* Codes are reproducible via [Codespaces](https://github.com/features/codespaces) or [Docker Desktop](https://www.docker.com/products/docker-desktop/) through Compose plugin. Docker build may take up to about 30 mins.
* Use [Jupytext](https://jupytext.readthedocs.io/en/latest/) commands to generate Jupyter Notebooks from the source files in /notebooks. Main analysis may take up to about 30 mins.
* If using VSCode-based Codespaces, install [Julia](https://marketplace.visualstudio.com/items?itemName=julialang.language-julia) and [Jupyter](https://marketplace.visualstudio.com/items?itemName=ms-toolsai.jupyter) extensions before opening the notebook.
* If Jupyter produces a kernel error, launch Julia from the terminal and try `]activate .` and `]instantiate` before reopening the notebook.
* See ./devcontainer/Dockerfile and Manifest.toml for environments and software dependencies
## Author
Akira Endo, Hiroaki Murayama, Shihui Jin, Borame L Dickens
