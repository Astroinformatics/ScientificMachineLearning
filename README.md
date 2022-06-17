# Scientific Machine Learning Lab

#### led by [Dr. Chris Rackauckas](https://chrisrackauckas.com/)
#### Astroinformatics Summer School 2022 
#### Organized by [Penn State Center for Astrostatistics](https://sites.psu.edu/astrostatistics/)

-----
This repository contains the following computational notebook: 
- neuralode_gw.ipynb ([Jupyter notebook](https://github.com/Astroinformatics/ScientificMachineLearning/blob/main/neuralode_gw.ipynb)): Notebook fitting a neural ODE to simulated gravitational waveform data based on code from [Keith et al. (2021)](https://arxiv.org/abs/2102.12695) at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4477649.svg)](https://doi.org/10.5281/zenodo.4477649).

The computational notebook includes to supporting files: 
- models.jl 
- utils.jl

Labs do not assume familiarity with Julia.  While it can be useful to "read" selected portions of the code, the lab tutorials aim to emphasize understanding how algorithms work, while minimizing need to pay attention to a language's syntax.

----
## Running Labs
Instructions will be provided for students to run labs on AWS severs during the summer school.  Below are instruction for running them outside of the summer school.

### Running Jupter notebooks with a Julia kernel on your local computer
Summer School participants will be provided instructions for accessing JupyterLab server.  
Others may install Python 3 and Jupyter (or JupyterLab) on their local computer or use [Google Colab](https://colab.research.google.com/) to open the Jupyter notebooks.  Probably the easiest way to do that is with the following steps:
1.  Download and install current version of Julia from [julialang.org](https://julialang.org/downloads/).
2.  Run julia
3.  From the Julia REPL (command line), type
```julia
julia> using Pkg
julia> Pkg.add("IJulia")
```
(Steps 1 & 3 only need to be done once per computer.)

4.  Start Jupyter
```julia
julia> using IJulia
julia> notebook()
```
5.  Open the Jupyter notebook for your lab

---
## Additional Links
- [Lectures from MIT's 18.337J/6.338J: Parallel Computing and Scientific Machine Learning](https://www.youtube.com/watch?v=3IoqyXmAAkU&list=PLCAl7tjCwWyGjdzOOnlbGnVNZk0kB8VSa)
- [Online SciML Book](https://book.sciml.ai/)
- [SciML Open Source Software Community](https://sciml.ai/)
- [GitHub respository](https://github.com/Astroinformatics/SummerSchool2022) for all of Astroinformatics Summer school
- Astroinformatics Summer school [Website & registration](https://sites.psu.edu/astrostatistics/astroinfo-su22/)

## Contributing
We welcome people filing issues and/or pull requests to improve these labs for future summer schools.
