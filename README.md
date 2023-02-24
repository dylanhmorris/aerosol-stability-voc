# [Comparative aerosol and surface stability of SARS-CoV-2 Variants of Concern](https://doi.org/10.1101/2022.11.21.517352)
Trenton Bushmaker, Kwe Claude Yinda, [Dylan H. Morris](https://dylanhmorris.com), Myndi G. Holbrook, Amandine Gamble, Danielle Adney, Cara Bushmaker, Neeltje van Doremalen, Robert J. Fischer, Raina K. Plowright, James O. Lloyd-Smith, Vincent J. Munster


## Repository information
This repository accompanies the manuscript "Comparative aerosol and surface stability of SARS-CoV-2 Variants of Concern" (T Bushmaker et al, 2023). It provides code for reproducing data analysis and recreating display figures from the paper

## License and citation information
If you use the code or data provided here, please make sure to do so in light of the project [license](LICENSE.txt) and please cite our work:

- Bushmaker, Trenton, et al. "Comparative aerosol and surface stability of SARS-CoV-2 Variants of Concern." bioRxiv (2022): 2022-11. https://doi.org/10.1101/2022.11.21.517352

Bibtex record:
```
@article{bushmaker2022comparative,
  title={Comparative aerosol and surface stability of SARS-CoV-2 Variants of Concern},
  author={Bushmaker, Trenton and Yinda, Kwe Claude and Morris, Dylan and Holbrook, Myndi and Gamble, Amandine and Adney, Danielle and Bushmaker, Cara and van Doremalen, Neeltje and Fischer, Robert and Plowright, Raina and others},
  journal={bioRxiv},
  pages={2022--11},
  year={2022},
  publisher={Cold Spring Harbor Laboratory},
  doi={10.1101/2022.11.21.517352}
}
```

## Article abstract 
SARS-CoV-2 is transmitted principally via air; contact and fomite transmission may also occur. Variants-of-concern (VOCs) are more transmissible than ancestral SARS-CoV-2. We find that early VOCs show greater aerosol and surface stability than the early WA1 strain, but Delta and Omicron do not. Stability changes do not explain increased transmissibility.


## Directories
- [``varstab``](varstab): R package containing useful functions and styling called in the analysis scripts.
- [``src``](tree/main/src): Analysis code including data preprocessing, Bayesian model definition and fitting, and results post-processing and figure generation. 
  - [``src/stan``](/src/stan) contains Bayesian model specification in Stan code.
  - [``src/parameters``](src/parameters) specifies parameters for the analysis.



## Reproducing analysis

A guide to reproducing the analysis from the paper follows. If you encounter issues, see the **Troubleshooting** section at the end of this README.

### Getting the code
First download this repository. The recommended way is to ``git clone`` it from the command line:

    git clone https://github.com/dylanhmorris/sars-cov-2-temp-humidity.git

Downloading it manually via Github's download button should also work.

### Dependency installation
The analysis can be auto-run via the project ``Makefile``, but you may need to install some external dependencies first. See the **Dependency installation guide** below for a complete walkthrough. In the first instance, you'll need a working installation of the statistical programming language R, a working C++ compiler, and a working installation of Gnu Make or similar. A few external R packages can then be installed from the command line by typing.

    make install

from within the project directory. This will install any missing external R packages as well as the project package ``varstab``.

### Running the analysis

The simplest approach is simply to type ``make`` at the command line, which should produce a full set of figures, tables, MCMC output (saved as R Dataset ``.Rds`` files in the ``out/chains/`` directory as ``<model_name>-chains.Rds``), etc. The MCMC output can be loaded in any working R installation, as long as the package ``rstan`` is also installed.

If you want to do things piecewise, typing ``make <filename>`` for any of the files listed in the ``dat/cleaned`` or ``out`` directories below should run the steps needed to produce that file.

Some shortcuts are available:

- ``make install`` installs dependencies and ``varstab``
- ``make uninstall`` removes ``varstab``
- ``make data`` produces cleaned data files.
- ``make chains`` produces all MCMC output, including prior predictive checks
- ``make figures`` produces all figures
- ``make tables`` produces all tables
- ``make clean`` removes all generated files, leaving only source code (though it does not uninstall packages)
- ``make rebuild`` runs ``make clean``, ``make uninstall`` and ``make install``, and then ``make``, to regenerate the whole analysis from scratch.

### Examining code

Examining the raw Stan code is the place to start to understand how Bayesian inference models have been specified. But note that parameters for the prior distributions are set at runtime rather than hard-coded into the ``.stan`` files, so that recompilation is not required when parameter choices are changed (this makes it easier to try the models using different priors, e.g. for sensitivity analysis).

Prior parameter choices are specified in the [``src/parameters``](src/parameters) directory.

Stan model code is modularized. Top level model files (found in [``src/stan``](src/stan) are built in part by including functions found in [``src/stan/functions``](src/stan/functions).

## Project structure when complete

Once the full analysis has been run, you should be able to find a full set of mcmc output in [``out/chains``](out/chains), figures in [``out/figures``](out/figures), and tables in [``out/tables``](out/tables).

## Dependency installation guide
You will need a working R installation with the command line interpreter ``Rscript`` (macOS and Linux) or ``Rscript.exe`` (Windows). On mac and Linux, you can check that you have an accessible ``Rscript`` by typing ``which Rscript``at the command line and seeing if one is found.

If you do not have an R installation, you can install it from [the R project website](https://www.r-project.org/) or from the command line using a package manager such as [Homebrew](https://brew.sh/) on macOS or ``apt-get`` on Linux. macOS users may also need to install the macOS "command line tools" by typing ``xcode-select --install`` at a command prompt.

Once R and ``make`` are installed, you can automatically install all other dependencies (including the Hamiltonian Monte Carlo software Stan and its R interface rstan) on most systems using ``make``. In the top level project directory, type the following at the command line:

    make install

Alternatively, you can run the script ``src/install.R`` manually. 

Note that installing Stan and RStan can be time-consuming. Some of the packages in the very valuable [tidyverse](https://www.tidyverse.org/) may also take some time to install.
