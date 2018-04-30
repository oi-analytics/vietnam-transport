# Vietnam Transport World Bank project

Data preprocessing and model assembly scripts for Vietnam Transport Risk Analysis

## Notebooks in git

Make sure not to commit data inadvertently if working with jupyter notebooks.

Get nbstripout python package:

    conda install -c conda-forge nbstripout

Install git hooks to filter notebooks when committing to git:

    cd /path/to/vietnam-transport
    nbstripout --install
