# Getting Started

Here, we have provided a quick start guide that will allow you to get up and running quickly and able to project PTMs on any splicing dataset you might have

## Installation

KSTAR can be installed via `pip`, tarball, and directly from the Git repository. We recommend using pip to install the most well-tested version of the package, but check our development branch on the GitHub repository to see some of the additional analyses/data we are adding to the package!

#### Pip

To install via pip, execute `pip install ptm-pose`.


#### Tarball

To install via a tarball, head over to the [Releases page](https://github.com/NaegleLab/PTM-POSE/releases/) and download the latest stable tar release.

Afterwards, navigate to your downloads directory and execute the following commands, substituting <version> for the release's version number:
```
tar -xvf PTM-POSE-<version>.tar.gz
cd PTM-POSE-<version>
python setup.py install
```

#### Git

If you want to try out the latest commit, you can install directly from the Git repository by executing the following commands:
```
git clone https://github.com/NaegleLab/PTM-POSE
cd PTM-POSE
python setup.py install
```

## Configuring your PTM-POSE environment

After installing PTM-POSE, you will need to download:
1. ptm_coordinates dataframe, which contains the genomic coordinates of all post-translational modifications in the proteome based on data from ProteomeScout and PhosphoSitePlus. This is available through the GitHub repository large file storage.
2. Translator file, which allows for quick conversion between UniProt IDs and other database IDs (namely Ensembl Gene ID and gene names). This is important for adding annotations to your data. This file is created using UniProt's Rest API.


In both cases, we have provided quick functions for downloading these files. While they do not automatically get saved to your machine, as we know that for some this may not be the desired behavior, we recommended running this functions with `save=True`, which will save each file in the directory where PTM-POSE is stored and keep you from needed to run these functions in the future.

```python
from ptm_pose import pose_config

ptm_coordinates = pose_config.download_ptm_coordinates(save = True)
translator = pose_config.download_translator(save = True)
```


## Obtaining annotation information from various database sources

Beyond projecting PTMs onto your data, we have also provided additional functions for appending information on the function, relationships, and interactions of each post-translational modification that have been recorded in various databases. These annotations include information from:

- [PhosphoSitePlus](https://www.phosphosite.org/homeAction.action)
- [DEPOD](https://depod.bioss.uni-freiburg.de/)
- *[RegPhos](http://140.138.144.141/~RegPhos/index.php)
- *[ELM](http://elm.eu.org/)
- *[PTMcode](https://ptmcode.embl.de/)

* indicate database information that can be downloaded directly with the PTM-POSE package, although for fastest processing its recommended that you download the information manually.

We are continuing to work on adding functions to append more contextual information for individual PTMs. If you have suggestions for what information you would like to be added, please let us know!


## How to use PTM-POSE

To get a feel for the structure of PTM-POSE and the data it requires, its recommend you first read through the [general instructions on how to use PTM-POSE](). We have also provided additional examples on how this looks in practice using data obtained from MATS or from SpliceSeq, which will give you an idea of the types of analysis that can be down with the information outputted by PTM-POSE.