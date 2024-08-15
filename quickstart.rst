===============
Getting Started
===============

Here, we have provided a quick start guide that will allow you to get up and running quickly and able to project PTMs on any splicing dataset you might have

Installation
------------

KSTAR can be installed via `pip`, tarball, and directly from the Git repository. We recommend using pip to install the most well-tested version of the package, but check our development branch on the GitHub repository to see some of the additional analyses/data we are adding to the package!

==================================== ================================================================================
Install Method                       Code
==================================== ================================================================================
pip                                  .. code-block:: bash

                                        pip install ptm-pose

Github Release                       .. code-block:: bash

                                        wget https://github.com/NaegleLab/PTM-POSE/archive/refs/tags/<version>.tar.gz
                                        tar -xvf PTM-POSE-<version>.tar.gz
                                        cd PTM-POSE-<version>
                                        python setup.py install

Git Clone (for development version)  .. code-block:: bash

                                        git clone https://github.com/NaegleLab/PTM-POSE
                                        cd PTM-POSE
                                        git checkout dev
                                        python setup.py install

==================================== ================================================================================



Configuring your PTM-POSE environment
-------------------------------------

After installing PTM-POSE, you will need to download:

ptm_coordinates
    contains the genomic coordinates of all post-translational modifications in the proteome based on data from ProteomeScout and PhosphoSitePlus. This is available through the GitHub repository large file storage.
translator
    allows for quick conversion between UniProt IDs and other database IDs (namely Ensembl Gene ID and gene names). This is important for adding annotations to your data. This file is created using UniProt's Rest API.


In both cases, we have provided quick functions for downloading these files. While they do not automatically get saved to your machine, as we know that for some this may not be the desired behavior, we recommended running this functions with `save=True`, which will save each file in the directory where PTM-POSE is stored and keep you from needed to run these functions in the future.

.. code-block:: python

    from ptm_pose import pose_config

    ptm_coordinates = pose_config.download_ptm_coordinates()
    translator = pose_config.download_translator()


Once these have been downloaded, you are ready to start projecting PTMs onto your data!

