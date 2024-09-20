===============
Installation
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

conda                                .. code-block:: bash

                                        conda install -c naeglelab ptm-pose
                                        conda install -c bioconda gseapy   

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




