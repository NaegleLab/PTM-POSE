from setuptools import setup, find_packages

# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

# Run setup
setup(
    name="ptm_pose",
    version="0.3.1",
    author="Naegle Lab",
    author_email="kmn4mj@virginia.edu",
    url="https://github.com/NaegleLab/PTM-POSE/tree/main",
    install_requires=['pandas>=2.2', 'numpy>=1.26', 'biopython>=1.78', 'tqdm>=4.66', 'scipy>=1.13', 'xlrd', 'matplotlib','seaborn', 'requests'],
    extras_require={
        'full_analysis': ['omnipath>=1.0.0', 'pyliftover>=0.4', 'nease>=1.3', 'gseapy>=1.1','kinase-library>=1.3', 'networkx>=3.3']
    },
    keywords=['PTM', 'splicing', 'bioinformatics', 'splicing events'],
    license='GNU General Public License v3',
    license_files=('LICENSE',),
    description='PTM-POSE: PTM Projection onto Splice Events',
    long_description=long_description,
    long_description_content_type="text/markdown",
    project_urls = {'Homepage':'https://github.com/NaegleLab/PTM-POSE/', 'Issues': 'https://github.com/NaegleLab/PTM-POSE/issues', 'Documentation': 'https://naeglelab.github.io/PTM-POSE/'},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
    ],
    packages=find_packages(),
    include_package_data = True,
    python_requires=">=3.10,<3.13"
)

