from setuptools import setup, find_packages


# Run setup
setup(
    name="ptm_pose",
    version="0.1.0",
    author="Naegle Lab",
    author_email="kmn4mj@virginia.edu",
    url="https://github.com/NaegleLab/PTM-POSE/tree/main",
    install_requires=['pandas==2.2.*', 'numpy==1.26.*', 'scipy==1.13.*', 'biopython==1.83.*', 'tqdm==4.66.*', 'xlrd', 'matplotlib', 'requests'],
    license='GNU General Public License v3',
    description='PTM-POSE: PTM Projection onto Splice Events',
    long_description="""PTM-POSE is an open-source tool for annotating splice event quantification data with post-translational modifications (PTMs) and their functional consequences.""",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
    ],
    packages=find_packages(),
    include_package_data = True,
    python_requires=">=3.9"
)

