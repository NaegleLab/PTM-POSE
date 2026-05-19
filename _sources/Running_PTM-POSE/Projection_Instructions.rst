==================================
Projecting PTMs onto Splice Events
==================================

Once you have PTM-POSE installed, you simply need to make sure your splice event data is in a compatible format, indicate which columns contain the necessary information, and then let PTM-POSE do the rest of the work. PTM-POSE allows you to assess two main potential impacts of splicing to PTMs: 

* Differential inclusion 
    lost or gained from the isoform as a result of a splice event 
* Altered flanking sequences
    the PTM site is present in both isoforms, but the adjacent residues around a PTM are changed in one isoform such that its linear motif that drives many protein interactions is unique 
    
Follow this tutorial to learn how to run PTM-POSE on your data, depending on the tool you used to quantify your splicing events. 

Formatting Data
---------------

Tools with Built-in Compatibility
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Depending on the tool you used to quantify your splicing events, you may be able to use one of PTM-POSE's built-in modules for processing and projecting PTMs onto your data. The following tools currently have built-in compatibility with PTM-POSE:
* MATS
* SpliceSeq
* MAJIQ (differential inclusion only)

If you have data from one of these tools, you can skip the formatting step and go straight to running the built-in functions for these tools, which will format your data for you. See the sections below for instructions on how to run PTM-POSE with these tools.


Other Quantification Tools
~~~~~~~~~~~~~~~~~~~~~~~~~~

If you have data from a tool other than those listed above, you can still use PTM-POSE, but you will need to make sure your data contains the information expected by PTM-POSE. To run PTM-POSE, your data should have the following information:

- Separate rows for each splice event
- Columns with the genomic coordinates of the splice event (chromosome, strand, and the bounds of the spliced region)
- A unique identifier for each splice event, which allows for tracking the splice event across different analyses and outputs
- (Optional but recommended) Columns with the gene name, delta PSI, and significance value for each splice event
- (For flanking sequence analysis) Columns with the coordinates of the flanking exonic regions next to the spliced region



At a minimum, the dataframe should look something like this (optional but recommended parameters indicated):

+---------------------+-----------------------+------------+--------+--------------+------------+-----------------+-------------------------+
| event id            | Gene name             | chromosome | strand | region start | region end | dPSI            | significance            |
| (optional)          | (recommended)         |            |        |              |            | (optional)      | (optional)              |
+=====================+=======================+============+========+==============+============+=================+=========================+
| first_event         | CSTN1                 |1           |  \-1   | 9797555      | 9797612    | 0.362           | 0.032                   |
+---------------------+-----------------------+------------+--------+--------------+------------+-----------------+-------------------------+

.. note::
    The strand information can be denoted either using integers (1 for positive strand, -1 for negative strand) or using characters ('+' for positive strand, '-' for negative strand).  

Should you wish to identify PTMs with altered flanking sequences, the dataframe should also contain the following columns with the coordinates of the flanking exonic regions next to the spliced region:

+---------------------+-------------------------+------------+--------+--------------+------------+-------------------+-----------------+--------------------+------------------+-----------------+-------------------------+
| event id            | Gene name               | chromosome | strand | region start | region end | first flank start | first flank end | second flank start | second flank end | dPSI            | significance            |
| (optional)          | (recommended)           |            |        |              |            |                   |                 |                    |                  | (recommended)   | (recommended)           |
+=====================+=========================+============+========+==============+============+===================+=================+====================+==================+=================+=========================+
| first event         |  CSTN1                  | 1          |  \-1   | 9797555      | 9797612    | 9687655           | 9688446         | 9811223            | 9811745          | 0.362           | 0.032                   |
+---------------------+-------------------------+------------+--------+--------------+------------+-------------------+-----------------+--------------------+------------------+-----------------+-------------------------+


Running PTM-POSE
-----------------

Regardless of the tool you used to quantify your splice events, the PTM-POSE workflow will be the same:

1. Load your data into pandas dataframe(s)
2. Initialize the appropriate class for your dataset, indicating any filters for dPSI/significance
3. Run PTM-POSE with the `run_pose()` function. You can choose whether or not to identify altered flanking sequences as well by setting the `identify_flanking_sequences` parameter to True or False.
4. (optional) Run NEASE analysis on the same splice event data with the `run_nease()` function, which provides complementary but non-overlapping information about the impact of splicing on protein interactions and pathways (domains, linear motifs, etc.)
5. Save the outputs to csv files for downstream analysis

See the following sections for instructions on how to run PTM-POSE with different types of datasets, including those from tools with built-in compatibility and those from other tools.


MATS Datasets
~~~~~~~~~~~~~

MATS outputs several different tab-separated files for five different splice events (SE, MXE, A5SS, A3SS, RI). Each of these contain different columns so must be processed separately and then combined into the final dataframes. The :class:`MATS_Dataset <ptm_pose.splicing_tools.MATS.MATS_Dataset>` class can be used to load and process these datasets. You can choose to process as many of the event files as you want.

.. code-block:: python

    from ptm-pose.splicing_tools.MATS import MATS_Dataset
    import pandas as pd

    #load MATS files
    SE = pd.read_csv('SE.MATS.JCEC.txt', sep = '\t')
    A3SS = pd.read_csv('A3SS.MATS.JCEC.txt', sep = '\t')
    A5SS = pd.read_csv('A5SS.MATS.JCEC.txt', sep = '\t')
    MXE = pd.read_csv('MXE.MATS.JCEC.txt', sep = '\t')
    RI = pd.read_csv('RI.MATS.JCEC.txt', sep = '\t')

    #initialize the MATS_Dataset object with the dataframes for each event type, specifying the coordinate type and column names for dPSI and significance
    mats_dataset = MATS_Dataset(SE = SE, A3SS = A3SS, A5SS = A5SS, MXE = MXE, RI = RI, coordinate_type = 'hg19', dpsi_col = 'meanDeltaPSI', sig_col = 'FDR', min_dpsi = 0.2, alpha = 0.05)

    #run PTM-POSE, including flanking sequence analysis
    mats_dataset.run_pose(identify_flanking_sequences = True)

    #run NEASE analysis on the same dataset
    mats_dataset.run_nease()

    mats_dataset.save_outputs(odir = './output_dir/')


SpliceSeq Datasets
~~~~~~~~~~~~~~~~~~

SpliceSeq outputs a single tab-separated file with quantification for all events that are associated with the SpliceSeq splicegraph. The output of SpliceSeq only outputs event information relative to SpliceSeq exons, so you will also need to provide the SpliceSeq splicegraph file. The :class:`SpliceSeq_Dataset <ptm_pose.splicing_tools.SpliceSeq.SpliceSeq_Dataset>` class can be used to load and process these datasets.

.. code-block:: python

    from ptm-pose.splicing_tools.SpliceSeq import SpliceSeq_Dataset
    import pandas as pd

    #load SpliceSeq data and splicegraph files
    SpliceSeq_data = pd.read_csv('SpliceSeq_output.txt', sep = '\t')
    splicegraph = pd.read_csv('SpliceSeq_splicegraph.txt', sep = '\t')

    #initialize the SpliceSeq_Dataset object with the dataframes for the SpliceSeq output and splicegraph, specifying the coordinate type and column names for dPSI and significance
    spliceseq_dataset = SpliceSeq_Dataset(data = SpliceSeq_data, splicegraph = splicegraph, coordinate_type = 'hg19')

    #run PTM-POSE, including flanking sequence analysis
    spliceseq_dataset.run_pose(identify_flanking_sequences = True)

    #run NEASE analysis on the same dataset
    spliceseq_dataset.run_nease()

    spliceseq_dataset.save_outputs(odir = './output_dir/')


MAJIQ Datasets (voila tsv files)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When using MAJIQ for splice event quantification and differential comparison, it outputs a tab-separated file with delta PSI and significance values for different LSVs (also called voila files). The :class:`MAJIQ_Dataset <ptm_pose.splicing_tools.MAJIQ.MAJIQ_Dataset>` class can be used to load and process these datasets.

.. code-block:: python

    from ptm-pose.splicing_tools.MAJIQ import MAJIQ_Dataset
    import pandas as pd

    #load MAJIQ voila file and indicate the sample names used
    voila_tsv_file = 'voila_output.tsv'
    samp1_name = 'sample1'
    samp2_name = 'sample2'

    #initialize the MAJIQ_Dataset object with the voila dataframe, specifying the coordinate type and column names for dPSI and significance
    majiq_dataset = MAJIQ_Dataset(voila_tsv_file = voila_tsv_file, coordinate_type = 'hg19', dpsi_col = 'delta_PSI', sig_col = 'pval', min_dpsi = 0.2, alpha = 0.05)

    #run PTM-POSE, including flanking sequence analysis
    majiq_dataset.run_pose(identify_flanking_sequences = True)

    #run NEASE analysis on the same dataset
    majiq_dataset.run_nease()

    majiq_dataset.save_outputs(odir = './output_dir/')





Other Datasets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For any splice event data other those described above, you can use the 'base.GenericDataset' class.

When using MAJIQ for splice event quantification and differential comparison, it outputs a tab-separated file with delta PSI and significance values for different LSVs (also called voila files). The :class:`GenericDataset <ptm_pose.splicing_tools.base.GenericDataset>` class can be used to load and process these datasets.


To run PTM-POSE with this class, first load your data into a pandas dataframe (or dictionary of dataframes for multiple event types) and then create a GenericDataset object with the appropriate column names. Then, you can run the :func:`run_pose_generic() <ptm_pose.splicing_tools.base.GenericDataset.run_pose_generic>` function within the GenericDataset object to project PTMs onto your splice events. If you wish to identify altered flanking sequences as well, add the `identify_flanking_sequences = True` parameter to the function. See below for an example of how to run PTM-POSE with a generic dataset:

.. code-block:: python

    import pandas as pd
    from ptm-pose.splicing_tools.base import GenericDataset

    my_splice_data = pd.read_csv('my_splice_data.csv')

    #initialize the object specifying the column names containing splice event information
    my_dataset = GenericDataset(my_splice_data, 
            chromosome_col = 'chromosome',
            strand_col = 'strand',
            region_start_col = 'region start',
            region_end_col =  'region end',
            event_id_col = 'event id',
            gene_col = 'Gene name',
            dpsi_col='dPSI',
            sig_col = 'significance',
            coordinate_type = 'hg19')

    my_dataset.run_pose_generic(identify_flanking_sequences = True, extra_cols = None)

    my_dataset.run_nease_generic()

    my_dataset.save_outputs(odir = './output_dir/')


You may need to provide additional information depending on your specific use case. See below for additional options you may wish to adjust when running PTM-POSE with a generic dataset.

+-------------------------+------------------------------------------------------------------------------------------+---------------+
| Parameter               | Description                                                                              | Default Value |
+=========================+==========================================================================================+===============+
| min_dpsi                | Minimum delta PSI required to be considered for analysis (ignored if dpsi_col is None)   | 0             |  
+-------------------------+------------------------------------------------------------------------------------------+---------------+
| alpha                   | Significance threshold for analysis (ignored if sig_col is None)                         | 0.05          |
+-------------------------+------------------------------------------------------------------------------------------+---------------+
| start_coordinate_system | Whether start coordinate is 0- or 1-based (MATS is 0-based, for example)                 | '1-based'     |
+-------------------------+------------------------------------------------------------------------------------------+---------------+


Downstream Analysis
-------------------

Once you have run PTM-POSE and generated the output dataframes, you can use the various analysis modules described in the analysis gallery to explore the impacts of splicing on PTMs in your dataset. You can also use the `annotate` module to append additional information about your PTMs from various databases, which may be useful for some of the analysis modules. See the following section for instructions on how to use the `annotate` module to append additional information about your PTMs.