================
Running PTM-POSE
================

PTM-POSE is an easily implementable tool to project PTM sites onto splice event data generated from RNA sequencing data and is compatible with any splice event quantification tool that outputs genomic coordinates of different splice events (MATS, SpliceSeq, etc.). PTM-POSE harnesses PTMs that have been mapped to their genomic location by a sister package, [ExonPTMapper](https://github.com/NaegleLab/ExonPTMapper). It also contains functions for annotating these PTMs with information from various databases, like PhosphoSitePlus and ELM.

Formatting Data
---------------

To run PTM-POSE, you first need to process your data such that each row corresponds to a unique splice event with the genomic location of that splice event (chromosome, strand, and the bounds of the spliced region). Strand can be indicated using either '+'/'-' or 1/-1. If desired, you can also provide a delta PSI and significance value which will be included in the final PTM dataframe. Any additional columns will be kept. At a minimum, the dataframe should look something like this (optional but recommended parameters indicated):

+---------------------+-----------------------+------------+--------+--------------+------------+-----------------+-------------------------+
| event id            | Gene name             | chromosome | strand | region start | region end | dPSI            | significance            |
| (optional)          | (recommended)         |            |        |              |            | (optional)      | (optional)              |
+=====================+=======================+============+========+==============+============+=================+=========================+
| first_event         | CSTN1                 |1           |  \-1   | 9797555      | 9797612    | 0.362           | 0.032                   |
+---------------------+-----------------------+------------+--------+--------------+------------+-----------------+-------------------------+


PTM-POSE allows you to assess two potential impacts of splicing on PTMs: 

Differential inclusion 
    lost or gained from the isoform as a result of a splice event 
Altered flanking sequences
    the PTM site is present in both isoforms, but the adjacent residues around a PTM are changed in one isoform such that its linear motif that drives many protein interactions is unique 

Identifying differentially included PTMs
----------------------------------------

Once the data is in the correct format, simply run the `project_ptms_onto_splice_events()` function, indicating the column names corresponding each data element. By default, PTM-POSE assumes the provided coordinates are in hg38 coordinates, but you can use older coordinate systems with the `coordinate_type` parameter. If you have saved ptm_coordinates locally, you can set this parameter to None.

.. code-block:: python

    from ptm-pose import project

    my_splice_data_annotated, spliced_ptms = project.project_ptms_onto_splice_events(my_splice_data, 
            chromosome_col = 'chromosome',
            strand_col = 'strand',
            region_start_col = 'region start',
            region_end_col =  'region end',
            event_id_col = 'event id',
            gene_col = 'Gene name',
            dPSI_col='dPSI',
            coordinate_type = 'hg19')

Altered Flanking Sequences
--------------------------

In addition to the previously mentioned columns, we will need to know the location of the flanking exonic regions next to the spliced region. Make sure your dataframe contains the following information prior to running flanking sequence analysis:

+---------------------+-------------------------+------------+--------+--------------+------------+-------------------+-----------------+--------------------+------------------+-----------------+-------------------------+
| event id            | Gene name               | chromosome | strand | region start | region end | first flank start | first flank end | second flank start | second flank end | dPSI            | significance            |
| (optional)          | (recommended)           |            |        |              |            |                   |                 |                    |                  | (recommended)   | (recommended)           |
+=====================+=========================+============+========+==============+============+===================+=================+====================+==================+=================+=========================+
| first event         |  CSTN1                  | 1          |  \-1   | 9797555      | 9797612    | 9687655           | 9688446         | 9811223            | 9811745          | 0.362           | 0.032                   |
+---------------------+-------------------------+------------+--------+--------------+------------+-------------------+-----------------+--------------------+------------------+-----------------+-------------------------+

Then, as with differentially included PTMs, you only need to run `get_flanking_changes_from_splice_data()` function:

.. code-block:: python

    from ptm-pose import project

    altered_flanks = project.get_flanking_changes_from_splice_data(my_splice_data, 
            chromosome_col = 'chromosome',
            strand_col = 'strand',
            region_start_col = 'region_start',
            region_end_col =  'region_end',
            first_flank_start_col = 'first_flank_start',
            first_flank_end_col = 'first_flank_end',
            second_flank_start_col = 'second_flank_start',
            second_flank_end_col = 'second_flank_start',
            event_id_col = 'event_id',
            gene_col = 'Gene name',
            dPSI_col='dPSI',
            coordinate_type = 'hg19')



Using PTM-POSE with MATS or SpliceSeq datasets
-----------------------------------------------

If you have quantified your splicing events using MATS or SpliceSeq, we have provided custom functions for processing and projecting PTMs onto these datasets.

MATS datasets
~~~~~~~~~~~~~

MATS outputs several different tab-separated files for five different splice events (SE, MXE, A5SS, A3SS, RI). Each of these contain different columns so must be processed separately and then combined into the final dataframes. We have provided a function `project_ptms_onto_MATS()` that will process each of the different dataframes and combine them into a single dataframe. You can choose whether or not you would like to identify flanking sequences as well by setting the `identify_flanking_sequences` parameter to True or False. To run this function, first load each event-specific dataset and then run the function:

.. code-block:: python

    from ptm-pose import project

    SE_data = pd.read_csv('SE.MATS.JunctionCountOnly.txt', sep = '\t')
    MXE_data = pd.read_csv('MXE.MATS.JunctionCountOnly.txt', sep = '\t')
    A5SS_data = pd.read_csv('A5SS.MATS.JunctionCountOnly.txt', sep = '\t')
    A3SS_data = pd.read_csv('A3SS.MATS.JunctionCountOnly.txt', sep = '\t')
    RI_data = pd.read_csv('RI.MATS.JunctionCountOnly.txt', sep = '\t')

    splice_data, spliced_ptms, altered_flanks = project.project_ptms_onto_MATS(SE_events = SE_data, MXE_events = MXE_data, A5SS_events = A5SS_data, A3SS_events = A3SS_data, RI_events = RI_data, 
            coordinate_type = 'hg19',
            identify_flanking_sequences = True)

If `identify_flanking_sequences = True`, the function will return three outputs, the original annotated splice data, the differentially included PTMs, and the altered flanking sequences. If `identify_flanking_sequences = False`, the function will only return the original annotated splice data and the differentially included PTMs.

.. note::
    The `project_ptms_onto_MATS()` function assumes that the columns in the MATS files are the same described in the MATS documentation, but if using older versions of MATS or customized outputs, you may need to adjust the column names, such as those for dPSI and significance.

SpliceSeq datasets
~~~~~~~~~~~~~~~~~~

SpliceSeq outputs a single tab-separated file with quantification for all events that are associated with the SpliceSeq splicegraph. The output of SpliceSeq only outputs event information relative to SpliceSeq exons, so you will also need to provide the SpliceSeq splicegraph file. We have provided a function `project_ptms_onto_SpliceSeq()` that will process the SpliceSeq data and the splicegraph and return the differentially included PTMs and altered flanking sequences. To run this function, first load the SpliceSeq data and the splicegraph file. 

.. code-block:: python

    from ptm-pose import project

    SpliceSeq_data = pd.read_csv('SpliceSeq_output.txt', sep = '\t')
    splicegraph = pd.read_csv('SpliceSeq_splicegraph.txt', sep = '\t')

    spliced_ptms, altered_flanks = project.project_ptms_onto_SpliceSeq(psi_data = SpliceSeq_data, splicegraph = splicegraph, dPSI_col = 'dPSI', sig_col = 'Significance',
            coordinate_type = 'hg19')

.. note::
    To have dPSI or significance columns in the final output, you must provide the column names in the `dPSI_col` and `sig_col` parameters. If you do not have these columns, you can omit them and the final output will not include them.