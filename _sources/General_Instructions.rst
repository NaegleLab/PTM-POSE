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
            ptm_coordinates,
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
            ptm_coordinates,
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

Combining outputs
-----------------
In some cases you may wish to work with a combined file that indicates both differential inclusion and altered flanking sequence events. This can be done quickly by running:

.. code-block:: python

    from ptm_pose import analyze
    combined_output = analyze.combine_outputs(spliced_ptms, altered_flanks)

Annotating PTMs with Functional Information
-------------------------------------------
Beyond projecting PTMs onto your data, we have also provided additional functions for appending information on the function, relationships, and interactions of each post-translational modification that have been recorded in various databases. These annotations include information from:

+---------------------------------------------------------------------+------------------------+--------------------------------------------------------------------------------------------------------------+
| Database                                                            |  Annotation types      | PTM-POSE function                                                                                            | 
+=====================================================================+========================+==============================================================================================================+
| `PhosphoSitePlus <https://www.phosphosite.org/homeAction.action>`_  |- Function              |.. code-block:: python                                                                                        | 
|                                                                     |- Biological Process    |                                                                                                              |    
|                                                                     |- interactions          |    annotate.add_PSP_regulatory_site_data(spliced_ptms, file = "/path/to/file/Regulatory_sites.gz")           | 
|                                                                     +------------------------+--------------------------------------------------------------------------------------------------------------+
|                                                                     |- Kinase substrates     |.. code-block:: python                                                                                        |           
|                                                                     |                        |                                                                                                              |               
|                                                                     |                        |    annotate.add_PSP_kinase_substrate_data(spliced_ptms, file = "/path/to/file/Kinase_Substrate_Dataset.gz"   |             
+---------------------------------------------------------------------+------------------------+--------------------------------------------------------------------------------------------------------------+
| `DEPOD <https://depod.bioss.uni-freiburg.de/>`_                     |- Phosphatase substrates|.. code-block:: python                                                                                        |
|                                                                     |                        |                                                                                                              |
|                                                                     |                        |    annotate.add_DEPOD_data(spliced_ptms, file = "/path/to/file/")                                            | 
+---------------------------------------------------------------------+------------------------+--------------------------------------------------------------------------------------------------------------+
| `RegPhos <http://140.138.144.141/~RegPhos/index.php>`_              |- Kinase substrates     |.. code-block:: python                                                                                        |
|                                                                     |                        |                                                                                                              |
|                                                                     |                        |    annotate.add_RegPhos_data(spliced_ptms, file = "/path/to/file/")                                          | 
+---------------------------------------------------------------------+------------------------+--------------------------------------------------------------------------------------------------------------+
| `ELM <http://elm.eu.org/>`_                                         |- Interactions          |.. code-block:: python                                                                                        |
|                                                                     |                        |                                                                                                              |
|                                                                     |                        |    annotate.add_PTMcode_interprotein(spliced_ptms, file = "/path/to/file/")                                  |
|                                                                     +------------------------+--------------------------------------------------------------------------------------------------------------+ 
|                                                                     |- Linear motifs         |.. code-block:: python                                                                                        |
|                                                                     |                        |                                                                                                              |
|                                                                     |                        |    annotate.add_PTMcode_intraprotein(spliced_ptms, file = "/path/to/file/")                                  |
+---------------------------------------------------------------------+------------------------+--------------------------------------------------------------------------------------------------------------+
| `PTMcode <https://ptmcode.embl.de/>`_                               |- Interactions          |.. code-block:: python                                                                                        |
|                                                                     |                        |                                                                                                              |
|                                                                     |                        |    annotate.add_PTMcode_interprotein(spliced_ptms, file = "/path/to/file/")                                  |
|                                                                     +------------------------+--------------------------------------------------------------------------------------------------------------+
|                                                                     |- Intraprotein contacts |.. code-block:: python                                                                                        |
|                                                                     |                        |                                                                                                              |
|                                                                     |                        |   annotate.add_PTMcode_intraprotein(spliced_ptms, file = "/path/to/file/")                                   |
+---------------------------------------------------------------------+------------------------+--------------------------------------------------------------------------------------------------------------+





Rather than running each function individually, you can also use the master function `annotate_ptms()` to annotate with all desired information at once.

We are continuing to work on adding functions to append more contextual information for individual PTMs. If you have suggestions for what information you would like to be added, please let us know!

Downstream Analysis
-------------------

PTM-POSE also provides functions in the `annotate` module for annotating the above outputs with functional information from various databases: PhosphoSitePlus, RegPhos, PTMcode, PTMInt, ELM, DEPOD. You can then identify PTMs with specific functions, interaction, etc. with the `analyze` module. See an example on a real dataset [here](Examples/ESRP1_knockdown).


