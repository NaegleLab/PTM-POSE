# PTM-POSE (PTM Projection Onto Splice Events)

PTM-POSE is an easily implementable tool to project PTM sites onto splice event data generated from RNA sequencing data and is compatible with any splice event quantification tool that outputs genomic coordinates of different splice events (MATS, SpliceSeq, etc.). PTM-POSE harnesses PTMs that have been mapped to their genomic location by a sister package, [ExonPTMapper](https://github.com/NaegleLab/ExonPTMapper). It also contains functions for annotating these PTMs with information from various databases, like PhosphoSitePlus and ELM.

## Running PTM-POSE

To run PTM-POSE, you first need to process your data such that each row corresponds to a unique splice event with the genomic location of that splice event (chromosome, strand, and the bounds of the spliced region). Strand can be indicated using either '+'/'-' or 1/-1. If desired, you can also provide a delta PSI and significance value which will be included in the final PTM dataframe. Any additional columns will be kept. At a minimum, the dataframe should look something like this (optional but recommended parameters indicated):
| event_id (optional) | chromosome | strand | region_start | region_end | dPSI (optional) | significance (optional) |
|---------------------|------------|--------|--------------|------------|-----------------|-------------------------|
| first_event         | 1          |  -     | 9797555      | 9797612    | 0.362           | 0.032                   |

 Once the data is in the correct format, simply run the project_ptms_onto_splice_events() function. By default, PTM-POSE assumes the provided coordinates are in hg38 coordinates, but you can use older coordinate systems with the `coordinate_type` parameter.
```python
from ptm-pose import project

my_splice_data_annotated, spliced_ptms = project.project_ptms_onto_splice_events(my_splice_data, ptm_coordinates,
                                                                                  chromosome_col = 'chromosome',
                                                                                  strand_col = 'strand',
                                                                                  region_start_col = 'region_start',
                                                                                  region_end_col =  'region_end',
                                                                                  event_id_col = 'event_id',
                                                                                  dPSI_col='dPSI',
                                                                                  coordinate_type = 'hg19')
```

This will produce two dataframes:
1. Original splice data with additional columns indicating the number and which PTMs were found associated with that splice event. 'PTMs column denotes the UniProtKB accession, residue, site number, and modification type for PTM identified.
   
| event_id (if provided) | chromosome | strand | region_start | region_end | dPSI (if provided) | significance (if provided) | PTMs                          | Number of PTMs Affected |
|---------------------|------------|--------|--------------|------------|-----------------|-------------------------|-------------------------------|-------------------------|
| first_event         | 1          |  -     | 9797555      | 9797612    | 0.362           | 0.032                   | O94985_N515 (N-Glycosylation) | 1                       |

2. New dataframe where each row is a unique event-PTM pair. This is useful for downstream analyses of the important PTM changes that are occuring in your dataset, and many functions provided for further annotation and analyses of these PTMs (see rest of documentation for examples)
   
| event_id (if provided) | UniProtKB Accession | Residue | Modifications | PTM Info | dPSI (if provided) | significance (if provided) |
|----------|---------------------|---------|---------------|----------|--------------------|----------------------------|
| first_event | O94985 | N515 |  N-Glycosylation | O94985_N515 (N-Glycosylation) | 0.362 | 0.032 |

Alternatively, we have also provided tool-specific functions for MATS and SpliceSeq outputs, which are tailored to the specific output of these tools. These functions are `project_ptms_onto_MATS()` and `project_ptms_onto_SpliceSeq()`, respectively.

## Annotating and Analyzing PTMs

In addition to projection of PTMs, PTM-POSE offers a variety of functions for annotating PTMs with functional information from various databases including PhosphoSitePlus, RegPhos, PTMcode, PTMsigDB, and ELM. These functions are useful for understanding the functional implications of PTMs in your dataset. We also provide functions for identifying PTMs and functions that are enriched in your dataset, including enrichment of kinase-substrate relationships. 

For more information on this analysis, see the full documentation.

