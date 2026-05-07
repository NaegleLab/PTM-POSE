What is PTM-POSE?
=================

PTM-POSE is an easily implementable python tool to project PTM sites onto splice event data generated from RNA sequencing data. It also contains functions for annotating these PTMs with information from various databases, like PhosphoSitePlus and ELM. This allows you to explore questions like:
* What events impact PTMs, and what might this mean for the function of the protein?
* How does splicing impact protein interactions facilated by PTMs?
* Are splicing alterations coordinated around specific PTM enzymes, pathways, or functional classes of PTMs?


PTM-POSE allows you to assess two main potential impacts of splicing to PTMs: 

* Differential inclusion 
    lost or gained from the isoform as a result of a splice event 
* Altered flanking sequences
    the PTM site is present in both isoforms, but the adjacent residues around a PTM are changed in one isoform such that its linear motif that drives many protein interactions is unique 


Compatibility with Splice Event Quantification Tools
----------------------------------------------------

It is compatible with any splice event quantification tool that outputs genomic coordinates of different splice events, although it has largely been tested with events-based tools like MATS.

PTM-POSE also comes with built-in tools for analyzing data from the following tools:

* MATS
* MAJIQ
* SpliceSeq (such as from TCGA SpliceSeq database)

If you have a tool you would like to add built-in compatibility for, you can either submit an issue or, better yet, submit a pull request with the code to do so! 


More Information
----------------

For more details about PTM projection and how it can be used to understand the impacts of splicing on cell signaling and other processes, see the published manuscript: `PTM-POSE Paper`_. For details on how genomic information for PTM sites was generated, see PTM-POSE's sister package, `ExonPTMapper`_.


.. _ExonPTMapper: https://github.com/NaegleLab/ExonPTMapper
.. _PTM-POSE Paper: https://www.nature.com/articles/s41467-022-32017-5

