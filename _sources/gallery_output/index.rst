:orphan:

Gallery
=======

Below you will find different ways you might choose to analyze the PTMs identified by PTM-POSE:


.. raw:: html

    <div class="sphx-glr-thumbnails">

.. thumbnail-parent-div-open

.. thumbnail-parent-div-close

.. raw:: html

    </div>

General Overview of Spliced PTMs
---------------------------------



.. raw:: html

    <div class="sphx-glr-thumbnails">

.. thumbnail-parent-div-open

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="One of the first things you will want to do when analyzing spliced PTMs is to check how many and what type of PTMs were identified as being differentially included or altered due to splicing events. This can be done easily using the plotting module in PTM-POSE.">

.. only:: html

  .. image:: /gallery_output/A_Overview_of_Spliced_PTMs/images/thumb/sphx_glr_plot_NumberOfPTMs_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_output_A_Overview_of_Spliced_PTMs_plot_NumberOfPTMs.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Inspecting types of PTMs impacted by splicing</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="As described in Running PTM-POSE section, PTM-POSE provides various options for annotating functional information for PTMs, coming from various databases. However, PTM functional information is inherently sparse, and so most annotations will only provide information on a handful of PTMs. For this reason, it can be useful to probe how many PTMsTo better understand the types of annotations that are available, as well as the number of PTMs that have an annotation of that type. This can be done using the analyze function in PTM-POSE.">

.. only:: html

  .. image:: /gallery_output/A_Overview_of_Spliced_PTMs/images/thumb/sphx_glr_plot_PTM_annotations_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_output_A_Overview_of_Spliced_PTMs_plot_PTM_annotations.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Inspecting number of PTMs with annotation information available</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="In some cases, you may wish to filter out PTMs with weaker evidence (only identified in small number of experiments), splice events with small changes, or PTMs that have not previously been observed to be impacted by splicing. This can be done using the filter_ptms function in the helpers module. We can also visualize the impact of filtering using the analyze.filter module. ">

.. only:: html

  .. image:: /gallery_output/A_Overview_of_Spliced_PTMs/images/thumb/sphx_glr_plot_filter_impact_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_output_A_Overview_of_Spliced_PTMs_plot_filter_impact.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Impact of filtering PTMs</div>
    </div>


.. thumbnail-parent-div-close

.. raw:: html

    </div>

Functional Impact of Spliced PTMs
---------------------------------



.. raw:: html

    <div class="sphx-glr-thumbnails">

.. thumbnail-parent-div-open

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="As described in Running PTM-POSE section, PTM-POSE provides various options for annotating functional information for PTMs, coming from various databases. However, PTM functional information is inherently sparse, and so most annotations will only provide information on a handful of PTMs. For this reason, it can be useful to probe how many PTMsTo better understand the types of annotations that are available, as well as the number of PTMs that have an annotation of that type. This can be done using the analyze function in PTM-POSE.">

.. only:: html

  .. image:: /gallery_output/FunctionEnrichment/images/thumb/sphx_glr_plot_num_annotations_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_output_FunctionEnrichment_plot_num_annotations.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Identifying available PTM-specific annotation information</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="As is commonly done for exon-centric analyses, we have provided the ability to perform gene set enrichment analysis for gene associated with spliced PTMs, using the EnrichR API from the gseapy module. By default, we include gene ontology terms, KEGG pathways, and Reactome pathways, but you can also provide your own gene sets listed in EnrichR.">

.. only:: html

  .. image:: /gallery_output/FunctionEnrichment/images/thumb/sphx_glr_plot_geneset_enrichment_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_output_FunctionEnrichment_plot_geneset_enrichment.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Identify enriched gene sets associated with differentially spliced PTMs</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="As described in Running PTM-POSE section, PTM-POSE provides various options for annotating functional information for PTMs, coming from various databases. Often, we will want to dig deeper into the specific functions, processes, interactions, etc. associated with the proteins in our dataset. First, we can look at the annotations currently available for analysis, based on annotations that have been appended using the annotate module:">

.. only:: html

  .. image:: /gallery_output/FunctionEnrichment/images/thumb/sphx_glr_plot_specific_annotation_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_output_FunctionEnrichment_plot_specific_annotation.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Assessing enriched PTM functions</div>
    </div>


.. thumbnail-parent-div-close

.. raw:: html

    </div>

Impact on protein interaction and regulatory networks
-----------------------------------------------------




.. raw:: html

    <div class="sphx-glr-thumbnails">

.. thumbnail-parent-div-open

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="In the realm of kinase activity analysis, Kinase Substrate Enrichment Analysis (KSEA) is a commonly used tool to identify kinases with substrates that are significantly changing between two conditions (based on phosphoproteomic data). See the publication here  for what the original algorithm looks like. We have adapted KSEA for use with splicing data, which uses the dPSI values to calculate a z-score to identify kinases with substrates that are either excluded or included from isoforms. This is a useful tool, but you should treat these results with caution as they will often be on only a couple kinase-substrates due to the limited kinase-substrate information available. First, let&#x27;s load the example data and run KSEA on the differentially included PTMs dataframe">

.. only:: html

  .. image:: /gallery_output/Networks/images/thumb/sphx_glr_plot_ksea_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_output_Networks_plot_ksea.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Kinase substrate enrichment analysis (KSEA)</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Post translational modifications (PTMs) often facilitate protein interactions, either through direct binding of domains specific to that particular modification (e.g. SH2 domains binding to phosphorylated tyrosines) or through allosteric effects that change the conformation of the protein to either enhance or disrupt interactions. We provide functions to annotate spliced PTMs with relevant protein interactions and to identify key PTMs that may disrupt protein interaction networks.">

.. only:: html

  .. image:: /gallery_output/Networks/images/thumb/sphx_glr_plot_protein_interactions_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_output_Networks_plot_protein_interactions.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Protein interactions network</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Given that phosphorlaiton are one of the most commonly impacted modifications, there is potential for kinases targeting these sites to be indirectly impacted by alternative splicing through changes in the availability of their substrates. While we provide functions for performing enrichment of known kinase substrates from databases like PhosphoSitePlus, RegPhos, and PTMsigDB, these resources are limited by the overall number of validated substrates (&lt;5%). For this purpose, we have adapted a previously developed algorithm called KSTAR (Kinase Substrate to Activity Relationships) for use with spliced PTM data, which harnesses kinase-substrate predictions to expand the overall number of phosphorylation sites that can be used as evidence. This particularly important as you may find many of the spliced PTMs in your dataset are less well studied and may not have any annotated kinases.">

.. only:: html

  .. image:: /gallery_output/Networks/images/thumb/sphx_glr_plot_kstar_enrichment_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_output_Networks_plot_kstar_enrichment.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Identify kinases with enriched substrates in differentially included exons, using an adapted version of KSTAR</div>
    </div>


.. thumbnail-parent-div-close

.. raw:: html

    </div>

Analyzing altered flanking sequences
------------------------------------




.. raw:: html

    <div class="sphx-glr-thumbnails">

.. thumbnail-parent-div-open

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Given a PTM identified as having an altered flanking sequence in the alternative isoform and/or due to a splice event, it can be usefult o visualize how the sequences differ. We provide functions to compare the flanking sequences of PTMs based on whether an exonic region is included or excluded">

.. only:: html

  .. image:: /gallery_output/flanking_sequences/images/thumb/sphx_glr_plot_sequence_comparison_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_output_flanking_sequences_plot_sequence_comparison.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Inspect difference of flanking sequence due to splice event</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="A potential consequence of altered flanking sequences is that the short linear motifs important for driving PTM-specific domain interactions are disrupted, such as for the phosphotyrosine binding domains SH2 or 14-3-3 proteins (which bind to phosphoserines and threonines). Using PTM-POSE and motif data from [ELM](http://elm.eu.org/searchdb.html), we can identify and visualize the altered 14-3-3 domain motifs due to splicing events.">

.. only:: html

  .. image:: /gallery_output/flanking_sequences/images/thumb/sphx_glr_plot_sh2_domain_motifs_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_output_flanking_sequences_plot_sh2_domain_motifs.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Identify altered SH2 domain motifs</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="In order to understand how PTMs may be altered due to splicing events, it is useful to identify the flanking sequences of the PTMs and how they may be altered due to nearby splice events (as identified by flanking sequence module). Once we have, this information we can analyze and visualize where the alterations in the flanking sequences occur. First, we need to compare the flanking sequences of PTMs based on whether an exonic region is included or excluded using the compare_flanking_sequences function in PTM-POSE.">

.. only:: html

  .. image:: /gallery_output/flanking_sequences/images/thumb/sphx_glr_plot_location_altered_flanks_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_output_flanking_sequences_plot_location_altered_flanks.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Probing where and how PTM flanking sequences are altered</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Kinases are known to recognize specific motifs surrounding Y, S, and T residues, and these motifs can be used to predict kinase interactions using a tool called Kinase Library. Here, we will use Kinase Library to identify kinases who may have a stronger preference for the flanking sequences of PTMs in one isoform vs. another. For this example, let&#x27;s focus on phosphorylation sites that have larger changes &gt;0.4 in the splicing data.">

.. only:: html

  .. image:: /gallery_output/flanking_sequences/images/thumb/sphx_glr_plot_kinase_library_affinity_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_output_flanking_sequences_plot_kinase_library_affinity.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Kinase affinity due to PTM flanking sequence alterations</div>
    </div>


.. thumbnail-parent-div-close

.. raw:: html

    </div>


.. toctree::
   :hidden:
   :includehidden:


   /./gallery_output/A_Overview_of_Spliced_PTMs/index.rst
   /./gallery_output/FunctionEnrichment/index.rst
   /./gallery_output/Networks/index.rst
   /./gallery_output/flanking_sequences/index.rst


.. only:: html

  .. container:: sphx-glr-footer sphx-glr-footer-gallery

    .. container:: sphx-glr-download sphx-glr-download-python

      :download:`Download all examples in Python source code: gallery_output_python.zip </gallery_output/gallery_output_python.zip>`

    .. container:: sphx-glr-download sphx-glr-download-jupyter

      :download:`Download all examples in Jupyter notebooks: gallery_output_jupyter.zip </gallery_output/gallery_output_jupyter.zip>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
