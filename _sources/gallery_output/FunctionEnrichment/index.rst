

.. _sphx_glr_gallery_output_FunctionEnrichment:

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


.. toctree::
   :hidden:

   /gallery_output/FunctionEnrichment/plot_num_annotations
   /gallery_output/FunctionEnrichment/plot_geneset_enrichment
   /gallery_output/FunctionEnrichment/plot_specific_annotation

