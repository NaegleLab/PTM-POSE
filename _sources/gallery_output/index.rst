:orphan:

Types of Analysis Performed with PTM-POSE
=========================================

Below you will find different ways you might choose to analyze the PTMs identified by PTM-POSE:


.. raw:: html

    <div class="sphx-glr-thumbnails">

.. thumbnail-parent-div-open

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Post translational modifications (PTMs) often facilitate protein interactions, either through direct binding of domains specific to that particular modification (e.g. SH2 domains binding to phosphorylated tyrosines) or through allosteric effects that change the conformation of the protein to either enhance or disrupt interactions. We provide functions to annotate spliced PTMs with relevant protein interactions and to identify key PTMs that may disrupt protein interaction networks.">

.. only:: html

  .. image:: /gallery_output/images/thumb/sphx_glr_plot_protein_interactions_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_output_plot_protein_interactions.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Identify protein interactions that may be impacted by splicing of PTMs</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="As described in Running PTM-POSE section, PTM-POSE provides various options for annotating functional information for PTMs, coming from various databases. However, PTM functional information is inherently sparse, and so most annotations will only provide information on a handful of PTMs. For this reason, it can be useful to probe how many PTMsTo better understand the types of annotations that are available, as well as the number of PTMs that have an annotation of that type. This can be done using the analyze function in PTM-POSE.">

.. only:: html

  .. image:: /gallery_output/images/thumb/sphx_glr_plot_num_annotations_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_output_plot_num_annotations.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Inspecting number of PTMs with annotation information available</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Given that phosphorlaiton are one of the most commonly impacted modifications, there is potential for kinases targeting these sites to be indirectly impacted by alternative splicing through changes in the availability of their substrates. While we provide functions for performing enrichment of known kinase substrates from databases like PhosphoSitePlus, RegPhos, and PTMsigDB, these resources are limited by the overall number of validated substrates (&lt;5%). For this purpose, we have adapted a previously developed algorithm called KSTAR (Kinase Substrate to Activity Relationships) for use with spliced PTM data, which harnesses kinase-substrate predictions to expand the overall number of phosphorylation sites that can be used as evidence. This particularly important as you may find many of the spliced PTMs in your dataset are less well studied and may not have any annotated kinases.">

.. only:: html

  .. image:: /gallery_output/images/thumb/sphx_glr_plot_kstar_enrichment_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_output_plot_kstar_enrichment.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Identify kinases with enriched substrates in differentially included exons, using an adapted version of KSTAR</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="In order to understand how PTMs may be altered due to splicing events, it is useful to identify the flanking sequences of the PTMs and how they may be altered due to nearby splice events (as identified by flanking sequence module). Once we have, this information we can analyze and visualize where the alterations in the flanking sequences occur. First, we need to compare the flanking sequences of PTMs based on whether an exonic region is included or excluded using the compare_flanking_sequences function in PTM-POSE.">

.. only:: html

  .. image:: /gallery_output/images/thumb/sphx_glr_plot_location_altered_flanks_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_output_plot_location_altered_flanks.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Probing where and how PTM flanking sequences are altered</div>
    </div>


.. thumbnail-parent-div-close

.. raw:: html

    </div>


.. toctree::
   :hidden:

   /gallery_output/plot_protein_interactions
   /gallery_output/plot_num_annotations
   /gallery_output/plot_kstar_enrichment
   /gallery_output/plot_location_altered_flanks


.. only:: html

  .. container:: sphx-glr-footer sphx-glr-footer-gallery

    .. container:: sphx-glr-download sphx-glr-download-python

      :download:`Download all examples in Python source code: gallery_output_python.zip </gallery_output/gallery_output_python.zip>`

    .. container:: sphx-glr-download sphx-glr-download-jupyter

      :download:`Download all examples in Jupyter notebooks: gallery_output_jupyter.zip </gallery_output/gallery_output_jupyter.zip>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
