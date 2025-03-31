

.. _sphx_glr_gallery_output_Networks:

Impact on protein interaction and regulatory networks
-----------------------------------------------------




.. raw:: html

    <div class="sphx-glr-thumbnails">

.. thumbnail-parent-div-open

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Given that phosphorlaiton are one of the most commonly impacted modifications, there is potential for kinases targeting these sites to be indirectly impacted by alternative splicing through changes in the availability of their substrates. While we provide functions for performing enrichment of known kinase substrates from databases like PhosphoSitePlus, RegPhos, and PTMsigDB, these resources are limited by the overall number of validated substrates (&lt;5%). For this purpose, we have adapted a previously developed algorithm called KSTAR (Kinase Substrate to Activity Relationships) for use with spliced PTM data, which harnesses kinase-substrate predictions to expand the overall number of phosphorylation sites that can be used as evidence. This particularly important as you may find many of the spliced PTMs in your dataset are less well studied and may not have any annotated kinases.">

.. only:: html

  .. image:: /gallery_output/Networks/images/thumb/sphx_glr_plot_kstar_enrichment_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_output_Networks_plot_kstar_enrichment.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Identify kinases with enriched substrates in differentially included exons, using an adapted version of KSTAR</div>
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


.. thumbnail-parent-div-close

.. raw:: html

    </div>


.. toctree::
   :hidden:

   /gallery_output/Networks/plot_kstar_enrichment
   /gallery_output/Networks/plot_protein_interactions

