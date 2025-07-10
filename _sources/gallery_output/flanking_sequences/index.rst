

.. _sphx_glr_gallery_output_flanking_sequences:

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

   /gallery_output/flanking_sequences/plot_sequence_comparison
   /gallery_output/flanking_sequences/plot_sh2_domain_motifs
   /gallery_output/flanking_sequences/plot_location_altered_flanks
   /gallery_output/flanking_sequences/plot_kinase_library_affinity

