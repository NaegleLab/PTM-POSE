
==================
PTM-POSE Reference
==================

#############
Configuration
#############

.. automodule:: ptm_pose.pose_config
	:members: download_translator

##############
PTM Projection
##############

.. automodule:: ptm_pose.project
	:members: find_ptms_in_region, project_ptms_onto_splice_events

##################
Flanking Sequences
##################

.. automodule:: ptm_pose.flanking_sequences
	:members: get_flanking_changes, get_flanking_changes_from_splice_data

#####################
Tool-Specific Modules
#####################

.. autoclass:: ptm_pose.splicing_tools.base.GenericDataset
	:members:

.. autoclass:: ptm_pose.splicing_tools.MATS.MATS_Dataset
	:members:

.. autoclass:: ptm_pose.splicing_tools.MAJIQ.MAJIQ_Dataset
	:members:

.. autoclass:: ptm_pose.splicing_tools.SpliceSeq.SpliceSeq_Dataset
	:members:
	 
###############
Annotating PTMs
###############

.. automodule:: ptm_pose.annotate
	:members:

###############
Analyze Modules
###############

---------
Summaries
---------

.. automodule:: ptm_pose.analyze.summarize
	:members:


-------------------------
Filtering PTMs and Events
-------------------------

.. automodule:: ptm_pose.analyze.filter
	:members:

-----------
Annotations
-----------

.. automodule:: ptm_pose.analyze.annotations
	:members:

--------------------
Protein Interactions
--------------------

.. automodule:: ptm_pose.analyze.interactions
	:members:

-----------------
Enzyme Regulation
-----------------

.. automodule:: ptm_pose.analyze.enzyme
	:members:


------------------
Flanking Sequences
------------------

.. automodule:: ptm_pose.analyze.flank_analysis
	:members:
	
