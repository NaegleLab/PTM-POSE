# Frequently Asked Questions
If you do not see your question here, feel free to ask your question with this [form]()

## Data Preparation

**Q: What type of data can be used with PTM-POSE? **

A: Any splicing quantification tools that return the chromosome, DNA strand, and start and stop regions of a given splice event can use PTM-POSE. We have tested this approach predominantly using short-read RNA sequencing analysis from tools like MATS, MAJIQ, and SpliceSeq tools, but PTM-POSE does not require any special dataset.

**Q: Are specific quantification values required for use with PTM-POSE**

A: No, quantification, such as a delta PSI between two experimental groups is not needed for PTM-POSE. However, it can be useful to include this information and indicate using the dPSI_col parameter when running PTM-POSE, as this provides more flexibilty for downstream analysis.




## Interpreting Results

**Q: Why do some PTMs have multiple entries and/or PSI values in the output?**

A: Depending on the quantification tool, you might find cases where there are multiple events related to the same region. PTM-POSE does not filter these out, as it relies on the splicing quantification tool and user to determine relevant splice events. We generally recommend removing PTM sites for which there are multiple conflicting events (different directional change in PSI), as these may be misleading or a consequence of noise in the data.

