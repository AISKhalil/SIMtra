# SIMtra: A tool for generating realistic-simulated WGS data with artificially-introduced CNVs

**SIMtra** is a Python-based software for generating realistic simulated WGS data where CNAs were artificially introduced to serve as ‘ground truths’ for performance evaluation of CNV detection tools. SIMtra utilizes a novel approach to manipulate the original WGS reads of a cancer genome to randomly introduce FAs embedded within the LCVs maintaining the inherent features and complexities of the RD signal. We used SIMtra for generating CNA profiles for the evaluation of [CNAtra tool] (https://github.com/AISKhalil/CNAtra).

**For a full description of the method and applications, please visit [CNAtra Manuscript](https://rdcu.be/b3Cki).**
