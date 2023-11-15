# ln-scaling-covenants
Protocols that use simple covenants to scale Lightning. See scalingcovenants_v1.2.pdf.

Also, the small python3 program tt_analysis.py analyzes scalability
given the requirement that all timeout-tree leaves can be put onchain.
That file includes a description of the program's inputs and outputs.

Sample input files in_tt_analysis00.csv, in_tt_analysis01.csv and
in_tt_analysis02.csv assume an annual cost of capital of 0.1%, 1% and
10%, respectively. The corresponding outputs are given as both
unformatted (.csv) files and formatted (.xlsx).
