#=GENOME_DIFF	1.0
#=CREATED	08:40:29 15 Sep 2026
#=PROGRAM	breseq 0.40.2 
#=COMMAND	breseq -r /hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/testdata/parity/synthetic_del/reference.fa -j 4 -o /hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/results/synthetic_del_2646780/breseq_raw /hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/testdata/parity/synthetic_del/reads.fq
#=REFSEQ	/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/testdata/parity/synthetic_del/reference.fa
#=READSEQ	/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/testdata/parity/synthetic_del/reads.fq
#=CONVERTED-BASES	75000
#=CONVERTED-READS	1500
#=INPUT-BASES	75000
#=INPUT-READS	1500
#=MAPPED-BASES	72850
#=MAPPED-READS	1457
DEL	1	2,3	ref_del	501	50
MC	2	.	ref_del	501	550	0	0	left_inside_cov=0	left_outside_cov=45	right_inside_cov=0	right_outside_cov=44
JC	3	.	ref_del	500	-1	ref_del	551	1	0	alignment_overlap=0	coverage_minus=22	coverage_plus=21	flanking_left=50	flanking_right=50	frequency=1	junction_possible_overlap_registers=49	key=ref_del__500__-1__ref_del__551__1__0____50__50__0__0	max_left=49	max_left_minus=49	max_left_plus=46	max_min_left=23	max_min_left_minus=23	max_min_left_plus=23	max_min_right=25	max_min_right_minus=20	max_min_right_plus=25	max_pos_hash_score=98	max_right=49	max_right_minus=49	max_right_plus=49	neg_log10_pos_hash_p_value=0.3	new_junction_coverage=0.85	new_junction_read_count=43	polymorphism_frequency=1.000e+00	pos_hash_score=36	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_overlap=0	side_1_possible_overlap_registers=49	side_1_read_count=0	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_overlap=0	side_2_possible_overlap_registers=49	side_2_read_count=0	side_2_redundant=0	total_non_overlap_reads=43
UN	4	.	ref_del	1	2
UN	5	.	ref_del	501	550
UN	6	.	ref_del	1497	1500
