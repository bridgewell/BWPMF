all :

2015-01-01.txt :
	Rscript read.R

2015-01-01.encode.Rds : 2015-01-01.txt
	Rscript encode.R


