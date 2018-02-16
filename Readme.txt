This program takes total 5 command line arguments as input
	1st arg : input graph file name (in Galib format)
	2nd arg : fraction of shuffle, range (0.0, 1.0], t = 0.5, 1.0 means 50% and 100% edge switch respectively
	3rd arg : output graph file name
	4th arg : number of steps (it is an integer number >= 1), edges are switched in these number of steps
	5th arg : write the output file?
		0 means: No
		1 means: Yes, in Galib format
		2 means: Yes, in edge list format
