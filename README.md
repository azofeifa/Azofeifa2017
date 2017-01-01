# Azofeifa2017
This is a collection of analysis scripts required to recreate the figures in the (hopefully upcoming) Azofeifa 2017 Paper

I will update this README more thoroughly as we get closer to pub. but these python scripts are really in two categories:

1. Generating tables; these are always named with ST_ prefix
2. Generating svgs; these are always with M_ or S_ prefix; in almost all cases, an svg requires .csv file or table output from some ST_ script. 

In all cases, each python script tells your what is required to sucessfully run the script and what will be outputted. This is always in the top portion of the script. I do not specifiy modules however. Off the top of my head, you will likely need:

1. matplotlib
2. numpy
3. scipy
4. scikitlearn
5. seaborn
6. venn
7. pandas
