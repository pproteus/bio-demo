$muscle_location = "./muscle.exe"
$iqtree = "../iqtree/bin/iqtree3.exe"
$buildfile = "testout.txt"
$alignfile = "testalign.txt"
$treefile = "trees/testtree.txt"

& $muscle_location -align $buildfile -output $alignfile -threads 2
& $iqtree -s $alignfile -pre $treefile -redo -B 1000
