$muscle_location = "./muscle.exe"
$iqtree = "../iqtree/bin/iqtree3.exe"
$inputfile = "testlist.txt"
$prealignfolder = "aggregates"
$alignfolder = "aligns"
$treefolder = "trees/test"
$genes = @('COX1','COX2','COX3','ATP6','ATP8','CYTB','ND1','ND2','ND3','ND4','ND4L','ND5','ND6')


mkdir $alignfolder -Force | out-null


py -m downloader $inputfile $prealignfolder $genes
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
foreach ($gene in $genes){
    & $muscle_location -align "$prealignfolder/$gene.txt" -output "$alignfolder/$gene.txt" -threads 2
}
$outgroup = get-childitem -path $prealignfolder -file | select-object -first 1 | get-content -first 1 | ForEach-Object { $_.Substring(1) }
& $iqtree -p $alignfolder -m MFP+MERGE --prefix $treefolder -redo -B 1000 -o $outgroup -bnni -alrt 1000
