if ( $#argv < 3 ) then
    echo "Usage:\n  source run_ttree.csh <name> <input_directory> <number of parallel jobs>"
    exit               
endif

set name=$1
set edmdir="/data/gridout/`whoami`/SusyAnalysis/B2G/edmtree/"$name
set ttreedir="/data/gridout/`whoami`/SusyAnalysis/B2G/ttree/"$name
set input_dir=$2
set Nparallel=$3
# check if input dir contains MINIAOD (should be in subdir name)
set is_miniaod=`echo $input_dir | grep "MINIAOD" | wc -l`

ls $input_dir | grep ".root" | awk '{ print "'$input_dir'/"$1 }' >! files.txt
set n=`cat files.txt | wc -l`
if ( -f make_ttrees_$name.csh ) then
    mv make_ttrees_$name.csh make_ttrees_"$name"_prev.csh
endif
if ( -f make_edmtrees_$name.csh ) then
    mv make_edmtrees_$name.csh make_edmtrees_"$name"_prev.csh
endif
foreach i ( `seq 1 $n` )
    set infile=`sed -n "$i"p files.txt`
    set i=`echo $i | awk '{ printf "%03d\n",$1 }'`
    if ( $is_miniaod == 1 ) then
	set edmfile="$edmdir/edmTree_$i.root"
	echo 'cmsRun ../../../B2GAnaFW/B2GAnaFW/test/b2gedmntuples_cfg.py maxEvts=-1 sample="file:'$infile'" outputLabel="file:'$edmfile'" LHE=False' >>! make_edmtrees_$name.csh
	set ttreefile="$ttreedir/b2gTree_$i.root"
	echo 'cmsRun b2gTrees_cfg.py maxEvts=-1 sample="file:'$edmfile'" outputLabel="file:'$ttreefile'"' >>! make_ttrees_$name.csh
    else
	set ttreefile="$ttreedir/b2gTree_$i.root"
	echo 'cmsRun b2gTrees_cfg.py maxEvts=-1 sample="file:'$infile'" outputLabel="file:'$ttreefile'"' >>! make_ttrees_$name.csh
    endif
end
rm files.txt

if ( $is_miniaod == 1 ) then
    mkdir -p $edmdir
    #source source_parallel.csh make_edmtrees_$name.csh $Nparallel
    #rm make_edmtrees_$name.csh
endif
mkdir -p $ttreedir
#source source_parallel.csh make_ttrees_$name.csh $Nparallel
#rm make_ttrees_$name.csh
