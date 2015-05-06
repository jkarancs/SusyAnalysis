set out="plots.root"
set log="plots.log"
set found=0
set args=""
foreach n ( `seq 1 $#argv` )
    if ( $found == 1) then
	set out=$argv[$n]
	set log=`echo $argv[$n] | sed "s;.root;.log;"`
    endif
    set found=`echo $argv[$n] | grep "\-o" | wc -l`
    set args="$args "$argv[$n]
end

make clean
make B2GPlotter
./B2GPlotter $args | tee $log
root 'show_result.C("'$out'")'
