num=${2:-5}
exp=$(calc "print 2^$num")
for (( i=1; i<$exp; i++))
do
	sel=""
	tmp=$i
	index=1
	while [ $tmp -gt 0 ]
	do
		mod=$(($tmp%2))
		tmp=$(($tmp/2))
		if [ $mod -eq 1 ]
		then
			sel=$sel" "$index
		fi
		index=$(($index+1))
	done
	#echo $sel
	./cross.prob.sh s $sel
done
