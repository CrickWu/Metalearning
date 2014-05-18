#rm r*.pr r*.spr r*.roc
# first convert the corresponding file to unix in case of ending up with ^M
# which is not compatible with auc.jar
#dos2unix test$i.txt > null
size=5
array=('WeightedProfile' 'NBI' 'LapRLS' 'GaussianKernel' 'RLS-Kron')
for ((pt=1; pt <= $size; pt++))
do
	for ((i=0; i < 10; i++))
	do
		t=$(($pt+1))
		cat test$i.txt | awk -v val=$t '{gsub(/.[:\r]/, "", $val);print $val}' >> result.single.txt
	done
	paste -d ' ' 'result.single.txt' 'record.txt' > single$pt.txt
	java -jar Lib/auc.jar single$pt.txt list > list$pt.txt
	echo 'Number' $pt ${array[$(($pt-1))]} ':'
	tail -n 2 list$pt.txt
	rm "result.single.txt"
done
