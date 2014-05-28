#rm r*.pr r*.spr r*.roc
# first convert the corresponding file to unix in case of ending up with ^M
# which is not compatible with auc.jar
#dos2unix test$i.txt > null
#size=${1:-10}
size=${1:-5}
array=('WeightedProfile' 'NBI' 'LapRLS' 'GaussianKernel' 'RLS-Kron' '' '' '' '' '')
#array=('' '' '' '' '' '' '' '' '' '')
for ((pt=1; pt <= $size; pt++))
do
	rm "result.single$pt.txt"
	for ((i=0; i < 10; i++))
	do
		t=$(($pt+1))
		cat test$i.txt | awk -F ' ' -v val=$t '{gsub(/([0-9]*:|\r)/, "", $val);print $val}' >> result.single$pt.txt
	done
	paste -d ' ' "result.single$pt.txt" 'record.txt' > single$pt.txt
	java -jar Lib/auc.jar single$pt.txt list > list$pt.txt
	echo 'Number' $pt ${array[$(($pt-1))]} ':'
	tail -n 2 list$pt.txt
done
./ada_eva.sh
