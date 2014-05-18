#for ((i=0; i<10; i++))
#do
	#java -classpath Lib/libsvm.jar svm_train data$i.txt
	## why can't I generate result$i.txt files??
	#java -classpath Lib/libsvm.jar svm_predict "test$i.txt" data$i.txt.model o.txt 
	#rm data$i.txt.model
	##rm data$i.txt
	##rm "test$i.txt"
	#rm o.txt

	#mv result.txt result$i.txt
#done
t="1 3 4"
./select.py data0.txt $t
