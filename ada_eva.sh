echo 'Adaboost:'
java -jar Lib/auc.jar test_int.txt list > list.ada.txt 
tail -n 2 list.ada.txt
#./single_eva.sh
