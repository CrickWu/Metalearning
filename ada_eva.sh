echo 'Adaboost:'
java -jar Lib/auc.jar test_int.txt list > result.ada.txt
tail -n 2 result.ada.txt
#./single_eva.sh
