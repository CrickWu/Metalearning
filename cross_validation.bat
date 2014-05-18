rem Format of this file should be converted to .bat so as to run this program.

REM java -jar Lib\meta.jar 2

for /l %%i in (0,1,9) do (
	java -classpath Lib\libsvm.jar svm_train data%%i.txt
	java -classpath Lib\libsvm.jar svm_predict test%%i.txt data%%i.txt.model o.txt
	del data%%i.txt.model
	del data%%i.txt
	del test%%i.txt
	del o.txt
	ren result.txt result%%i.txt
)

java -jar Lib\ResultsList.jar

java -jar Lib\auc.jar result.txt list >> list.txt

REM for /l %%i in (0,1,9) do (
   REM del result%%i.txt
REM )

REM del result.txt
del result.txt.pr
del result.txt.roc
del result.txt.spr
del record.txt
