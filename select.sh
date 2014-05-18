# Usage: ./select.sh file ..
# where .. stands for the columns that is used for an ensemble learning
#file=${*##" *"};
#echo $file
# note here we can not replace $* with $@
echo "1.2 2.2 3.2 4.2" | awk -v val="$*" 'BEGIN{
split(val, arr, " ");
j = 0;
for (i in arr) j++;
}
{
	for (i = 1; i < j + 1; i++)
		printf("%g ", $i);
	print "";
}
END{}'

