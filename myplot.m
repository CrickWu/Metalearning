cd /home/crick/tmp/java/src
a = dlmread('single1.txt.pr');
x1 = a(:,1);
x2 = a(:,2);
plot(x1, x2);


b = dlmread('test_int.txt.pr');
y1 = b(:,1);
y2 = b(:,2);
plot(y1,y2)
