package EnsembleLearning;
import Jama.*;
public class Test {
	public static void main(String[] args) {
//		double[][] array = {{1, 2}, {3, 4}}; 
//		Matrix a = new Matrix(array);
//		Matrix b = a.plus(a.transpose());
//		printMatrix(b);
		
		double[] a = {5, 4, 3, 2, 0, 1};
		double[] b = {0, 1, 0, 0, 0, 1};
		qsort(a, b, 0, 5);
		printArray(a);
		printArray(b);
	}
	public static void printArray(double[] a) {
		for (int i = 0; i < a.length; i++)
			System.out.print(a[i] + " ");
		System.out.println();
	}
	public static void printMatrix(Matrix ma) {
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++)
				System.out.print(ma.get(i, j) + " ");
			System.out.println();
		}
		
	}
	private static void qsort(double[] array, double[] real, int s, int t) {
		int i = s, j = t;
		if (s >= t) return;
		double comp = array[i];
		double value = real[i];
		System.out.print(s + " " + t + ": ");
		while (i < j) {
			while (array[j] > comp && i < j) j--;
			if (i < j) {
				array[i] = array[j];
				real[i] = real[j];
				i++;
			}
			System.out.print("\ti: " + i + " " + j + ": ");
			printArray(array);
			while (array[i] < comp && i < j) i++;
			if (i < j) {
				array[j] = array[i];
				real[j] = real[i];
				j--;
			}
			System.out.print("\tj: " + i + " " + j + ": ");
			printArray(array);
		}
		array[i] = comp;
		real[i] = value;
		printArray(array);
		qsort(array, real, s, i - 1);
		qsort(array, real, i + 1, t);
	}
}