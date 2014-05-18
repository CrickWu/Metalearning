package EnsembleLearning;
import Jama.*;
public class Test {
	public static void main(String[] args) {
		double[][] array = {{1, 2}, {3, 4}}; 
		Matrix a = new Matrix(array);
		Matrix b = a.plus(a.transpose());
		printMatrix(b);
	}
	public static void printMatrix(Matrix ma) {
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++)
				System.out.print(ma.get(i, j) + " ");
			System.out.println();
		}
	}
}