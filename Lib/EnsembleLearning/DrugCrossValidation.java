package EnsembleLearning;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import Jama.Matrix;

public class DrugCrossValidation {
	// F: Array of matrices.
	// The first dimension i represents the i-th method.
	// The second dimension j represents the j-th matrix in the 10-fold cross-validation.
	// Namely in each round of 10-fold cross-validation, method i gives 10 score matrices, labeled from F[i][0] to F[i][9].
	int size = 5; // There are totally 5 methods incorporated in our meta approach.
	int foldNum = 10; // 10-fold cross-validation.
	Matrix F[][] = new Matrix[size][foldNum];
	Matrix Y, drugSim, targetSim;
	Matrix[] integrate = new Matrix[foldNum];
	
	// the matrix that base learner can see
	Matrix Y0;
	int testArray[][] = new int[foldNum][];
	
	// weight for positives and negatives
	double pos, neg;
	double realNum;
	
	// totalSet: We divide the dataset of DTIs randomly into 10 small sets, from totalSet[0] to totalSet[9], in nearly equal size.
	// Hence the first parameter ranges from 0 to 9.
	// The second parameter represents the sequence number of an interaction. So totalSet[i][j] represents the j-th interaction in set i.
	// The third parameter keeps the sequence number of the drug (target) corresponding to that interaction.
	// Namely totalSet[i][j][0] = Sequence number of drug, totalSet[i][j][1] = Sequence number of target.
	int[][][] totalSet;
	int[] sizeSet; // sizeSet[i] represents the number of interactions in totalSet[i].
	
	DrugCrossValidation (Matrix Adj, Matrix Sd, Matrix St) throws IOException {
		Y = Adj;
		drugSim = Sd;
		targetSim = St;
		Initialization();
	}
	void Initialization() {
		int m = Y.getRowDimension();
		int n = Y.getColumnDimension();
		//-------------Create ``foldNum'' cross validation assignment------------------------
		// here is the drug's testSize in the ``foldNum''-fold cross validation
		int testSize = m / foldNum;
		int numberOfplusOneArray = m - testSize * foldNum;
		int[] sizeArray = new int[foldNum];
		int[] pointerArray = new int[foldNum];
		// allocate number of drugs per cross validation 
		for (int i = 0; i < numberOfplusOneArray; i++) {
			testArray[i] = new int[testSize + 1];
			sizeArray[i] = testSize + 1;
		}
		for (int i = numberOfplusOneArray; i < foldNum; i++) {
			testArray[i] = new int[testSize];
			sizeArray[i] = testSize;
		}
		// initialize the pointerArray (current number of drugs in each set in the random allocation procedure)
		for (int i = 0; i < foldNum; i++)
			pointerArray[i] = 0;
		
		for (int i = 0; i < m; i++) {
			int ran;
			do {
				ran = (int)(Math.random() * foldNum);
			} while (ran == foldNum || pointerArray[ran] >= sizeArray[ran]);
			testArray[ran][pointerArray[ran]++] = i;
		}
		//------------Test a ``foldNum'' cross validation---------------
		/*
		for (int i = 0; i < foldNum; i++)
			printArray(testArray[i]);
		*/
		
		for (int currentCross = 0; currentCross < foldNum; currentCross++) { // totalSet[i] for the test.
			Y0 = Y.copy();
			double[][] Y0Array = Y0.getArray();
			// delete the testing lines
			for (int i = 0; i < testArray[currentCross].length; i++)
				for (int j = 0; j < n; j++) {
					int delLine = testArray[currentCross][i];
					Y0Array[delLine][j] = 0.0;
				}
			
			WeightedProfileOnDrug wm = new WeightedProfileOnDrug(Y0, Y0, drugSim, 
					targetSim, testArray[currentCross]);
			NBI nm = new NBI(Y0, Y0, drugSim, targetSim);
			LapRLS lm = new LapRLS(Y0, Y0, drugSim, targetSim);
			GaussianKernel km = new GaussianKernel(Y0, Y0, drugSim, targetSim);
			RLSKron rm = new RLSKron(Y0, Y0, drugSim, targetSim);
			
			F[0][currentCross] = wm.F;
			F[1][currentCross] = nm.F;
			F[2][currentCross] = lm.F;
			F[3][currentCross] = km.F;
			F[4][currentCross] = rm.F;
			
			integrate[currentCross] = new Matrix(m, n, 2.0);
		}
		
	}
	
	// Write the information into files for cross-validation, including data-sets and test-sets.
	void WriteToFiles() throws IOException {		
		int n = Y0.getColumnDimension();
		int m = Y0.getRowDimension();
		BufferedWriter bw2 = new BufferedWriter(new FileWriter("record.txt"));	
		BufferedWriter bwi = new BufferedWriter(new FileWriter("test_int.txt"));
		
		for (int i = 0; i < foldNum; i++) {	

			// Files of data-sets.
			// all the written things are done through proper scaling of adaboost
			BufferedWriter bw = new BufferedWriter(new FileWriter("data"+i+".txt"));
			
			// print out all the testArray[j != i] as the test set
			for (int j = 0; j < foldNum; j++) {
				if (i == j) continue;
				for (int k = 0; k < testArray[j].length; k++) {
					int x = testArray[j][k];
					int y;
					for (y = 0; y < n; y++) {
						if (Y.get(x, y) == 1) bw.write("1");
						else bw.write("0");
						for (int l = 0; l < size; l++) {
							int l2 = l + 1;
							bw.write(" "+l2+":"+F[l][j].get(x, y));
						}
						bw.newLine();
					}
				}
			}
			bw.close();
			
			// Files of test-sets.
			bw = new BufferedWriter(new FileWriter("test"+i+".txt"));
			
			for (int k = 0; k < testArray[i].length; k++) {
				int x = testArray[i][k];
				int y;
				for (y = 0; y < n; y++) {
					if (Y.get(x, y) == 1) bw.write("1");
					else bw.write("0");
					for (int l = 0; l < size; l++) {
						int l2 = l + 1;
						bw.write(" "+l2+":"+F[l][i].get(x, y));
					}
					bw.newLine();
					
					if (Y.get(x, y) == 1) bw2.write("1");
					else bw2.write("0");
					bw2.newLine();
					
					bwi.write(integrate[i].get(x, y) + " ");
					if (Y.get(x, y) == 1) bwi.write("1");
					else bwi.write("0");
					bwi.newLine();
				}
			}
			bw.close();
			
		}			
		bw2.close();
		bwi.close();
	}

	static public void printArray(int[] a) {
		for (int i = 0; i < a.length; i++)
			System.out.print(a[i] + " ");
		System.out.println();
	}

}
