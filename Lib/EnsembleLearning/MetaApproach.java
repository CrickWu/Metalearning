package EnsembleLearning;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import Jama.Matrix;

public class MetaApproach {
	
	Map<String, Integer> drug = new HashMap<String, Integer>(); // Keep the ID of each drug.
	Map<String, Integer> target = new HashMap<String, Integer>(); // Keep the ID of each target protein.
	int drugSize = 0, targetSize = 0, adjSize = 0; // Number of drugs, targets and interactions.
	Matrix drugSim, targetSim, adj; // Store all the information in matrices.
	
	
	// Read the file of drug similarity scores.
	void Read1(String src) throws IOException, FileNotFoundException {
		BufferedReader br = new BufferedReader(new FileReader(src));
		String pattern = "([A-Z][0-9]+)\\s";
		Pattern p = Pattern.compile(pattern);
		String s = br.readLine();
		Matcher m = p.matcher(s);
		while (m.find()) {
			drug.put(m.group(1), drugSize);
			drugSize++;
		}
		drugSim = new Matrix(drugSize, drugSize);
		br.close();
		br = new BufferedReader(new FileReader(src));
		s = br.readLine();
		for (int i = 0; i < drugSize; i++) {
			s = br.readLine();
			int label = 0;
			while ((s.charAt(label) >= 'A' && s.charAt(label) <= 'Z') || (s.charAt(label) >= '0' && s.charAt(label) <= '9')) 
				label++;
			while (s.charAt(label) < '0' || s.charAt(label) > '9')
				label++;
			for (int j = 0;; j++) {
				int x = label;
				while (s.charAt(label) == '.' || (s.charAt(label) >= '0' && s.charAt(label) <= '9')) {
					label++;
					if (label >= s.length()) break;
				}
				drugSim.set(i, j, Double.parseDouble(s.substring(x, label)));
				if (j == drugSize - 1) break;
				while (s.charAt(label) < '0' || s.charAt(label) > '9')
					label++;
			}			
		}
		br.close();
	}
	
	// Read the file of target protein similarity scores.
	void Read2(String src) throws IOException, FileNotFoundException { // target.
		BufferedReader br = new BufferedReader(new FileReader(src));
		String pattern = "hsa([0-9]+)\\s";
		Pattern p = Pattern.compile(pattern);
		String s = br.readLine();
		Matcher m = p.matcher(s);
		while (m.find()) {
			target.put(m.group(1), targetSize);
			targetSize++;
		}
		targetSim = new Matrix(targetSize, targetSize);
		br.close();
		br = new BufferedReader(new FileReader(src));
		s = br.readLine();
		for (int i = 0; i < targetSize; i++) {
			s = br.readLine();
			int label = 0;
			while ((s.charAt(label) >= 'a' && s.charAt(label) <= 'z') || (s.charAt(label) >= '0' && s.charAt(label) <= '9')) 
				label++;
			while (s.charAt(label) < '0' || s.charAt(label) > '9')
				label++;
			for (int j = 0;; j++) {
				int x = label;
				while (s.charAt(label) == '.' || (s.charAt(label) >= '0' && s.charAt(label) <= '9')) {
					label++;
					if (label >= s.length()) break;
				}
				targetSim.set(i, j, Double.parseDouble(s.substring(x, label)));
				if (j == targetSize - 1) break;
				while (s.charAt(label) < '0' || s.charAt(label) > '9')
					label++;
			}			
		}
		br.close();
	}
	
	// Read the file of drug-target interactions.
	void Read3(String src) throws IOException, FileNotFoundException {
		adj = new Matrix(drugSize, targetSize, 0);
		BufferedReader br = new BufferedReader(new FileReader(src));
		String pattern = "hsa:([0-9]+)\\s*([A-Z][0-9]+)";
		Pattern p = Pattern.compile(pattern);
		String s;
		while ((s = br.readLine()) != null) {
			Matcher m = p.matcher(s);
			if (m.find()) {
				adj.set(drug.get(m.group(2)), target.get(m.group(1)), 1);
				adjSize++;
			}
		}
		br.close();
	}
	
	// Call methods Read1, Read2 and Read3 for initialization.
	void Initialization(int Class) throws IOException, FileNotFoundException {
		drug.clear();
		target.clear();
		Read1("Data/drug_similarity_"+Class+".txt");
		Read2("Data/target_similarity_"+Class+".txt");
		Read3("Data/interaction_list_"+Class+".txt");
	}
	
	public static void main(String args[]) {
		int classNum = (Integer.parseInt(args[0]));
		MetaApproach meta = new MetaApproach();
		try {
			meta.Initialization(classNum);
		} catch (IOException e) {
			System.out.println("Error: IOException. Source files not found or other reasons.");
		}		
		try {
			CrossValidation cv = new CrossValidation(meta.adj, meta.drugSim, meta.targetSim);
			cv.WriteToFiles();
		} catch (IOException e) {
			System.out.println("Error: IOException. Writing to files failed.");
		}
	}
	
}

class CrossValidation {
	
	// F: Array of matrices.
	// The first dimension i represents the i-th method.
	// The second dimension j represents the j-th matrix in the 10-fold cross-validation.
	// Namely in each round of 10-fold cross-validation, method i gives 10 score matrices, labeled from F[i][0] to F[i][9].
	int size = 5; // There are totally 5 methods incorporated in our meta approach.
	int foldNum = 10; // 10-fold cross-validation.
	Matrix F[][] = new Matrix[size][foldNum];
	Matrix Y, drugSim, targetSim;
	Matrix integrate;
	
	// the matrix that base learner can see
	Matrix Y0;
	
	double realNum;
	
	// totalSet: We divide the dataset of DTIs randomly into 10 small sets, from totalSet[0] to totalSet[9], in nearly equal size.
	// Hence the first parameter ranges from 0 to 9.
	// The second parameter represents the sequence number of an interaction. So totalSet[i][j] represents the j-th interaction in set i.
	// The third parameter keeps the sequence number of the drug (target) corresponding to that interaction.
	// Namely totalSet[i][j][0] = Sequence number of drug, totalSet[i][j][1] = Sequence number of target.
	int[][][] totalSet;
	int[] sizeSet; // sizeSet[i] represents the number of interactions in totalSet[i].	
	
	CrossValidation (Matrix Adj, Matrix Sd, Matrix St) throws IOException {
		Y = Adj;
		drugSim = Sd;
		targetSim = St;
		Initialization();
	}
	
	// By calculation, we get F, totalSet, sizeSet, before cross-validation.
	void Initialization() {
		int m = Y.getRowDimension();
		int n = Y.getColumnDimension();
		int[][] app1 = new int[m][n];
		int[][] app2 = new int[m][n];
		int upSize = m * n / foldNum + 1;
		int trueSize = 0;
		totalSet = new int[foldNum][m * n / foldNum + foldNum][2];
		sizeSet = new int[foldNum];
		for (int i = 0; i < foldNum; i++)
			sizeSet[i] = 0;
		
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++) {
				app1[i][j] = (int)(Y.get(i, j));
				if (app1[i][j] == 1) trueSize++;
			}
		trueSize = trueSize / foldNum + 1;
		
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++) {
				int ran;
				if (Y.get(i, j) == 1) {
					do {
						ran = (int)(Math.random() * foldNum);
					} while (ran == foldNum || sizeSet[ran] >= trueSize);
					totalSet[ran][sizeSet[ran]][0] = i;
					totalSet[ran][sizeSet[ran]][1] = j;
					sizeSet[ran]++;
				}
			}
		
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++) {
				int ran;
				if (Y.get(i, j) == 0) {
					do {
						ran = (int)(Math.random() * foldNum);
					} while (ran == foldNum || sizeSet[ran] >= upSize);
					totalSet[ran][sizeSet[ran]][0] = i;
					totalSet[ran][sizeSet[ran]][1] = j;
					sizeSet[ran]++;
				}
			}
		
		
		for (int i = 0; i < foldNum; i++) { // totalSet[i] for the test.
			for (int i1 = 0; i1 < m; i1++)
				for (int i2 = 0; i2 < n; i2++)
					app2[i1][i2] = app1[i1][i2];
			for (int j = 0; j < sizeSet[i]; j++) {
				if (app1[totalSet[i][j][0]][totalSet[i][j][1]] == 1) 
					app2[totalSet[i][j][0]][totalSet[i][j][1]] = 0;
			}
			Y0 = new Matrix(m, n);
			for(int p = 0; p < m; p++)
				for (int q = 0; q < n; q++)
					Y0.set(p, q, app2[p][q]);
			
//			the initial array is all 1's
			/*
			double[][] weightArray = new double[m][n];
			for (int p = 0; p < m; p++)
				for (int q = 0; q < n; q++)
					// change of weight
//					weightArray[p][q] = 1.0 / (m * n);
					weightArray[p][q] = 1.0;
			Matrix weight = new Matrix(weightArray);*/
			
			//Adaboost test
			Matrix weight = transferMatrix(Y0);
//			printMatrix(weight);
			
			//calculate realNum
			double[][] YArray = Y0.getArray();
			realNum = 0;
			for (int p = 0; p < m; p++)
				for (int q = 0; q < n; q++)
					realNum += YArray[p][q];
			System.out.print(realNum + " ");
			
//			System.out.println("-----"+realNum+"-----");
			//alpha 
			double[] alpha = new double[size];
//			System.out.print(sum(weight) + " ");
			
			WeightedProfile wm = new WeightedProfile(Y0, weight, drugSim, targetSim);	
			alpha[0] = updateWeight(wm.F, weight);
//			System.out.print(sum(weight) + " ");
			
			NBI nm = new NBI(Y0, weight, drugSim, targetSim);
			alpha[1] = updateWeight(nm.F, weight);
//			System.out.print(sum(weight) + " ");
			
			LapRLS lm = new LapRLS(Y0, weight, drugSim, targetSim);
			alpha[2] = updateWeight(lm.F, weight);
//			System.out.print(sum(weight) + " ");
			
			GaussianKernel km = new GaussianKernel(Y0, weight, drugSim, targetSim);	
			alpha[3] = updateWeight(km.F, weight);
			System.out.print(sum(weight) + " ");
			System.out.println();
			
//			RLSKron rm = new RLSKron(Y0, weight, drugSim, targetSim);
//			alpha[4] = updateWeight(rm.F, weight);
			
			F[0][i] = wm.F;
			F[1][i] = nm.F;
			F[2][i] = lm.F;
			F[3][i] = km.F;
//			F[4][i] = rm.F;
			
			integrate = F[0][i].times(alpha[0]);
			for (int p = 1; p < size - 1; p++)
				integrate = integrate.plus(F[p][i].times(alpha[p]));
			
			
			// non-ada compare
			wm = new WeightedProfile(Y0, Y0, drugSim, targetSim);	
			nm = new NBI(Y0, Y0, drugSim, targetSim);
			lm = new LapRLS(Y0, Y0, drugSim, targetSim);
			km = new GaussianKernel(Y0, Y0, drugSim, targetSim);	
			
			RLSKron rm = new RLSKron(Y0, Y0, drugSim, targetSim);
			
			F[0][i] = wm.F;
			F[1][i] = nm.F;
			F[2][i] = lm.F;
			F[3][i] = km.F;
			F[4][i] = rm.F;
		}
	}
	
	// Write the information into files for cross-validation, including data-sets and test-sets.
	void WriteToFiles() throws IOException {		
		BufferedWriter bw2 = new BufferedWriter(new FileWriter("record.txt"));	
		BufferedWriter bwi = new BufferedWriter(new FileWriter("test_int.txt"));
		
		for (int i = 0; i < foldNum; i++) {	

			// Files of data-sets.
			// all the written things are done through proper scaling of adaboost
			BufferedWriter bw = new BufferedWriter(new FileWriter("data"+i+".txt"));
			
			for (int j = 0; j < foldNum; j++) {
				if (i == j) continue;
				for (int k = 0; k < sizeSet[j]; k++) {
					int x = totalSet[j][k][0], y = totalSet[j][k][1];
					if (Y.get(x, y) == 1) bw.write("1");
					else bw.write("0");
					for (int l = 0; l < size; l++) {
						int l2 = l + 1;
						bw.write(" "+l2+":"+F[l][j].get(x, y));
					}
					bw.newLine();
				}
			}
			bw.close();
			
			// Files of test-sets.
			bw = new BufferedWriter(new FileWriter("test"+i+".txt"));
			
			for (int j = 0; j < sizeSet[i]; j++) {
				int x = totalSet[i][j][0], y = totalSet[i][j][1];
				if (Y.get(x, y) == 1) bw2.write("1");
				else bw2.write("0");
				bw2.newLine();
//				bw.write("0");
				if (Y.get(x, y) == 1) bw.write("1");
				else bw.write("0");
				
				for (int l = 0; l < size; l++) {
					int l2 = l + 1;
					bw.write(" "+l2+":"+F[l][i].get(x, y));
				}
				bw.newLine();
				
				// abs?
				bwi.write(Math.abs(integrate.get(x, y)) + " ");
				if (Y.get(x, y) == 1) bwi.write("1");
				else bwi.write("0");
				bwi.newLine();
			}
			bw.close();
	
		}			
		bw2.close();
		bwi.close();
	}
	public static Matrix pairwiseProduct(Matrix a, Matrix b) {
		double[][] aArray = a.getArray();
		double[][] bArray = b.getArray();
		int n = a.getColumnDimension();
		int m = a.getRowDimension();
		
		double[][] retArray = new double[m][n];
		
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				retArray[i][j] = aArray[i][j] * bArray[i][j];
		return new Matrix(retArray);
	}
	
	public Matrix transferMatrix(Matrix ma) {
		int m = ma.getRowDimension();
		int n = ma.getColumnDimension();
		double[][] tArray = new double[m][n];
		double[][] maArray = ma.getArray();
		
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				if (maArray[i][j] == 1)
					
					tArray[i][j] = 1;
				else
					tArray[i][j] = -1.0 / (m * n);
		Matrix tMatrix = new Matrix(tArray);
		return tMatrix;
	}
	
	public static double[] serilize(Matrix Y) {
		double[][] YArray =Y.getArray();
		int m = Y.getRowDimension();
		int n = Y.getColumnDimension();
		
		double[] ret = new double[m * n];
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				ret[i * n + j] = YArray[i][j];
		return ret;
	}
	
	public static double[] serilize(double[][] YArray) {
		int m = YArray.length;
		int n = YArray[0].length;
		double[] ret = new double[m * n];
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				ret[i * n + j] = YArray[i][j];
		return ret;
	}
	
	public double sum(double[] array) {
		double sum = 0;
		for (int i = 0; i < array.length; i++)
			sum += array[i];
		return sum;
	}
	
	public double sum(Matrix ma) {
		double[][] maArray = ma.getArray();
		double sum = 0;
		for (int i = 0; i < ma.getRowDimension(); i++)
			for (int j = 0; j < ma.getColumnDimension(); j++)
				sum += maArray[i][j];
		return sum;
	}
	
	/**
	 * @return double[0] alpha, double[1] cutoff*/
	public double[] classifyError(double[][] pred) {
		double[] predSerial = serilize(pred);
		int length = predSerial.length;
//		System.out.println("length" + length);
//		double[] errArray = new double[length];
		double[] realSerial = serilize(Y0);
		
//		System.out.print("real:" + sum(realSerial) + "pred:" + sum(predSerial));
		
		qsort(predSerial, realSerial, 0, length - 1);
		// since the function only searches for the minError's cutoff
		// we could assume each prediction value is different and we could get the optimal value
		
		// number of positive predictions
		//TODO only have to do once
		double count = 0;
		
//		for (int i = 0; i < length; i++)
//			count += realSerial[i];
		count = realNum;
		
		//TODO the meaning of err is ``right value''
		// that's why we need to finally let the prediction be the absolute value
		// since we actually decrease the weight of the wrong prediction
		
		//TODO but NaN should not appear in this situation
		
		double err = count;
		// for value ">=" cutoff, predict true
		// in normal case, predict each instance to be false does not help improve "Accuracy"
		// thus, we do not check this solution
		int cutoffIndex = 0;
//		System.out.print("real:" + sum(realSerial) + "pred:" + sum(predSerial));
		
		for (int i = 1; i < length; i++) {
//			if (Math.abs(realSerial[i - 1] - 1.0) < 1e-4)
			if ((realSerial[i - 1] - 1.0) < 1e-4)
//			if (realSerial[i - 1]  == 1)
				count = count - 1;
			else
				count = count + 1;
			if (count > err) {
				err = count;
				cutoffIndex = i;
			}
		}
		System.out.println(" err:" + err + " index:" + cutoffIndex);
		err = err / length;
		double alpha = 0.5 * Math.log((1 - err) / err);
		
		double[] retArray = {alpha, predSerial[cutoffIndex]};
		return retArray;
	}
	
	// why real does not hold the 
	private void qsort(double[] array, double[] real, int s, int t) {
		int i = s, j = t;
		if (s >= t) return;
		double comp = array[i];
		double value = real[i];
		while (i < j) {
			while (array[j] > comp && i < j) j--;
			if (i < j) {
				array[i] = array[j];
				real[i] = real[j];
				i++;
			}
			while (array[i] < comp && i < j) i++;
			if (i < j) {
				array[j] = array[i];
				real[j] = real[i];
				j--;
			}
		}
		array[i] = comp;
		real[i] = value;
		qsort(array, real, s, i - 1);
		qsort(array, real, i + 1, t);
	}
	
	
	public double updateWeight(Matrix ma, Matrix weight) {
		double[][] YArray = Y0.getArray();
		double[][] maArray = ma.getArray();
		double[][] weightArray = weight.getArray();
		int m = ma.getRowDimension();
		int n = ma.getColumnDimension();
		double[] errArray = classifyError(maArray);
		
		double alpha = errArray[0];
		// note for value ">=" cutoff, predict true
		double cutoff = errArray[1];
		
		double timesFactor = Math.exp(alpha);
//		System.out.println("times:" + alpha);
		
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++) {
				// if predict true
				if ((maArray[i][j] >= cutoff && YArray[i][j] == 1)
						|| (maArray[i][j] < cutoff && YArray[i][j] != 1))
					weightArray[i][j] *= timesFactor;
				// if incorrect
				else
					weightArray[i][j] /= timesFactor;
			}
		
		// calculate sum
		double sum = 0.0;
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				sum += weightArray[i][j];
		
//		double totalFactor = realNum + (m * n - realNum) / (m * n);
		double totalFactor = realNum - (m * n - realNum) / (m * n);
		// normalize
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				weightArray[i][j] *= totalFactor / sum;
		
		return alpha;
	}
	
	public static void printMatrix(Matrix ma) {
		/*
		double[][] maArray = ma.getArray();
		for (int i = 0; i < maArray.length; i++) {
			for (int j = 0; j < maArray[0].length; j++)
				System.out.print(maArray[i][j] + " ");
			System.out.println();
		}
		*/
		System.out.println(ma.get(0, 0));
	}
}
