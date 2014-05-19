package EnsembleLearning;

import Jama.Matrix;

// The network-based inference (NBI) method.

// Input: (Y, Sd, St). ---> Here Sd and St are not used.
// m: # of drugs. n: # of target proteins.
// Y(m * n) is the DTI adjacency matrix.
// Sd(m * m) is the drug similarity matrix.
// St(n * n) is the target protein similarity matrix.

//Output: F(m * n) as the score matrix of DTI predictions.

class NBI { 
	
	Matrix F;
	
	NBI (Matrix Y0, Matrix weight, Matrix Sd, Matrix St) {
		Matrix Y = CrossValidation.pairwiseProduct(Y0, weight);
		
		int m = Y.getRowDimension();
		int n = Y.getColumnDimension();
		double[] kd = new double[m];
		double[] kt = new double[n];
		for (int i = 0; i < m; i++) 
			kd[i] = 0;
		for (int i = 0; i < n; i++) 
			kt[i] = 0;
		for (int i = 0; i < m; i++) 
			for (int j = 0; j < n; j++)
				if (Y.get(i, j) == 1) {
					kd[i]++;
					kt[j]++;
				}
		F = new Matrix(m, n, 0);
		double[][] sum = new double[n][n];
		for (int j = 0; j < n; j++) 
			for (int l = 0; l < n; l++) {
				if (kt[l] == 0) continue;
				double x = 0;
				for (int o = 0; o < m; o++) 
					if (kd[o] != 0)	x += Y.get(o, l) * Y.get(o, j) / kd[o];
				x /= kt[l];
				sum[j][l] = x;
			}
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++) {
				double x = 0;
				for (int l = 0; l < n; l++)
					x += Y.get(i, l) * sum[j][l];
				F.set(i, j, x);
			}
	}
	
}