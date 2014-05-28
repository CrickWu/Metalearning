package EnsembleLearning;

import Jama.Matrix;

// The weighted profile (Naive) method.

// Input: (Y, Sd, Sp).
// m: # of drugs. n: # of target proteins.
// Y(m * n) is the DTI adjacency matrix.
// Sd(m * m) is the drug similarity matrix.
// Sp(n * n) is the target protein similarity matrix.

//Output: F(m * n) as the score matrix of DTI predictions.

class WeightedProfile {
	
	Matrix F;
	
	WeightedProfile (Matrix Y0, Matrix weight, Matrix Sd, Matrix Sp) {
		Matrix Y = CrossValidation.pairwiseProduct(Y0, weight);
		
		int m = Y.getRowDimension();
		int n = Y.getColumnDimension();
		
		Matrix K_GIP_d = new Matrix(m, m);
		Matrix K_GIP_t = new Matrix(n, n);
		double rd = 0, rt = 0;
		for (int i = 0; i < m; i++) 
			for (int j = 0; j < n; j++)
				rd += Y.get(i, j);
		rt = 1.0 * n / rd;
		rd = 1.0 * m / rd;
		for (int i = 0; i < m; i++)
			for (int j = i; j < m; j++) {
				if (i == j) K_GIP_d.set(i, j, 1);
				else {
					double delta = 0;
					for (int k = 0; k < n; k++) 
						delta += (Y.get(i, k) - Y.get(j, k)) * (Y.get(i, k) - Y.get(j, k));
					K_GIP_d.set(i, j, Math.exp(-rd * delta));
					K_GIP_d.set(j, i, Math.exp(-rd * delta));
				}
			}
		for (int i = 0; i < n; i++)
			for (int j = i; j < n; j++) {
				if (i == j) K_GIP_t.set(i, j, 1);
				else {
					double delta = 0;
					for (int k = 0; k < m; k++) 
						delta += (Y.get(k, i) - Y.get(k, j)) * (Y.get(k, i) - Y.get(k, j));
					K_GIP_t.set(i, j, Math.exp(-rt * delta));
					K_GIP_t.set(j, i, Math.exp(-rt * delta));
				}
			}
		double alpha = 0.5;
		Matrix Wd = K_GIP_d.times(1 - alpha).plus(Sd.times(alpha));
		Matrix Wp = K_GIP_t.times(1 - alpha).plus(Sp.times(alpha));
		
		F = new Matrix(m, n);
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++) {
				double x = 0;
				double y = 0;
				for (int k = 0; k < m; k++) {
					x += Sd.get(i, k) * Y.get(k, j);
					y += Sd.get(i, k);
//					x += Wd.get(i, k) * Y.get(k, j);
//					y += Wd.get(i, k);
				}
				double z = x / y;
				x = 0;
				y = 0;
				for (int k = 0; k < n; k++) {
					x += Sp.get(j, k) * Y.get(i, k);
					y += Sp.get(j, k);
//					x += Wp.get(j, k) * Y.get(i, k);
//					y += Wp.get(j, k);
				}
				z += x / y;
				F.set(i, j, z / 2);
			}
	}
	
}
