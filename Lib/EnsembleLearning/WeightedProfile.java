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
		F = new Matrix(m, n);
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++) {
				double x = 0;
				double y = 0;
				for (int k = 0; k < m; k++) {
					x += Sd.get(i, k) * Y.get(k, j);
					y += Sd.get(i, k);
				}
				double z = x / y;
				x = 0;
				y = 0;
				for (int k = 0; k < n; k++) {
					x += Sp.get(j, k) * Y.get(i, k);
					y += Sp.get(j, k);
				}
				z += x / y;
				F.set(i, j, z / 2);
			}
	}
	
}
