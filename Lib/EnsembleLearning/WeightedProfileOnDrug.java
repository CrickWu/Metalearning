package EnsembleLearning;

import Jama.Matrix;

// The weighted profile (Naive) method.

// Input: (Y, Sd, Sp).
// m: # of drugs. n: # of target proteins.
// Y(m * n) is the DTI adjacency matrix.
// Sd(m * m) is the drug similarity matrix.
// Sp(n * n) is the target protein similarity matrix.

//Output: F(m * n) as the score matrix of DTI predictions.

class WeightedProfileOnDrug {
	
	Matrix F;
	
	WeightedProfileOnDrug (Matrix Y0, Matrix weight, Matrix Sd, Matrix Sp, int[] testLines) {
		Matrix Y = CrossValidation.pairwiseProduct(Y0, weight);
		
		int m = Y.getRowDimension();
		int n = Y.getColumnDimension();
		
		F = new Matrix(m, n);
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++) {
				double x = 0;
				double y = 0;
				for (int k = 0; k < m; k++) {
					// skip the unknown drugs
					for (int l = 0; l < testLines.length; l++)
						if (k == testLines[l])
							continue;
					x += Sd.get(i, k) * Y.get(k, j);
					y += Sd.get(i, k);
				}
				double z = x / y;
				F.set(i, j, z);
			}
	}
	
}
