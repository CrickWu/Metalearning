package EnsembleLearning;

import Jama.Matrix;

// The Laplacian regularized least squares (LapRLS) method.

// Input: (Y, Sd, Sp).
// m: # of drugs. n: # of target proteins.
// Y0(m * n) is the DTI adjacency matrix.
// weight(m * n) is the DTI weight matrix.
// Sd(m * m) is the drug similarity matrix.
// Sp(n * n) is the target protein similarity matrix.

// Output: F(m * n) as the score matrix of DTI predictions.

class LapRLS {
	
	Matrix F;
	
	LapRLS (Matrix Y0, Matrix weight, Matrix Sd, Matrix Sp) {
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

		/*
		Matrix Kd = new Matrix(m, m);
		Matrix Kp = new Matrix(n, n);
		for (int i = 0; i < m; i++)
			for (int j = 0; j < m; j++) {
				double x = 0;
				for (int k = 0; k < n; k++)
					x += Y.get(i, k) * Y.get(j, k);
				Kd.set(i, j, x);
			}
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++) {
				double x = 0;
				for (int k = 0; k < m; k++)
					x += Y.get(k, i) * Y.get(k, j);
				Kp.set(i, j, x);
			}
		double r1 = 1, r2 = 0.01;
		Matrix Wd = Sd.times(r1).plus(Kd.times(r2)).times(1.0 / (r1 + r2));
		Matrix Wp = Sp.times(r1).plus(Kp.times(r2)).times(1.0 / (r1 + r2));
		*/
		Matrix Dd = new Matrix(m, m, 0);
		Matrix Dp = new Matrix(n, n, 0);
		Matrix I1 = new Matrix(m, m, 0);
		Matrix I2 = new Matrix(n, n, 0);
		for (int i = 0; i < m; i++) {
			I1.set(i, i, 1);
			double x = 0;
			for (int j = 0; j < m; j++)
				x += Wd.get(i, j);
			Dd.set(i, i, 1.0 / Math.sqrt(x));
		}
		for (int i = 0; i < n; i++) {
			I2.set(i, i, 1);
			double x = 0;
			for (int j = 0; j < n; j++)
				x += Wp.get(i, j);
			Dp.set(i, i, 1.0 / Math.sqrt(x));
		}
		Double beta = 0.3;
		Double sigma = 0.3;
		Matrix asst1 = Wd.plus((I1.minus(Dd.times(Wd).times(Dd))).times(Wd).times(beta)).plus(I1.times(sigma));
		Matrix asst2 = Wp.plus((I2.minus(Dp.times(Wp).times(Dp))).times(Wp).times(beta)).plus(I2.times(sigma));
		Matrix Fd = Wd.times(asst1.inverse()).times(Y);
		Matrix Fp = Wp.times(asst2.inverse()).times(Y.transpose());
		F = Fd.plus(Fp.transpose()).times(0.5);
	}
}