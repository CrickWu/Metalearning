package EnsembleLearning;

import Jama.Matrix;

// The Gaussian kernel method.

// Input: (Y, Sd, St).
// m: # of drugs. n: # of target proteins.
// Y(m * n) is the DTI adjacency matrix.
// Sd(m * m) is the drug similarity matrix.
// St(n * n) is the target protein similarity matrix.

// Output: F(m * n) as the score matrix of DTI predictions.

class GaussianKernel {
	
	Matrix F;
	
	GaussianKernel (Matrix Y, Matrix Sd, Matrix St) { 
//		System.out.println();
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
		
//		System.out.println(K_GIP_d.get(0, 0) + " " + K_GIP_d.get(0, 1));
		
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
//		System.out.println(K_GIP_t.get(0, 0) + " " + K_GIP_t.get(0, 1));
		double alpha = 0.5;
		Matrix Kd = K_GIP_d.times(1 - alpha).plus(Sd.times(alpha));
		Matrix Kt = K_GIP_t.times(1 - alpha).plus(St.times(alpha));
//		System.out.println(Kd.get(0, 0) + " " + Kd.get(0, 1));
//		System.out.println(Kt.get(0, 0) + " " + Kt.get(0, 1));
		double sigma = 1;
		Matrix I1 = new Matrix(m, m, 0);
		Matrix I2 = new Matrix(n, n, 0);
		for (int i = 0; i < m; i++)
			I1.set(i, i, sigma);
		for (int i = 0; i < n; i++)
			I2.set(i, i, sigma);
		double beta = 11;
		double G_corr_drug[][] = new double[m][m];
		double G_corr_target[][] = new double[n][n];
		double com_drug[] = new double[m];
		double com_target[] = new double[n];
		double comdrugmax = 0;
		double comtarmax = 0;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				double dis = 0;
				for (int k = 0; k < n; k++) 
					if(Y.get(i, k) != Y.get(j, k)) dis++;
				G_corr_drug[i][j] = Math.exp(-dis / 7);
			}
			double com = 0;
			for (int k = 0; k < n; k++) 
				if (Y.get(i, k) == 1) com++;
			com_drug[i] = com;
			if (comdrugmax < com) comdrugmax = com;
		}
		/*
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				double dis = 0;
				for (int k = 0; k < m; k++) 
					if(Y.get(k, i) != Y.get(k, j)) dis++;
				G_corr_target[i][j] = Math.exp(-dis / 5);
			}
			double com = 0;
			for (int k = 0; k < m; k++) 
				if(Y.get(k, i) == 1) com++;
			com_target[i] = com;
			if (comtarmax < com) comtarmax = com;
		}*/
		Matrix tmp = Kd.times(Kd.plus(I1).inverse().times(Y));
//		System.out.println(tmp.get(0, 0) + " " + tmp.get(0, 1));
		Matrix tmp2 = Kt.times(Kt.plus(I2).inverse()).times(Y.transpose());
//		System.out.println(tmp2.get(0, 0) + " " + tmp2.get(0, 1));
		Matrix tmp3 = Kt.times(Kt.plus(I2).inverse());
//		System.out.println(tmp3.get(0, 0) + " " + tmp3.get(0, 1));
		F = tmp.plus(tmp2.transpose()).times(0.5);
//		System.out.println(F.get(0, 0) + " " + F.get(0, 1));
	}
	
}