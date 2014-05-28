package EnsembleLearning;
import Jama.*;
// The Gaussian kernel method.

// Input: (Y, Sd, St).
// m: # of drugs. n: # of target proteins.
// Y(m * n) is the DTI adjacency matrix.
// Sd(m * m) is the drug similarity matrix.
// St(n * n) is the target protein similarity matrix.

// Output: F(m * n) as the score matrix of DTI predictions.

class GKOriginal {
	
	Matrix F;
	double sigma = 1;
	int m;
	int n;
	
	GKOriginal (Matrix Y0, Matrix weight, Matrix Sd0, Matrix St0) {
		Matrix Y = CrossValidation.pairwiseProduct(Y0, weight);
		Matrix Sd = Sd0.copy(), St = St0.copy();
		
		m = Y.getRowDimension();
		n = Y.getColumnDimension();
		
		matrixModify(Y, Sd, St);
		
		
//		System.out.println(m + " " + n);
		
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
		Matrix Kd = K_GIP_d.times(1 - alpha).plus(Sd.times(alpha));
		Matrix Kt = K_GIP_t.times(1 - alpha).plus(St.times(alpha));
		Matrix I1 = new Matrix(m, m, 0);
		Matrix I2 = new Matrix(n, n, 0);
		for (int i = 0; i < m; i++)
			I1.set(i, i, sigma);
		for (int i = 0; i < n; i++)
			I2.set(i, i, sigma);
		
		F = (Kd.times(Kd.plus(I1).inverse()).times(Y)).plus(Kt.times(Kt.plus(I2).inverse()).times(Y.transpose()).transpose()).times(0.5);
	}
	
	private void matrixModify(Matrix Y, Matrix Sd, Matrix St) {
		
		double beta = 11;
		double G_corr_drug[][] = new double[m][m];
		double G_corr_target[][] = new double[n][n];
		double com_drug[] = new double[m];
		double com_target[] = new double[n];
		double comdrugmax = 0;
		double comtarmax = 0;
		// compute the difference between two drugs
		// apply a gaussian on it
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
		// compute the difference between two targets
		// apply a gaussian on it
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
		}
		
		//TODO why these lines influence the computation of AUC and AUPR????
		//This is because we run 10-fold crossvalidation, the first round would acts on the influence of later ones
//		/*
		// not cross the bound????
//		System.out.println("m:" +  m + " " + Sd.getRowDimension() + " n:" + n + " " + St.getColumnDimension());
		for (int i = 0; i < m; i++)
			for(int j = 0; j < m; j++)
				Sd.set(i, j, (Sd.get(i, j) + (com_drug[i] + com_drug[j]) * beta * G_corr_drug[i][j] / 2 / comdrugmax) / 2);
//				Sd.set(i, j, 0.8);
		for(int i = 0; i < n ; i ++ )
			for(int j = 0 ; j < n ; j ++ )
				St.set(i, j, (St.get(i, j) + (com_target[i] + com_target[j]) * beta * G_corr_target[i][j] / 2 / comtarmax) / 2);
//				Sd.set(i, j, 0.8);
//			*/	
	}
	
}

