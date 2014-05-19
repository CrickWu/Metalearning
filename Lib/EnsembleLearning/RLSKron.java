package EnsembleLearning;
import Jama.*;

//The network-based inference (NBI) method.

//Input: (Y, Sd, St). ---> Here Sd and St are not used.
//m: # of drugs. n: # of target proteins.
//Y0(m * n) is the DTI adjacency matrix.
//weight(m * n) is the DTI weight matrix.
//Sd(m * m) is the drug similarity matrix.
//St(n * n) is the target protein similarity matrix.

//Output: F(m * n) as the score matrix of DTI predictions.

//TODO
// Actually, we first implement RLS-avg since it performs similarly
// with RLSKron and requires less running and storing resources
public class RLSKron {
	public Matrix F;
	private int m;
	private int n;
	private double bandwidth = 1.0;
	private double alphaD = 0.5;
	private double alphaT = 0.5;
	private double sigma = 1.0;
	private double gammaD, gammaT;
	
	RLSKron(Matrix Y0, Matrix weight, Matrix Sd, Matrix St) {
		
		Matrix Y = CrossValidation.pairwiseProduct(Y0, weight);
		
		m = Y.getRowDimension();
		n = Y.getColumnDimension();
		
		double[][] YArray = Y.getArray();
		double[][] KGIPdArray = new double[m][m];
		double[][] KGIPtArray = new double[n][n];
		
		computeGamma(YArray);
		// compute K_GIP,d and K_GIP,t (Gaussian Interaction Profile)
		for (int i = 0; i < m; i++)
			for (int j = i; j < m; j++) {
				KGIPdArray[i][j] = KGIPdArray[j][i]
						= GIPKernel(YArray[i], YArray[j], gammaD);
			}		
		
		// different values of gamma
		Matrix Yt = Y.transpose();
		double[][] YtArray = Yt.getArray();
		for (int i = 0; i < n; i++)
			for (int j = i; j < n; j++) {
				KGIPtArray[i][j] = KGIPtArray[j][i]
						= GIPKernel(YtArray[i], YtArray[j], gammaT);
			}
		
		Matrix KGIPd = new Matrix(KGIPdArray);
		Matrix KGIPt = new Matrix(KGIPtArray);
		
		// TODO change back to symKernel
		Matrix Kd = symKernel(Sd).times(alphaD).plus(KGIPd.times(1.0 - alphaD));
		Matrix Kt = symKernel(St).times(alphaT).plus(KGIPt.times(1.0 - alphaT));
		
//		F = computeAvg(Kd, Kt, Y);
		F = computeKron(Kd, Kt, Y);
//		System.out.println(F.get(0, 0) + " " + F.get(0, 1));
	}
	private Matrix computeKron(Matrix Kd, Matrix Kt, Matrix Y) {
		EigenvalueDecomposition eigD = Kd.eig();
		EigenvalueDecomposition eigT = Kt.eig();
		double[] ld = eigD.getRealEigenvalues(),
				lt = eigT.getRealEigenvalues();
		// manually computation of the new l matrix
		double[][] l = new double[m][n];
		
		//l = ld * lt^T
		//l = l ./ (l + sigma)
		for (int i = 0; i < m; i ++)
			for (int j = 0; j < n; j++) {
				l[i][j] = ld[i] * lt[j];
				l[i][j] = l[i][j] / (l[i][j] + sigma);
			}
		
		Matrix Vd = eigD.getV();
		Matrix Vt = eigT.getV();
		
		Matrix m1 = Vd.transpose().times(Y).times(Vt);
		// m2 = m1 .* l (here let m2 = m1)
		double[][] m1Array = m1.getArray();
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				m1Array[i][j] = m1Array[i][j] * l[i][j];
		
		return Vd.times(m1).times(Vt.transpose());
		
	}
	
	private Matrix computeAvg(Matrix Kd, Matrix Kt, Matrix Y) {
		Matrix left = Kd.times(Kd.plus(Matrix.identity(m, m).times(sigma)).inverse()).times(Y).times(0.5);
		Matrix right = Kt.times(Kt.plus(Matrix.identity(n, n).times(sigma)).inverse()).times(Y.transpose()).times(0.5);
		return left.plus(right.transpose());
	}
	// gamma = bandwidthParameter / (1/n_d sum |y_di|^2)
	// n_d <- m, sum |y_di|^2 = sum(YArray), bandwidthParameter = 1.0 in the paper
	private void computeGamma(double[][] YArray) {
		double ret = 0.0;
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				ret += YArray[i][j];
		gammaD = bandwidth * m / ret;
		gammaT = bandwidth * n / ret;
	}
	
	// this function computes the GIPKernel
	// i.e. exp(-gamma_d|y_di-y_dj|^2)
	static private double GIPKernel(double[] di, double[] dj, double gamma) {
		double ret = 0.0;
		for (int i = 0; i < di.length; i++)
			ret += (di[i] - dj[i]) * (di[i] - dj[i]);
		ret = Math.exp(-gamma * ret);
		return ret;
	}
	
	static private Matrix tsymKernel(Matrix ma) {
		return ma;
	}
	
	static private Matrix symKernel(Matrix ma) {
		int dim = ma.getColumnDimension();
		Matrix s = ma.plus(ma.transpose()).times(0.5);
//		System.out.println(s.getRowDimension() + " " + s.getColumnDimension());
		Matrix tmp = s.eig().getD();
//		System.out.println(tmp.getRowDimension() + " " + tmp.getColumnDimension());
		double min = tmp.get(dim - 1, dim - 1);
		// plus a small multiple of identity matrix to enforce positive definite property
		if (min > 0.0)
			return s;
		else {
			double plusFactor = Math.ceil(-min);
			return s.plus(Matrix.identity(dim, dim).times(plusFactor));
		}
	}
	
}
