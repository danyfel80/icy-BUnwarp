package algorithms.danyfel80.registration.bunwarp.old;

/**
 * @author Daniel Felipe Gonzalez Obando
 */
public class MathTools {

	private static final int MAX_SVD_ITERATIONS = 1000;
	private static final double FLT_EPSILON = (double) Float.intBitsToFloat((int) 0x33FFFFFF);;

	/**
	 * B-spline 03.
	 * 
	 * @param x
	 */
	public static double Bspline03(double x) {
		x = Math.abs(x);
		if (x < 1.0F) {
			return (0.5F * x * x * (x - 2.0F) + (2.0F / 3.0F));
		} else if (x < 2.0F) {
			x -= 2.0F;
			return (x * x * x / (-6.0F));
		} else {
			return (0.0F);
		}
	}

	/**
	 * N choose K.
	 *
	 * @param n
	 * @param k
	 */
	public static double nChooseK(int n, int k) {
		if (k > n)
			return 0;
		if (k == 0)
			return 1;
		if (k == 1)
			return n;
		if (k > n / 2)
			k = n - k;
		double prod = 1;
		for (int i = 1; i <= k; i++)
			prod *= (n - k + i) / i;
		return prod;
	}

	/**
	 * Singular Value Decomposition.
	 *
	 * @param U
	 *          input matrix
	 * @param W
	 *          vector of singular values
	 * @param V
	 *          untransposed orthogonal matrix
	 */
	public static void singularValueDecomposition(double[][] U, double[] W, double[][] V) {
		final int lines = U.length;
		final int columns = U[0].length;
		final double[] rv1 = new double[columns];
		double norm, scale;
		double c, f, g, h, s;
		double x, y, z;
		int l = 0;
		int nm = 0;
		boolean flag;
		g = scale = norm = 0.0F;
		for (int i = 0; (i < columns); i++) {
			l = i + 1;
			rv1[i] = scale * g;
			g = s = scale = 0.0F;
			if (i < lines) {
				for (int k = i; (k < lines); k++) {
					scale += Math.abs(U[k][i]);
				}
				if (scale != 0.0) {
					for (int k = i; (k < lines); k++) {
						U[k][i] /= scale;
						s += U[k][i] * U[k][i];
					}
					f = U[i][i];
					g = (0.0 <= f) ? (-Math.sqrt((double) s)) : (Math.sqrt((double) s));
					h = f * g - s;
					U[i][i] = f - g;
					for (int j = l; (j < columns); j++) {
						s = 0.0F;
						for (int k = i; (k < lines); k++) {
							s += U[k][i] * U[k][j];
						}
						f = s / h;
						for (int k = i; (k < lines); k++) {
							U[k][j] += f * U[k][i];
						}
					}
					for (int k = i; (k < lines); k++) {
						U[k][i] *= scale;
					}
				}
			}
			W[i] = scale * g;
			g = s = scale = 0.0F;
			if ((i < lines) && (i != (columns - 1))) {
				for (int k = l; (k < columns); k++) {
					scale += Math.abs(U[i][k]);
				}
				if (scale != 0.0) {
					for (int k = l; (k < columns); k++) {
						U[i][k] /= scale;
						s += U[i][k] * U[i][k];
					}
					f = U[i][l];
					g = (0.0 <= f) ? (-Math.sqrt(s)) : (Math.sqrt(s));
					h = f * g - s;
					U[i][l] = f - g;
					for (int k = l; (k < columns); k++) {
						rv1[k] = U[i][k] / h;
					}
					for (int j = l; (j < lines); j++) {
						s = 0.0F;
						for (int k = l; (k < columns); k++) {
							s += U[j][k] * U[i][k];
						}
						for (int k = l; (k < columns); k++) {
							U[j][k] += s * rv1[k];
						}
					}
					for (int k = l; (k < columns); k++) {
						U[i][k] *= scale;
					}
				}
			}
			norm = ((Math.abs(W[i]) + Math.abs(rv1[i])) < norm) ? (norm) : (Math.abs(W[i]) + Math.abs(rv1[i]));
		}
		for (int i = columns - 1; (0 <= i); i--) {
			if (i < (columns - 1)) {
				if (g != 0.0) {
					for (int j = l; (j < columns); j++) {
						V[j][i] = U[i][j] / (U[i][l] * g);
					}
					for (int j = l; (j < columns); j++) {
						s = 0.0F;
						for (int k = l; (k < columns); k++) {
							s += U[i][k] * V[k][j];
						}
						for (int k = l; (k < columns); k++) {
							if (s != 0.0) {
								V[k][j] += s * V[k][i];
							}
						}
					}
				}
				for (int j = l; (j < columns); j++) {
					V[i][j] = V[j][i] = 0.0F;
				}
			}
			V[i][i] = 1.0F;
			g = rv1[i];
			l = i;
		}
		for (int i = (lines < columns) ? (lines - 1) : (columns - 1); (0 <= i); i--) {
			l = i + 1;
			g = W[i];
			for (int j = l; (j < columns); j++) {
				U[i][j] = 0.0F;
			}
			if (g != 0.0) {
				g = 1.0F / g;
				for (int j = l; (j < columns); j++) {
					s = 0.0F;
					for (int k = l; (k < lines); k++) {
						s += U[k][i] * U[k][j];
					}
					f = s * g / U[i][i];
					for (int k = i; (k < lines); k++) {
						if (f != 0.0) {
							U[k][j] += f * U[k][i];
						}
					}
				}
				for (int j = i; (j < lines); j++) {
					U[j][i] *= g;
				}
			} else {
				for (int j = i; (j < lines); j++) {
					U[j][i] = 0.0F;
				}
			}
			U[i][i] += 1.0F;
		}
		for (int k = columns - 1; (0 <= k); k--) {
			for (int its = 1; (its <= MAX_SVD_ITERATIONS); its++) {
				flag = true;
				for (l = k; (0 <= l); l--) {
					nm = l - 1;
					if ((Math.abs(rv1[l]) + norm) == norm) {
						flag = false;
						break;
					}
					if ((Math.abs(W[nm]) + norm) == norm) {
						break;
					}
				}
				if (flag) {
					c = 0.0F;
					s = 1.0F;
					for (int i = l; (i <= k); i++) {
						f = s * rv1[i];
						rv1[i] *= c;
						if ((Math.abs(f) + norm) == norm) {
							break;
						}
						g = W[i];
						h = EuclideanNorm(f, g);
						W[i] = h;
						h = 1.0F / h;
						c = g * h;
						s = -f * h;
						for (int j = 0; (j < lines); j++) {
							y = U[j][nm];
							z = U[j][i];
							U[j][nm] = y * c + z * s;
							U[j][i] = z * c - y * s;
						}
					}
				}
				z = W[k];
				if (l == k) {
					if (z < 0.0) {
						W[k] = -z;
						for (int j = 0; (j < columns); j++) {
							V[j][k] = -V[j][k];
						}
					}
					break;
				}
				if (its == MAX_SVD_ITERATIONS) {
					return;
				}
				x = W[l];
				nm = k - 1;
				y = W[nm];
				g = rv1[nm];
				h = rv1[k];
				f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0F * h * y);
				g = EuclideanNorm(f, 1.0F);
				f = ((x - z) * (x + z) + h * ((y / (f + ((0.0 <= f) ? (Math.abs(g)) : (-Math.abs(g))))) - h)) / x;
				c = s = 1.0F;
				for (int j = l; (j <= nm); j++) {
					int i = j + 1;
					g = rv1[i];
					y = W[i];
					h = s * g;
					g = c * g;
					z = EuclideanNorm(f, h);
					rv1[j] = z;
					c = f / z;
					s = h / z;
					f = x * c + g * s;
					g = g * c - x * s;
					h = y * s;
					y *= c;
					for (int jj = 0; (jj < columns); jj++) {
						x = V[jj][j];
						z = V[jj][i];
						V[jj][j] = x * c + z * s;
						V[jj][i] = z * c - x * s;
					}
					z = EuclideanNorm(f, h);
					W[j] = z;
					if (z != 0.0F) {
						z = 1.0F / z;
						c = f * z;
						s = h * z;
					}
					f = c * g + s * y;
					x = c * y - s * g;
					for (int jj = 0; (jj < lines); jj++) {
						y = U[jj][j];
						z = U[jj][i];
						U[jj][j] = y * c + z * s;
						U[jj][i] = z * c - y * s;
					}
				}
				rv1[l] = 0.0F;
				rv1[k] = f;
				W[k] = x;
			}
		}
	}

	/**
	 * Euclidean Norm.
	 *
	 * @param a
	 * @param b
	 */
	private static double EuclideanNorm(double a, double b) {
		final double absa = Math.abs(a);
		final double absb = Math.abs(b);
		if (absb < absa) {
			return (absa * Math.sqrt(1.0 + (absb * absb / (absa * absa))));
		} else {
			return ((absb == 0.0F) ? (0.0F) : (absb * Math.sqrt(1.0 + (absa * absa / (absb * absb)))));
		}
	}

	/**
	 * Invert a matrix by the Singular Value Decomposition method.
	 *
	 * @param Ydim
	 *          input, Y-dimension
	 * @param Xdim
	 *          input, X-dimension
	 * @param B
	 *          input, matrix to invert
	 * @param iB
	 *          output, inverted matrix
	 * @return under-constrained flag
	 */
	public static boolean invertMatrixSVD(int Ydim, int Xdim, double[][] B, double[][] iB) {
		boolean underconstrained = false;

		final double[] W = new double[Xdim];
		final double[][] V = new double[Xdim][Xdim];
		// B=UWV^t (U is stored in B)
		singularValueDecomposition(B, W, V);

		// B^-1=VW^-1U^t

		// Compute W^-1
		int Nzeros = 0;
		for (int k = 0; k < Xdim; k++) {
			if (Math.abs(W[k]) < FLT_EPSILON) {
				W[k] = 0.0F;
				Nzeros++;
			} else
				W[k] = 1.0F / W[k];
		}
		if (Ydim - Nzeros < Xdim)
			underconstrained = true;

		// Compute VW^-1
		for (int i = 0; i < Xdim; i++)
			for (int j = 0; j < Xdim; j++)
				V[i][j] *= W[j];

		// Compute B^-1
		// iB should have been already resized
		for (int i = 0; i < Xdim; i++) {
			for (int j = 0; j < Ydim; j++) {
				iB[i][j] = 0.0F;
				for (int k = 0; k < Xdim; k++)
					iB[i][j] += V[i][k] * B[j][k];
			}
		}
		return underconstrained;
	}

	/**
	 * Gives the least-squares solution to (A * x = b) such that (A^T * A)^-1 *
	 * A^T * b = x is a vector of size (column), where A is a (line x column)
	 * matrix, and where b is a vector of size (line). The result may differ from
	 * that obtained by a singular-value decomposition in the cases where the
	 * least-squares solution is not uniquely defined (SVD returns the solution of
	 * least norm, not QR).
	 *
	 * @param A
	 *          An input matrix A[line][column] of size (line x column)
	 * @param b
	 *          An input vector b[line] of size (line)
	 * @return An output vector x[column] of size (column)
	 */
	public static double[] linearLeastSquares(double[][] A, double[] b) {
		if (A == null || A.length == 0)
			return null;
		final int lines = A.length;
		final int columns = A[0].length;
		final double[][] Q = new double[lines][columns];
		final double[][] R = new double[columns][columns];
		final double[] x = new double[columns];
		double s;
		for (int i = 0; (i < lines); i++) {
			for (int j = 0; (j < columns); j++) {
				Q[i][j] = A[i][j];
			}
		}
		QRdecomposition(Q, R);
		for (int i = 0; (i < columns); i++) {
			s = 0.0F;
			for (int j = 0; (j < lines); j++) {
				s += Q[j][i] * b[j];
			}
			x[i] = s;
		}
		for (int i = columns - 1; (0 <= i); i--) {
			s = R[i][i];
			if ((s * s) == 0.0F) {
				x[i] = 0.0F;
			} else {
				x[i] /= s;
			}
			for (int j = i - 1; (0 <= j); j--) {
				x[j] -= R[j][i] * x[i];
			}
		}
		return (x);
	}

	/**
	 * Decomposes the (line x column) input matrix Q into an orthonormal output
	 * matrix Q of same size (line x column) and an upper-diagonal square matrix R
	 * of size (column x column), such that the matrix product (Q * R) gives the
	 * input matrix, and such that the matrix product (Q^T * Q) gives the
	 * identity.
	 *
	 * @param Q
	 *          An in-place (line x column) matrix Q[line][column], which expects
	 *          as input the matrix to decompose, and which returns as output an
	 *          orthonormal matrix
	 * @param R
	 *          An output (column x column) square matrix R[column][column]
	 */
	private static void QRdecomposition(double[][] Q, double[][] R) {
		final int lines = Q.length;
		final int columns = Q[0].length;
		final double[][] A = new double[lines][columns];
		double s;
		for (int j = 0; (j < columns); j++) {
			for (int i = 0; (i < lines); i++) {
				A[i][j] = Q[i][j];
			}
			for (int k = 0; (k < j); k++) {
				s = 0.0F;
				for (int i = 0; (i < lines); i++) {
					s += A[i][j] * Q[i][k];
				}
				for (int i = 0; (i < lines); i++) {
					Q[i][j] -= s * Q[i][k];
				}
			}
			s = 0.0F;
			for (int i = 0; (i < lines); i++) {
				s += Q[i][j] * Q[i][j];
			}
			if ((s * s) == 0.0F) {
				s = 0.0F;
			} else {
				s = 1.0F / Math.sqrt(s);
			}
			for (int i = 0; (i < lines); i++) {
				Q[i][j] *= s;
			}
		}
		for (int i = 0; (i < columns); i++) {
			for (int j = 0; (j < i); j++) {
				R[i][j] = 0.0F;
			}
			for (int j = i; (j < columns); j++) {
				R[i][j] = 0.0F;
				for (int k = 0; (k < lines); k++) {
					R[i][j] += Q[k][i] * A[k][j];
				}
			}
		}
	}

}
