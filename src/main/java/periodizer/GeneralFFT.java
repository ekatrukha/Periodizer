package periodizer;


/** this class/code (referenced above) was modified and extended here, thank you, Nayuki! 
 * 
 * mostly changed precision  to float (since input images almost never double)
 * **/

public final class GeneralFFT {
	
	// pre-computed tables of cos and sin 
	// for speed performance
	// Z-chirp cos/sin in Bluestein part
	double [] sinXB;
	double [] cosXB;
		
	//helpful stuff
	float[] yreal;
	float[] yimag;
	

	/**
	 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 * The vector can have any length. This is a wrapper function.
	 */
	public static void transform(float[] real, float[] imag) {
		int n = real.length;
		if (n != imag.length)
			throw new IllegalArgumentException("Mismatched lengths");
		if (n == 0)
			return;
		else if ((n & (n - 1)) == 0)  // Is power of 2
			transformRadix2(real, imag);
		else  // More complicated algorithm for arbitrary sizes
			transformBluestein(real, imag);
	}
	
	/** 
	 * Computes the discrete Fourier transform (DFT) of the given real vector, returning complex array [0]=real and [1]=imag.
	 * The vector can have any length. This is a wrapper function.
	 */
	public static float[][] transformR(float[] realin) {
		int n = realin.length;
		
		float[][] fin = new float [2][n];
		System.arraycopy(realin, 0, fin[0], 0, n);
		transform(fin[0],fin[1]);
	
		return fin;		
	}
	
	/** 
	 * Computes the discrete Fourier transform (DFT) of the given real vector, returning complex array [0]=real and [1]=imag.
	 * The vector can have any length. Uses pre-computed tables of cos and sin for defined n.
	 */
	public float[][] prC_transformR(float[] realin) {
		int n = realin.length;
		
		float[] real = new float [n];
		float[] imag = new float [n];
		float[][] fin = new float [2][n];
		System.arraycopy(realin, 0, real, 0, n);
		prC_transform(real,imag);
		fin[0]=real;
		fin[1]=imag;
		return fin;
		
	}
	/** 
	 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 * The vector can have any length. Uses pre-computed tables of cos and sin for defined n.
	 */
	public void prC_transform(float[] real, float[] imag) {
		int n = real.length;
		if (n != imag.length)
			throw new IllegalArgumentException("Mismatched lengths");
		if (n == 0)
			return;
		else if ((n & (n - 1)) == 0)  // Is power of 2
			prC_transformRadix2(real, imag);
		else  // More complicated algorithm for arbitrary sizes
			prC_transformBluestein(real, imag);
	}
	
	/** 
	 * Computes the inverse discrete Fourier transform (IDFT) of the given complex vector, storing the result back into the vector.
	 * The vector can have any length. This is a wrapper function. 
	 **/
	public static void inverseTransform(float[] real, float[] imag) {
		int n = real.length;
		for(int i=0;i<n;i++)
		{
			real[i]/=n;
			imag[i]/=n;
		}
		transform(imag, real);
	}
	
	/** 
	 * Computes the inverse discrete Fourier transform (IDFT) of the given complex vector, storing the result back into the vector.
	 * The vector can have any length. Uses pre-computed tables of cos and sin for defined n.
	 */
	public void prC_inverseTransform(float[] real, float[] imag) 
	{
		int n = real.length;
		for(int i=0;i<n;i++)
		{
			real[i]/=n;
			imag[i]/=n;
		}
		
		prC_transform(imag, real);
	}

	/**
	 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 * The vector's length must be a power of 2. Uses the Cooley-Tukey decimation-in-time radix-2 algorithm.
	 */
	public static void transformRadix2(float[] real, float[] imag) {
		// Length variables
		int n = real.length;
		if (n != imag.length)
			throw new IllegalArgumentException("Mismatched lengths");
		int levels = 31 - Integer.numberOfLeadingZeros(n);  // Equal to floor(log2(n))
		if (1 << levels != n)
			throw new IllegalArgumentException("Length is not a power of 2");
		
		// Trigonometric tables
		double[] cosTable = new double[n / 2];
		double[] sinTable = new double[n / 2];
		double valI;
		for (int i = 0; i < n / 2; i++) {
			valI = 2.0 * i* Math.PI/ n;
			cosTable[i] = Math.cos(valI);
			sinTable[i] = Math.sin(valI);
		}
		
		// Bit-reversed addressing permutation
		for (int i = 0; i < n; i++) {
			int j = Integer.reverse(i) >>> (32 - levels);
			if (j > i) {
				float temp = real[i];
				real[i] = real[j];
				real[j] = temp;
				temp = imag[i];
				imag[i] = imag[j];
				imag[j] = temp;
			}
		}
		
		// Cooley-Tukey decimation-in-time radix-2 FFT
		for (int size = 2; size <= n; size *= 2) {
			int halfsize = size / 2;
			int tablestep = n / size;
			for (int i = 0; i < n; i += size) {
				for (int j = i, k = 0; j < i + halfsize; j++, k += tablestep) {
					int l = j + halfsize;
					double tpre =  real[l] * cosTable[k] + imag[l] * sinTable[k];
					double tpim = -real[l] * sinTable[k] + imag[l] * cosTable[k];
					real[l] = real[j] - (float)tpre;
					imag[l] = imag[j] - (float)tpim;
					real[j] += tpre;
					imag[j] += tpim;
				}
			}
			if (size == n)  // Prevent overflow in 'size *= 2'
				break;
		}
	}
	
	
	public void NR_fourier(float[] real, float[] imag)
	{
		int n = real.length;
		int i;
		float [] datain = new float [2*n];
		for (i = 0; i < n; i++)
		{
			datain[i*2]=real[i];
			datain[i*2+1]=imag[i];
		}
		four1(datain, 1);
		for (i=0;i<n;i++)
		{
			real[i]=datain[i*2];
			imag[i]=datain[i*2+1];
		}
	}
	/** 
	 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 * The vector's length must be a power of 2. Uses the Cooley-Tukey decimation-in-time radix-2 algorithm.
	 */
	public void prC_transformRadix2(float[] real, float[] imag) 
	{
		
		NR_fourier(real, imag);
	}
	
	/** 
	 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 * The vector can have any length. This requires the convolution function, which in turn requires the radix-2 FFT function.
	 * Uses Bluestein's chirp z-transform algorithm.
	 */
	public static void transformBluestein(float[] real, float[] imag) {
		// Find a power-of-2 convolution length m such that m >= n * 2 + 1
		int n = real.length;
		if (n != imag.length)
			throw new IllegalArgumentException("Mismatched lengths");
		if (n >= 0x20000000)
			throw new IllegalArgumentException("Array too large");
		int m = Integer.highestOneBit(n) * 4;
		
		// Trigonometric tables
		double[] cosTable = new double[n];
		double[] sinTable = new double[n];
		double valI;
		for (int i = 0; i < n; i++) {
			//int j = (int)((long)i * i % (n * 2));  // This is more accurate than j = i * i
			valI = ((int)((long)i * i % (n * 2)))*Math.PI / n;
			cosTable[i] = Math.cos(valI);
			sinTable[i] = Math.sin(valI);
		}
		
		// Temporary vectors and preprocessing
		float[] areal = new float[m];
		float[] aimag = new float[m];
		for (int i = 0; i < n; i++) {
			areal[i] = (float)(real[i] * cosTable[i] + imag[i] * sinTable[i]);
			aimag[i] = (float)(-real[i] * sinTable[i] + imag[i] * cosTable[i]);
		}
		float[] breal = new float[m];
		float[] bimag = new float[m];
		breal[0] = (float)cosTable[0];
		bimag[0] = (float)sinTable[0];
		for (int i = 1; i < n; i++) {
			breal[i] = breal[m - i] = (float)cosTable[i];
			bimag[i] = bimag[m - i] = (float)sinTable[i];
		}
		
		// Convolution
		float[] creal = new float[m];
		float[] cimag = new float[m];
		convolve(areal, aimag, breal, bimag, creal, cimag);
		
		// Postprocessing
		for (int i = 0; i < n; i++) {
			real[i] = (float) (creal[i] * cosTable[i] + cimag[i] * sinTable[i]);
			imag[i] = (float)(-creal[i] * sinTable[i] + cimag[i] * cosTable[i]);
		}
	}
	
	/** 
	 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 * The vector can have any length. This requires the convolution function, which in turn requires the radix-2 FFT function.
	 * Uses Bluestein's chirp z-transform algorithm.
	 */
	public void prC_transformBluestein(float[] real, float[] imag) {
		// Find a power-of-2 convolution length m such that m >= n * 2 + 1
		int n = real.length;
		if (n != imag.length)
			throw new IllegalArgumentException("Mismatched lengths");
		if (n >= 0x20000000)
			throw new IllegalArgumentException("Array too large");
		int m = Integer.highestOneBit(n) * 4;
		
		
		// Temporary vectors and preprocessing
		float[] areal = new float[m];
		float[] aimag = new float[m];
		for (int i = 0; i < n; i++) {
			areal[i] = (float)(real[i] * cosXB[i] + imag[i] * sinXB[i]);
			aimag[i] = (float)(-real[i] * sinXB[i] + imag[i] * cosXB[i]);
		}

		// Convolution
		float[] creal = new float[m];
		float[] cimag = new float[m];
		prC_convolve(areal, aimag, creal, cimag);
		
		// Postprocessing
		for (int i = 0; i < n; i++) {
			real[i] = (float) (creal[i] * cosXB[i] + cimag[i] * sinXB[i]);
			imag[i] = (float)(-creal[i] * sinXB[i] + cimag[i] * cosXB[i]);
		}
	}

	
	/** 
	 * Computes the circular convolution of the given real vectors. Each vector's length must be the same.
	 */
	public static void convolve(float[] xvec, float[] yvec, float[] outvec) {
		int n = xvec.length;
		if (n != yvec.length || n != outvec.length)
			throw new IllegalArgumentException("Mismatched lengths");
		convolve(xvec, new float[n], yvec, new float[n], outvec, new float[n]);
	}
	
	
	/** 
	 * Computes the circular convolution of the given complex vectors. Each vector's length must be the same.
	 */
	public static void convolve(float[] xreal, float[] ximag,
			float[] yreal, float[] yimag, float[] outreal, float[] outimag) {
		
		int n = xreal.length;
		if (n != ximag.length || n != yreal.length || n != yimag.length
				|| n != outreal.length || n != outimag.length)
			throw new IllegalArgumentException("Mismatched lengths");
		
		xreal = xreal.clone();
		ximag = ximag.clone();
		yreal = yreal.clone();
		yimag = yimag.clone();
		transform(xreal, ximag);
		transform(yreal, yimag);
		
		for (int i = 0; i < n; i++) {
			double temp = xreal[i] * yreal[i] - ximag[i] * yimag[i];
			ximag[i] = ximag[i] * yreal[i] + xreal[i] * yimag[i];
			xreal[i] = (float)temp;
		}
		//inverse non-scaled transform
		transform(ximag,xreal);
		
		for (int i = 0; i < n; i++) {  // Scaling (because this FFT implementation omits it)
			outreal[i] = xreal[i] / n;
			outimag[i] = ximag[i] / n;
		}
	}
	
	/**
	 * Computes the circular convolution of the given complex vectors. 
	 * Each vector's length must be the same.
	 */
	public void prC_convolve(float[] xreal, float[] ximag, float[] outreal, float[] outimag) 
	{
		
		int n = xreal.length;
		if (n != ximag.length || n != yreal.length || n != yimag.length
				|| n != outreal.length || n != outimag.length)
			throw new IllegalArgumentException("Mismatched lengths");
		
		xreal = xreal.clone();
		ximag = ximag.clone();
		//yreal = yreal.clone();
		//yimag = yimag.clone();
		prC_transform(xreal, ximag);
		//prC_transform(yreal, yimag);
		
		for (int i = 0; i < n; i++) {
			double temp = xreal[i] * yreal[i] - ximag[i] * yimag[i];
			ximag[i] = ximag[i] * yreal[i] + xreal[i] * yimag[i];
			xreal[i] = (float)temp;
		}
		//inverse non scaled transform
		prC_transform(ximag, xreal);
		//prC_inverseTransform(xreal, ximag);
		
		for (int i = 0; i < n; i++) {  // Scaling (because this FFT implementation omits it)
			outreal[i] = xreal[i] / n;
			outimag[i] = ximag[i] / n;
		}
	}
	
 	/** Function pre-computes sin and cos values
 	 * if FFT will happen on the same length arrays,
 	 * to speed up calculations **/
 	public void preComputeCosSin(int N)
 	{
 		int i;
		// Trigonometric tables
		cosXB = new double[N];
		sinXB = new double[N];
		double valI;
		for (i = 0; i < N; i++) 
		{
			valI= ((int)((long)i * i % (N * 2)))*Math.PI / N;
			cosXB[i] = Math.cos(valI);
			sinXB[i] = Math.sin(valI);
		}
		int m = Integer.highestOneBit(N) * 4;
		
		// Temporary vectors and preprocessing

		yreal = new float[m];
		yimag = new float[m];
		yreal[0] = (float)cosXB[0];
		yimag[0] = (float)sinXB[0];
		for (i = 1; i < N; i++) {
			yreal[i] = yreal[m - i] = (float)cosXB[i];
			yimag[i] = yimag[m - i] = (float)sinXB[i];
		}
		transform(yreal,yimag);
 	}
 	
	/** Calculates DFT transform of complex input 
	 *  check "Numerical Recipes for explanation  
	 *  @param data supposed to be 0,2,4,... real part, 1,3,5... imaginary part
	 *  @param isign specifies whether it is forward or inverse transform
	 * **/
	public static void four1(float[] data, int isign)
	{
		int nn, mmax, m, j, istep, i;
		float wtemp, wr, wpr, wpi, wi, theta, tempr, tempi, tempf;
		
		nn = data.length;
		
		int n=nn>>1;
		j=1;
		for(i=1;i<nn;i+=2)
		{
			if(j>i)
			{
				tempf=data[j-1];
				data[j-1]=data[i-1];
				data[i-1]=tempf;
				tempf=data[j];
				data[j]=data[i];
				data[i]=tempf;
			}
			m=n;
			while (m>=2 && j>m)
			{
				j-=m;
				m>>=1;
			}
			j+=m;
		}
		//Danielson-Lanczos routine section
		mmax=2;
		while(nn>mmax)
		{
			istep=mmax<<1;
			theta = (float) (isign*(2.0*Math.PI/mmax));
			wtemp = (float) Math.sin(0.5*theta);
			wpr=-2.0f*wtemp*wtemp;
			wpi=(float)Math.sin(theta);
			wr=1.0f;
			wi=0.0f;
			for(m=1;m<mmax;m+=2)
			{
				for(i=m;i<=nn;i+=istep)
				{
					j=i+mmax;
					tempr=wr*data[j-1]-wi*data[j];
					tempi=wr*data[j]+wi*data[j-1];
					data[j-1]=data[i-1] - tempr;
					data[j]=data[i]-tempi;
					data[i-1]+=tempr;
					data[i]+=tempi;
				}
				wtemp=wr;
				wr=wtemp*wpr-wi*wpi+wr;
				wi=wi*wpr+wtemp*wpi+wi;
			}
			mmax=istep;
			
		}
	}
 	
 	/*
	public static void main( String... args) throws Exception
	{
		
		
		float [] testArr = new float[]{20f,21f,22f,23f,24f,23f,22f,21f,45f};
		int size =testArr.length;
		float[][] FFTtr = GeneralFFT.transformR(testArr);
		
		
		//for(int i=0;i<FFTtr[0].length;i++)
		//{
		//	FFTtr[0][i]/=size;
		//	FFTtr[1][i]/=size;
		//}
		
		
		GeneralFFT.inverseTransform(FFTtr[0], FFTtr[1]);
		for(int i=0;i<FFTtr[0].length;i++)
		{
			System.out.println(FFTtr[0][i]);
		}
		System.out.print("done");
	}
*/
}