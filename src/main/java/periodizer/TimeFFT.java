package periodizer;

import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.complex.ComplexFloatType;
import net.imglib2.type.numeric.real.FloatType;

public class TimeFFT
{
	public static < T extends RealType< T > & NativeType< T > > Img<ComplexFloatType> getPerPixelFFT(final RandomAccessibleInterval< T > rai)
	{		
		long[] origDim = rai.dimensionsAsLongArray();
		int numDims = origDim.length;

		int timeDim = numDims - 1;
		int numFrames = ( int ) rai.dimension(timeDim);
	
		final Img<ComplexFloatType> fftResult = new ArrayImgFactory<>(new ComplexFloatType()).create(origDim);
		final RandomAccess< ComplexFloatType > rafftResult = fftResult.randomAccess();
		
		final GeneralFFT gfft = new GeneralFFT();
		gfft.preComputeCosSin(numFrames);
		float[] timeSeries = new float[numFrames];
		//double[] meanPower = new double[numFreqs];
		
		long numSpatialPixels = 1;
		for (int d = 0; d < timeDim; d++) {
			numSpatialPixels *= rai.dimension(d);
		}
		
		final long[] pos = new long[numDims];
		final RandomAccess<T> srcAccess = rai.randomAccess();
		for (long p = 0; p < numSpatialPixels; p++) 
		{
			// Unravel 1D spatial index p into multi-dimensional position
			long temp = p;
			for (int d = 0; d < timeDim; d++) {
				pos[d] = temp % rai.dimension(d);
				temp /= rai.dimension(d);
			}
			
			for (int t = 0; t < numFrames; t++) 
			{
				pos[timeDim] = t;
				srcAccess.setPosition(pos);
				timeSeries[t] = srcAccess.get().getRealFloat();
			}
			
			final float[][] cfft = gfft.prC_transformR( timeSeries );
			
			for (int t = 0; t < numFrames; t++) 
			{
				pos[timeDim] = t;
				rafftResult.setPosition( pos );
				float realv = cfft[0][t];
				float imv = cfft[1][t];
				
				rafftResult.get().set( new ComplexFloatType(realv,imv) );
			}
			
		}
		return fftResult;
	}
	
	public static Img<FloatType> getAmplitude(final RandomAccessibleInterval< ComplexFloatType > raiComplex)
	{
		final long[] pos = raiComplex.dimensionsAsLongArray();
		final Img<FloatType> ampIm = new ArrayImgFactory<>(new FloatType()).create(pos);
		final Cursor< ComplexFloatType > cursor = raiComplex.localizingCursor();
		final RandomAccess< FloatType > raAmp = ampIm.randomAccess();
		cursor.reset();
		while(cursor.hasNext())
		{
			cursor.fwd();
			cursor.localize( pos);
			raAmp.setPosition( pos );
			float realv = cursor.get().getRealFloat();
			float imv = cursor.get().getImaginaryFloat();
			float fAmp = ( float ) Math.sqrt( realv * realv + imv * imv );
			raAmp.get().set( fAmp );
		}
		return ampIm;
	}
	
	public static Img<FloatType> getPhase(final RandomAccessibleInterval< ComplexFloatType > raiComplex)
	{
		final long[] pos = raiComplex.dimensionsAsLongArray();
		final Img<FloatType> phaseImg = new ArrayImgFactory<>(new FloatType()).create(pos);
		final Cursor< ComplexFloatType > cursor = raiComplex.localizingCursor();
		final RandomAccess< FloatType > raPhase = phaseImg.randomAccess();
		cursor.reset();
		while(cursor.hasNext())
		{
			cursor.fwd();
			cursor.localize( pos);
			raPhase.setPosition( pos );
			float realv = cursor.get().getRealFloat();
			float imv = cursor.get().getImaginaryFloat();
			float phase = ( float ) Math.atan2( imv, realv );
			raPhase.get().set( phase );
		}
		return phaseImg;
	}
	
	public static RandomAccessibleInterval< ComplexFloatType > getInverseFFT(final RandomAccessibleInterval< ComplexFloatType > raiFFT)
	{
		final long[] origDim = raiFFT.dimensionsAsLongArray();
		int numDims = origDim.length;

		final Img<ComplexFloatType> out = new ArrayImgFactory<>(new ComplexFloatType()).create(origDim);

		int timeDim = numDims - 1;
		int numFrames = ( int ) raiFFT.dimension(timeDim);
		final float[][] realim = new float[2][numFrames];
		long numSpatialPixels = 1;
		for (int d = 0; d < timeDim; d++) {
			numSpatialPixels *= raiFFT.dimension(d);
		}
		
		final long[] pos = new long[numDims];
		final RandomAccess< ComplexFloatType > fftAccess = raiFFT.randomAccess();
		final RandomAccess< ComplexFloatType > raOut = out.randomAccess();
		final GeneralFFT gfft = new GeneralFFT();
		gfft.preComputeCosSin(numFrames);
		
		for (long p = 0; p < numSpatialPixels; p++) 
		{
			// Unravel 1D spatial index p into multi-dimensional position
			long temp = p;
			for (int d = 0; d < timeDim; d++) {
				pos[d] = temp % raiFFT.dimension(d);
				temp /= raiFFT.dimension(d);
			}
			
			for (int t = 0; t < numFrames; t++) 
			{
				pos[timeDim] = t;
				fftAccess.setPosition(pos);
				realim[0][t] = fftAccess.get().getRealFloat();
				realim[1][t] = fftAccess.get().getImaginaryFloat();			
			}
			gfft.prC_inverseTransform( realim[0], realim[1] );
			for (int t = 0; t < numFrames; t++) 
			{
				pos[timeDim] = t;
				raOut.setPosition( pos );
				raOut.get().set( new ComplexFloatType(realim[0][t],realim[1][t]) );
			}
		}
		return out;
	}
	
	public static Img<ComplexFloatType> gaussPositiveBandPass(final RandomAccessibleInterval< ComplexFloatType > raiFFT, final float fCenter, final float fSD)
	{
		final long[] pos = raiFFT.dimensionsAsLongArray();
		int numDims = pos.length;	
		
		int timeDim = numDims - 1;
		int numFrames = ( int ) raiFFT.dimension(timeDim);
		int nHalf = ( int ) ( numFrames * 0.5 ); 
		final Img<ComplexFloatType> fftBP = new ArrayImgFactory<>(new ComplexFloatType()).create(pos);
		final Cursor< ComplexFloatType > origC = raiFFT.localizingCursor();
		final Cursor< ComplexFloatType > outC = fftBP.localizingCursor();
		origC.reset();
		outC.reset();
		while(origC.hasNext())
		{
			origC.fwd();
			outC.fwd();
			origC.localize( pos );
			ComplexFloatType val = origC.get().copy();
			float fT = pos[timeDim];
			float coeff = 0.0f;
			if(fT < nHalf)
			{
				coeff = ( float ) ( 2.0f * Math.exp( -0.5 * (fT - fCenter)*(fT - fCenter)/(fSD*fSD) ) );

			}
			val.mul( coeff );
			outC.get().set( val );
		}
		return fftBP;
	}
	

}
