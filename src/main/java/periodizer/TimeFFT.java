package periodizer;

import java.util.ArrayList;
import java.util.List;

import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.complex.ComplexFloatType;
import net.imglib2.type.numeric.real.FloatType;

public class TimeFFT
{
	public static <R extends RealType< R >> Img<ComplexFloatType> getPerPixelFFT(final RandomAccessibleInterval< R > rai)
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
		final RandomAccess<R> srcAccess = rai.randomAccess();
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
	
	public static <R extends RealType< R >> Img<FloatType> getPowerSpectrum(final RandomAccessibleInterval< ComplexFloatType > raiFFT)
	{
		final long[] pos = raiFFT.dimensionsAsLongArray();
		final Img<FloatType> powerSpectrum = new ArrayImgFactory<>(new FloatType()).create(pos);
		final Cursor< ComplexFloatType > cursor = raiFFT.localizingCursor();
		final RandomAccess< FloatType > raPS = powerSpectrum.randomAccess();
		cursor.reset();
		while(cursor.hasNext())
		{
			cursor.fwd();
			cursor.localize( pos);
			raPS.setPosition( pos );
			float realv = cursor.get().getRealFloat();
			float imv = cursor.get().getImaginaryFloat();
			float power = ( float ) Math.sqrt( realv* realv + imv * imv );
			raPS.get().set( power );
		}
		return powerSpectrum;
	}
	
	public static <R extends RealType< R >> List<Img<FloatType>> getInverseFFT(final RandomAccessibleInterval< ComplexFloatType > raiFFT)
	{
		final long[] origDim = raiFFT.dimensionsAsLongArray();
		int numDims = origDim.length;

		
		final ArrayList<Img<FloatType>> out = new ArrayList<>();
		final ArrayList<RandomAccess< FloatType >> cursors = new ArrayList<>();
		
		for(int i = 0; i < 2; i++)
		{
			out.add( new ArrayImgFactory<>(new FloatType()).create(origDim) );
			cursors.add( out.get( i ).randomAccess());
		}
		int timeDim = numDims - 1;
		int numFrames = ( int ) raiFFT.dimension(timeDim);
		final float[][] realim = new float[2][numFrames];
		long numSpatialPixels = 1;
		for (int d = 0; d < timeDim; d++) {
			numSpatialPixels *= raiFFT.dimension(d);
		}
		
		final long[] pos = new long[numDims];
		RandomAccess< ComplexFloatType > fftAccess = raiFFT.randomAccess();
		
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
				for(int i = 0; i < 2; i++)
				{
					cursors.get( i ).setPosition( pos );
					cursors.get( i ).get().set( realim[i][t] );
				}
			}
		}
		
		
		return out;
	}
	
	public static  Img<ComplexFloatType> gaussPositiveBandPass(final RandomAccessibleInterval< ComplexFloatType > raiFFT, final float fCenter, final float fSD)
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
