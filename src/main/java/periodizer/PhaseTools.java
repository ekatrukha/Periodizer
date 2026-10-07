package periodizer;

import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Util;

public class PhaseTools
{
	public static Img<FloatType> unwrapPhasesPerPixel(final RandomAccessibleInterval< FloatType > phaseIn)
	{		
		final long[] pos = phaseIn.dimensionsAsLongArray();
		final Img<FloatType> phaseOut = new ArrayImgFactory<>(new FloatType()).create(pos);
		final RandomAccess< FloatType > raIn = phaseIn.randomAccess();
		final RandomAccess< FloatType > raOut = phaseOut.randomAccess();
		
		int numDims = pos.length;
		int timeDim = numDims - 1;
		int numFrames = ( int ) phaseIn.dimension(timeDim);
		long numSpatialPixels = 1;
		
		for (int d = 0; d < timeDim; d++) {
			numSpatialPixels *= phaseIn.dimension(d);
		}
		
		for (long p = 0; p < numSpatialPixels; p++) 
		{
			// Unravel 1D spatial index p into multi-dimensional position
			long temp = p;
			for (int d = 0; d < timeDim; d++) {
				pos[d] = temp % phaseIn.dimension(d);
				temp /= phaseIn.dimension(d);
			}
			
			pos[timeDim] = 0;
			raIn.setPosition( pos );
			float prevPhase = raIn.get().get();
			float cumCorrection = 0.0f;
			raOut.setPosition( pos );
			raOut.get().set( prevPhase );
			for (int t = 1; t < numFrames; t++) 
			{
				pos[timeDim] = t;
				raIn.setPosition( pos );
				float currentWrappedPhase = raIn.get().get();
				float diff = currentWrappedPhase - prevPhase;
				// Detect 2*PI phase jumps
                if (diff > Math.PI) {
                    cumCorrection -= 2.0f * (float) Math.PI;
                } else if (diff < -Math.PI) {
                    cumCorrection += 2.0f * (float) Math.PI;
                }
                raOut.setPosition( pos );
                raOut.get().set( currentWrappedPhase + cumCorrection );
                prevPhase = currentWrappedPhase;
			}						
		}
		
		return phaseOut;
	}
	
	public static < T extends RealType< T > & NativeType< T > > Img<T> remapPhase(final RandomAccessibleInterval< T > orig, final RandomAccessibleInterval< FloatType > unwrapPh, final double freq )
	{
		final long[] pos = orig.dimensionsAsLongArray();
		ImgFactory< T > factory = Util.getArrayOrCellImgFactory( orig, orig.getType() );
		Img< T > out = factory.create( pos );

		final RandomAccess< FloatType > raPh = unwrapPh.randomAccess();
		final RandomAccess< T > raOrig = orig.randomAccess();
		final RandomAccess< T > raOut = out.randomAccess();
		
		int numDims = pos.length;
		int timeDim = numDims - 1;
		int numFrames = ( int ) orig.dimension(timeDim);
		long numSpatialPixels = 1;
		
		for (int d = 0; d < timeDim; d++) {
			numSpatialPixels *= orig.dimension(d);
		}
		
		final double[] phaseUnr = new double[numFrames];
		final double[] phaseLin = new double[numFrames];
		final double[] timeFr = new double[numFrames];
		
		for (int t = 0; t < numFrames; t++) 
		{
			timeFr[t] = t;
		}
		for (long p = 0; p < numSpatialPixels; p++) 
		{
			// Unravel 1D spatial index p into multi-dimensional position
			long temp = p;
			for (int d = 0; d < timeDim; d++) {
				pos[d] = temp % orig.dimension(d);
				temp /= orig.dimension(d);
			}

			for (int t = 0; t < numFrames; t++) 
			{
				pos[timeDim] = t;
				raPh.setPosition( pos );
				//get phase
				phaseUnr[t] = raPh.get().get();
			}						
			//normalize phase
			double iniPh = phaseUnr[0];
			for (int t = 0; t < numFrames; t++) 
			{
				phaseUnr[t] = Math.abs( phaseUnr[t] -iniPh  );	
			}			
			for (int t = 0; t < numFrames; t++) 
			{
				phaseLin[t] = phaseUnr[0] + 2.0 * Math.PI * freq * t;
			}
			//create linear interpolation
			double[] resampledPos = LinearInterpolation.evalLinearInterp( phaseUnr, timeFr, phaseLin );
			for (int t = 0; t < numFrames; t++) 
			{
				int newT = ( int ) Math.round( resampledPos[t] );
				if(t>=0 && t<numFrames)
				{
					pos[timeDim] = t;
					raOut.setPosition( pos );
					pos[timeDim] = newT;
					raOrig.setPosition( pos );
					raOut.get().set( raOrig.get().copy() );
				}
			}
		}
		
		return out;
	}

}
