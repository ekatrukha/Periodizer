package periodizer;



import net.imglib2.RandomAccess;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.complex.ComplexFloatType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Cast;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.plugin.PlugIn;


public class Periodizer implements PlugIn 
{

	@Override
	public void run( String arg )
	{
    }

	/**
	 * Main method for debugging.
	 *
	 * For debugging, it is convenient to have a method that starts ImageJ, loads
	 * an image and calls the plugin, e.g. after setting breakpoints.
	 *
	 * @param args unused
	 */
	public static <R extends RealType< R >> void main(String[] args) throws Exception 
	{
			new ImageJ();
			
			// 1. Load image
			ImagePlus image = IJ.openImage("/home/eugene/Desktop/projects/Periodizer/05.tif");
			image.show();

			Img<R> img_in = Cast.unchecked(ImageJFunctions.wrap(image));

			IntervalView< R > rai = Views.permute( img_in, 0, 1 );
			// Dimensions: assuming [X, Y, T] order
			long numDims = rai.numDimensions();
			
			long[] origDim = rai.dimensionsAsLongArray();
			//make time dimension the last one'
			//int timeDim = (int) numDims - 1; // Assuming time axis is the last dimension (2)
			int timeDim = (int) numDims - 1;
			int numFrames = ( int ) rai.dimension(timeDim);
		
			// Number of positive frequency bins (up to Nyquist limit)
			//int numFreqs = numFrames / 2 + 1;
//			long[] specDims = new long[( int ) numDims];
//			for (int d = 0; d < timeDim; d++) {
//				specDims[d] = rai.dimension(d);
//			}
//			specDims[timeDim] = numFreqs;
			
			
			Img<FloatType> powerSpectrum = new ArrayImgFactory<>(new FloatType()).create(origDim);		
			Img<ComplexFloatType> fftResult = new ArrayImgFactory<>(new ComplexFloatType()).create(origDim);
			RandomAccess< ComplexFloatType > rafftResult = fftResult.randomAccess();
			RandomAccess< FloatType > raPS = powerSpectrum.randomAccess();
			
			GeneralFFT gfft = new GeneralFFT();
			gfft.preComputeCosSin(numFrames);
			float[] timeSeries = new float[numFrames];
			//double[] meanPower = new double[numFreqs];
			
			long numSpatialPixels = 1;
			for (int d = 0; d < timeDim; d++) {
				numSpatialPixels *= rai.dimension(d);
			}
			
			final long[] pos = new long[(int) numDims];
			RandomAccess<R> srcAccess = rai.randomAccess();
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
					raPS.setPosition( pos );
					float power = ( float ) Math.sqrt( realv* realv + imv * imv );
					raPS.get().set( power );
				}
				
			}
			ImageJFunctions.show( powerSpectrum, "PS" );
			ImageJFunctions.show( fftResult, "FFT" );

		}
}




