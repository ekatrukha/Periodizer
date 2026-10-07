package periodizer;



import java.util.List;

import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
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
	public static < T extends RealType< T > & NativeType< T > > void main(String[] args) throws Exception 
	{
			new ImageJ();
			
			// 1. Load image
			ImagePlus image = IJ.openImage("/home/eugene/Desktop/projects/Periodizer/05rot.tif");
			//ImagePlus image = IJ.openImage("/home/eugene/Desktop/projects/Periodizer/01_Concatenated_crop_descewed-05_TL.tif");
			//ImagePlus image = IJ.openImage("/home/eugene/Desktop/projects/Periodizer/crop.tif");		
			image.show();

			Img<T> img_in = Cast.unchecked(ImageJFunctions.wrap(image));
			
			final long[] dims = img_in.dimensionsAsLongArray();

			int numDims = dims.length;
			int timeDim = numDims - 1;
			int numFrames = ( int ) img_in.dimension(timeDim);
			Img< ComplexFloatType > fftResult = TimeFFT.getPerPixelFFT( img_in );
			int nVal = 17;
			Img< FloatType > powerSpectrum = TimeFFT.getAmplitude( fftResult );
			ImageJFunctions.show( powerSpectrum, "Amplitude" );

			Img< ComplexFloatType > bpFFT = TimeFFT.gaussPositiveBandPass( fftResult, nVal + 1, 2 );

			RandomAccessibleInterval< ComplexFloatType > inv = TimeFFT.getInverseFFT( bpFFT );

			Img< FloatType > wrappedPhaseImg = TimeFFT.getPhase( inv );
			ImageJFunctions.show( wrappedPhaseImg, "phaseswrapped" );
			Img< FloatType > unwrappedPhase = PhaseTools.unwrapPhasesPerPixel( wrappedPhaseImg );
			Img< T > remapped = PhaseTools.remapPhase(img_in, unwrappedPhase, 16/((double)numFrames));
			ImageJFunctions.show(remapped, "remapped");

			ImageJFunctions.show(  unwrappedPhase, "bpFFTPhaseUnW" );
			


		}
}




