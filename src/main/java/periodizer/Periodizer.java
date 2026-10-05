package periodizer;



import java.util.List;

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
			ImagePlus image = IJ.openImage("/home/eugene/Desktop/projects/Periodizer/05rot.tif");
			image.show();

			Img<R> img_in = Cast.unchecked(ImageJFunctions.wrap(image));

			//IntervalView< R > rai = Views.permute( img_in, 0, 1 );
			//Img< ComplexFloatType > fftResult = TimeFFT.getPerPixelFFT( rai );
			Img< ComplexFloatType > fftResult = TimeFFT.getPerPixelFFT( img_in );

			Img< FloatType > powerSpectrum = TimeFFT.getPowerSpectrum( fftResult );
			ImageJFunctions.show( powerSpectrum, "PS" );
			//ImageJFunctions.show( fftResult, "FFT" );
			
			Img< ComplexFloatType > bpFFT = TimeFFT.gaussPositiveBandPass( fftResult, 16, 2 );
			Img< FloatType > powerSpectrumBP = TimeFFT.getPowerSpectrum( bpFFT );
			ImageJFunctions.show( powerSpectrumBP, "bpPS" );

			//ImageJFunctions.show( bpFFT, "bpFFT" );
			List< Img< FloatType > > inv = TimeFFT.getInverseFFT( bpFFT );
			
//			for(Img< FloatType > im : inv)
//			{
//				ImageJFunctions.show(im); 
//			}

		}
}




