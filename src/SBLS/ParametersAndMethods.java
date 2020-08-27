package SBLS;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.Arrays;
import java.util.BitSet;

import org.apache.commons.math3.util.FastMath;


/**
 * Created by aviadpi on 23/11/2016.
 */
public abstract class ParametersAndMethods {
	public static final double EPS = 1e-10;
	
	public static final int iIndex          = 0;
	public static final int iTime           = 1;
	public static final int iValue          = 2;
	public static final int iError          = 3;
	public static final int iWeight         = 4;
	public static final int lineSize        = 5;

	// Duration parameters
    public static final double minDuration = 0.02;
    public static final double maxDuration = 0.323; // 0.15 * 10^(1/3)
    public static final double durationFactor = 0.15; // This factor times the Period^1/3 in days gives the maximal duration in days
    public static final int minPointsInTransit = 3;
    
 // Frequency parameters
 	public static final double fStep = 1e-4;
 	public static final double fStart = 0.01;
 	public static final double fEnd = 2.0 + EPS;
 	
    public static double[][] foldAndSortPhases(double[][] Data, double[][] phasesIndices, double period, double zeroTime, boolean usePrevData) throws IndexOutOfBoundsException {
        int N = Data.length;
        int j;
        
        // Folding
        try {
            for (int i = 0; i < N; i++) {
        	    j = usePrevData ? (int) phasesIndices[i][iIndex] : i;
                double phase 	= dMod(Data[j][iTime] - zeroTime , period)/period;
                phasesIndices[i]= new double[] {j, phase} ;
            }
        } catch (IndexOutOfBoundsException e){
            e.printStackTrace();
        }

        // Sorting
        Arrays.parallelSort(phasesIndices, (a, b) -> Double.compare(a[iTime], b[iTime]));
        return phasesIndices;
    }
    
    /**
	 * Modulo corrected for negative dividend values
	 * @param a The dividend
	 * @param n The divisor
	 * @return The remainder of a mod n, between 0 and n, like it should be.
	 */
	public static double dMod(double a, double n) {
		double rem = a%n;
		return (rem>=0. ? rem : rem+n);
	}
	public static int iMod(int a, int n) {
		int rem = a%n;
		return (rem>=0. ? rem : rem+n);
	}
	
	/**
	 * Centers the time series so that the weighted average is zero, with the option of normalizing the weights.
	 * @param data
	 * @return
	 */
	public static double[][] centerLightCurveData(double[][] data, boolean normalizeWeights){
		int N = data.length;
		double[][] dataClone = new double[N][lineSize];
		double weightedSum = 0.,sumOfWeights=0.;

		// Clone, calculate weights and the weighted sum
		for (int i=0; i<N; i++){
			dataClone[i] = data[i].clone();
			dataClone[i][iWeight] = FastMath.pow(data[i][iError], -2);
			sumOfWeights += dataClone[i][iWeight];
			weightedSum += dataClone[i][iValue]*dataClone[i][iWeight];
		}
		weightedSum /= sumOfWeights;
		for (int i=0; i<N; i++){
			if (normalizeWeights) dataClone[i][iWeight] /= sumOfWeights;
			dataClone[i][iValue] -= weightedSum;
		}
		return dataClone;
	}

	public static double[][] createDataFromTimesValuesErrors(double[] times, double[] values, double[] errors){
		if (times.length != values.length) return null;
		int N = times.length;
		double[][] data = new double[N][lineSize];
		for (int i=0; i< N; i++){
			data[i][iTime] = times[i];
			data[i][iValue] = values[i];
			data[i][iError] = errors[i];
			data[i][iIndex] = i;
		}
		return data;
	}

	public static double roundInverse(double f, double end, double step){
		int scale = (int) (2* FastMath.log10(end/step));
		BigDecimal dec = new BigDecimal(1.);
		try {
			dec = new BigDecimal(1. / f);
		} catch (Exception e) {
			System.out.println(String.format("Failed to calculate 1/%.5f",f));
			e.printStackTrace();
		}
		double p = dec.setScale(scale, RoundingMode.HALF_EVEN).doubleValue();
		return p;	
	}
	
	public static double[] createFrequencyGrid(double start, double end, double step){
		int nFrequencies = (int) (FastMath.abs(end-start)/step) + 1;
		double[] frequencies = new double[nFrequencies];
		int scale = (int) (2* FastMath.log10(nFrequencies));
		for (int i = 0; i < nFrequencies; i++) {
			frequencies[i] = new BigDecimal(start + i * step).setScale(scale,
					RoundingMode.HALF_EVEN).doubleValue();
		}
		return frequencies;
	}

	public static double[] createFrequencyGrid(double step){
		return createFrequencyGrid(fStart, fEnd, step);
	}
	
	public static int[] convertBitSetToIndices(BitSet bs) {
		int nIT = bs.cardinality();
		int[] indices = new int[nIT];
		int j=0;
		for (int i=0; i<bs.size(); i++){
			if (bs.get(i)) {
				indices[j] = i;
				j++;
			}
		}
		return indices;
	}
}
