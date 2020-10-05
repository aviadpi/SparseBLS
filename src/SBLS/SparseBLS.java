package SBLS;

import SBLS.ParametersAndMethods;
import BLS_Results.SBLS_Results;

import java.util.BitSet;
import org.apache.commons.math3.util.FastMath;


public class SparseBLS extends ParametersAndMethods {

	public static final String description = "Sparse BLS";
	
	/**
	 * My (A. Panahi) version of the BLS
	 * @param DataOriginal - Array of "LineArray"s
	 * @param nFrequencies - Number of test frequencies
	 * @param isFlux - true=flux, false=magnitude
	 * @param usePrevOrder - for faster sorting (Nearly-sorted array)
	 * @return
	 * TransitSearchResults with:
	 * Period (days), Phase (days), Duration (days), Depth (magnitude), OOT level (magnitude) etc.
	 */
	private static SBLS_Results SparseBLSmethod(double[][] DataOriginal, double fStep, boolean isFlux, boolean usePrevOrder){
		// Variables & constants declarations
		int j2, i2, transitDurationInPoints, j1ln, i1ln, i2rn, i1best=0, i2best=0;
		double P, s, r, sMinPrev, rMinPrev, SR, SRmax, TOP_SRmax, transitDurationInTime, maxTransitDuration=0., ingress, egress;

		// Create centered lightcurve with normalized weights
		boolean normalizeWeights = true;
		double[][] Data =  centerLightCurveData(DataOriginal, normalizeWeights);
		int N = Data.length;
		double[][] phasesIndices = new double[N][2];

		// Calculate some sums
		double timesAvg = 0.;

		for(int i=0; i<N; i++) {
			phasesIndices[i][iIndex] = i;
			phasesIndices[i][iTime]  = Data[i][iTime];

			timesAvg += DataOriginal[i][iTime];
		}
		timesAvg /= N;
		
		SRmax = 0.0; TOP_SRmax = 0.0;
		double[] frequencies = createFrequencyGrid(fStep);
		int nf = frequencies.length;

		double[] periods 					= new double[nf];
		double[] MaxSignalResidues 			= new double[nf];
		BitSet proposedInTransitIndices 	= new BitSet(N);

		for (int i = 0; i < nf; i++) {
			periods[i] = roundInverse(frequencies[i], fEnd, fStep);
		}

		// Main Loop over all periods
		for (int iMain=0; iMain<nf; iMain++) {
			P = periods[iMain];

			// Fold it up!
			try {
				phasesIndices = foldAndSortPhases(Data, phasesIndices, P, timesAvg, usePrevOrder);
			} catch (IndexOutOfBoundsException e) {
				e.printStackTrace();
			}
			int i2MinPrev = 0;
			s = r = sMinPrev = rMinPrev = 0.0;

			// i1 Loop - Iterate over all possible i1 values
			for (int i1 = 0; i1 < N; i1++) {
				int j1 	= (int) phasesIndices[i1][iIndex];

				// Load previous indices
				i2 = i2MinPrev;
				j2 = (int) phasesIndices[i2][iIndex];
				s = sMinPrev;
				r = rMinPrev;
				i1ln = iMod(i1-1,N);
				j1ln = (int) phasesIndices[i1ln][iIndex];
				
				if (i1 == 0) {
					ingress =  ((phasesIndices[i1ln][iTime]-1.) + phasesIndices[i1][iTime])/2.;
				} else {
					// Remove values from each end
					s -= (Data[j1ln][iValue] * Data[j1ln][iWeight]);
					r -= Data[j1ln][iWeight];
					s -= (Data[j2][iValue] * Data[j2][iWeight]);
					r -= Data[j2][iWeight];
					ingress = (phasesIndices[i1ln][iTime] + phasesIndices[i1][iTime])/2.;
				}
								
				// i2 Loop - Look for relevant i2 indices
				while (true) {
					// Calculate transit duration in time and in points
					transitDurationInPoints = i2 - i1 + 1;
					// Adjust in case transit falls between folds
					if (i2 < i1) transitDurationInPoints += N;				
					// Stoppage conditions
					if ((transitDurationInPoints < minPointsInTransit)){
						// Sum s,r
						s += (Data[j2][iValue] * Data[j2][iWeight]);
						r += Data[j2][iWeight];
						
						// Save data
						i2MinPrev = i2;
						sMinPrev = s;
						rMinPrev = r;
						
						i2 = (i2 + 1) % N;
						j2 = (int) phasesIndices[i2][iIndex];
						continue;
					}
					
					i2rn = (i2 + 1) % N; // j2 right neighbor
					
					egress  = (phasesIndices[i2][iTime] + phasesIndices[i2rn][iTime]) / 2.;
					if (i2rn<i2) egress += 0.5;
					
					transitDurationInTime = dMod(egress-ingress, 1.)*P;


					// Stoppage condition - Keplerian duration too long
					maxTransitDuration = ParametersAndMethods.durationFactor * FastMath.pow(P,1./3.) + EPS; // Maximal transit duration up to 0.15*P=^0.333 days
					if (transitDurationInTime > maxTransitDuration) break;

					s += (Data[j2][iValue] * Data[j2][iWeight]);
					r += Data[j2][iWeight];

					// For each i1,i2 calculate SR and update if maximal
					double depth = s / (r * (r - 1));
					SR = (s * s / (r * (1 - r))); // This is actually the square, the root was omitted.
					boolean bigEnoughTransit = isFlux ? (depth>0) : (depth<0);

					if ((SR > (SRmax + EPS)) && bigEnoughTransit){
						SRmax = SR;
						i1best = i1;
						i2best = i2;
					}
					// Allow transit to fall between folds
					i2 = (i2 + 1) % N;
					j2 = (int) phasesIndices[i2][iIndex];
				} // end of i2
			} // end of i1

			// Save maximum SR value
			MaxSignalResidues[iMain] = SRmax;

			if (SRmax > TOP_SRmax) {
				TOP_SRmax = SRmax;
				
				// Update InTransit indices
				int iCurrent;
				proposedInTransitIndices.clear();
				if (i1best<i2best){
					for (int q=i1best; q <= i2best; q++){
						iCurrent = (int) phasesIndices[q][iIndex];
						proposedInTransitIndices.set(iCurrent);
					}
				}else if (i1best>i2best){
					for (int q=0; q <= i2best; q++){
						iCurrent = (int) phasesIndices[q][iIndex];
						proposedInTransitIndices.set(iCurrent);
					}
					for (int q=i1best; q < N; q++){
						iCurrent = (int) phasesIndices[q][iIndex];
						proposedInTransitIndices.set(iCurrent);
					}
				}
			}
			// Reset the max variables
			SRmax = 0.0;
		} // end of Periods loop
		int[] foundInTransitIndices = convertBitSetToIndices(proposedInTransitIndices);
		SBLS_Results result = new SBLS_Results(frequencies, MaxSignalResidues, foundInTransitIndices);
		return result;
	}

	public static SBLS_Results SparseBLSmethod(double[] times, double[] values, double[] errors){
		double[][] Data = createDataFromTimesValuesErrors(times, values, errors);
		boolean isFlux = false;
		boolean usePrevOrder = true;
		return SparseBLSmethod(Data, fStep, isFlux, usePrevOrder);
	}
	
	
}
