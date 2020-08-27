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
		int i2, j2, transitDurationInPoints, i1ln, j1ln, j2rn, j1best=0, j2best=0;
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
			int j2MinPrev = 0;
			s = r = sMinPrev = rMinPrev = 0.0;

			// i1 Loop - Iterate over all possible i1 values
			for (int j1 = 0; j1 < N; j1++) {
				int i1 	= (int) phasesIndices[j1][iIndex];

				// Load previous indices
				j2 = j2MinPrev;
				i2 = (int) phasesIndices[j2][iIndex];
				s = sMinPrev;
				r = rMinPrev;
				j1ln = iMod(j1-1,N);
				i1ln = (int) phasesIndices[j1ln][iIndex];
				
				if (j1 == 0) {
					ingress =  ((phasesIndices[j1ln][iTime]-1.) + phasesIndices[j1][iTime])/2.;
				} else {
					// Remove values from each end
					s -= (Data[i1ln][iValue] * Data[i1ln][iWeight]);
					r -= Data[i1ln][iWeight];
					s -= (Data[i2][iValue] * Data[i2][iWeight]);
					r -= Data[i2][iWeight];
					ingress = (phasesIndices[j1ln][iTime] + phasesIndices[j1][iTime])/2.;
				}
								
				// i2 Loop - Look for relevant i2 indices
				while (true) {
					// Calculate transit duration in time and in points
					j2rn = (j2 + 1) % N; // j2 right neighbor
					
					egress  = (phasesIndices[j2][iTime] + phasesIndices[j2rn][iTime]) / 2.;
					if (j2rn<j2) egress += 0.5;
					
					transitDurationInTime = dMod(egress-ingress, 1.)*P;
					transitDurationInPoints = j2 - j1 + 1;

					// Adjust in case transit falls between folds
					if (j2 < j1) transitDurationInPoints += N;				

					// Stoppage conditions
					if ((transitDurationInPoints < minPointsInTransit)){
						// Sum s,r
						s += (Data[i2][iValue] * Data[i2][iWeight]);
						r += Data[i2][iWeight];

						// Save data
						j2MinPrev = j2;
						sMinPrev = s;
						rMinPrev = r;

						j2 = (j2 + 1) % N;
						i2 = (int) phasesIndices[j2][iIndex];
						continue;
					}

					// Stoppage condition - Keplerian duration too long
					maxTransitDuration = ParametersAndMethods.durationFactor * FastMath.pow(P,1./3.) + EPS; // Maximal transit duration up to 0.15*P=^0.333 days
					if (transitDurationInTime > maxTransitDuration) break;

					s += (Data[i2][iValue] * Data[i2][iWeight]);
					r += Data[i2][iWeight];

					// For each i1,i2 calculate SR and update if maximal
					double depth = s / (r * (r - 1));
					SR = (s * s / (r * (1 - r))); // This is actually the square, the root was omitted.
					boolean bigEnoughTransit = isFlux ? (depth>0) : (depth<0);

					if ((SR > (SRmax + EPS)) && bigEnoughTransit){
						SRmax = SR;
						j1best = j1;
						j2best = j2;
					}
					// Allow transit to fall between folds
					j2 = (j2 + 1) % N;
					i2 = (int) phasesIndices[j2][iIndex];
				} // end of i2
			} // end of i1

			// Save maximum SR value
			MaxSignalResidues[iMain] = SRmax;

			if (SRmax > TOP_SRmax) {
				TOP_SRmax = SRmax;
				
				// Update InTransit indices
				int iCurrent;
				proposedInTransitIndices.clear();
				if (j1best<j2best){
					for (int j=j1best; j <= j2best; j++){
						iCurrent = (int) phasesIndices[j][iIndex];
						proposedInTransitIndices.set(iCurrent);
					}
				}else if (j1best>j2best){
					for (int j=0; j <= j2best; j++){
						iCurrent = (int) phasesIndices[j][iIndex];
						proposedInTransitIndices.set(iCurrent);
					}
					for (int j=j1best; j < N; j++){
						iCurrent = (int) phasesIndices[j][iIndex];
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
