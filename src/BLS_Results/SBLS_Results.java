package BLS_Results;

public class SBLS_Results{
	
	private double[] frequencies=null;
	private double[] powers=null;
	private int[] foundInTransitIndices=null;
	
    // Constructors
    public SBLS_Results(double[] frequencies, double[] powers, int[] foundInTransitIndices) {
		this.frequencies = frequencies;
		this.powers = powers;
		this.foundInTransitIndices = foundInTransitIndices;
	}

	public double[] getFrequencies() {
		return frequencies;
	}

	public double[] getPowers() {
		return powers;
	}

	public int[] getFoundInTransitIndices() {
		return foundInTransitIndices;
	}

	public void setFrequencies(double[] frequencies) {
		this.frequencies = frequencies;
	}

	public void setPowers(double[] powers) {
		this.powers = powers;
	}

	public void setFoundInTransitIndices(int[] foundInTransitIndices) {
		this.foundInTransitIndices = foundInTransitIndices;
	}
	
	

	
}