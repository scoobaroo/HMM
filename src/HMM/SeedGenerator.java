package HMM;

public class SeedGenerator {

	public SeedGenerator() {

	}
	
	public double[] generate(int count) {
		double[] seed = new double[count];
		double total = 0.0;
		for(int i = 0 ; i <count ; i++) {
			double num = Math.random() * 100;
			seed[i] = num;
			total += num;
		}
		double[] probSeed = new double[count];
		for(int i = 0 ; i <count ; i++) {
			probSeed[i] = seed[i] / total;
		}
		return probSeed;
	}
}
