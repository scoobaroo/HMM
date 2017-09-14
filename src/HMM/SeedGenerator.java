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
	public double[] generate2(int count) {
		double[] seed = new double[count];
		double fraction = (double) 1/count;
		for(int i = 0; i< count; i++) {
			seed[i] = fraction;
		}
		seed[seed.length-1] = fraction - 0.001;
		seed[seed.length-2] = fraction - 0.002;
		seed[seed.length-3] = fraction - 0.003;
		seed[0] = fraction + 0.001;
		seed[1] = fraction + 0.002;
		seed[2] = fraction + 0.003;
		return seed;
	}
}
