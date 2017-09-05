package HMM;
import java.util.Hashtable;
import java.util.Vector;
import org.paukov.combinatorics.Factory;
import org.paukov.combinatorics.Generator;
import org.paukov.combinatorics.ICombinatoricsVector;
import org.paukov.combinatorics.util.Util;

public class HMM {
	double[][] directComputationResult = new double[81][16]; 
	public int numStates;
	public int sigmaSize;
	public double pi[];
	public double A[][];
	public double B[][];

//	public Hashtable<Double, Double> pi;
    private Vector<Double> states;
    private Vector<Double> observations;

    private Hashtable<Pair<Double, Double>, Double> transitionMatrix;
    private Hashtable<Pair<Double, Double>, Double> emissionMatrix;


	public HMM(int numStates, int sigmaSize) {
		this.numStates = numStates;
		this.sigmaSize = sigmaSize;
		double[] pi = new double[numStates];
		double[][] A = new double[numStates][numStates];
		double[][] B = new double[numStates][sigmaSize];
	}
	
//    public void setTransitionMatrix(Hashtable<Pair<Double, Double>, Double> transitionMatrix) {
//        this.transitionMatrix = transitionMatrix;
//    }
//
//    public Hashtable<Pair<Double, Double>, Double> getEmissionMatrix() {
//        return emissionMatrix;
//    }
//
//    public void setEmissionMatrix(Hashtable<Pair<Double, Double>, Double> emissionMatrix) {
//        this.emissionMatrix = emissionMatrix;
//    }
//
//    public Double getTransitionValue(Double firstState, Double secondState) {
//        return this.transitionMatrix.get(new Pair<Double, Double>(firstState, secondState));
//    }
//
//    public Double getEmissionValue(Double state, Double observation) {
//        return this.emissionMatrix.get(new Pair<Double, Double>(state, observation));
//    }
//
//    public void setInitialProbabilities(Hashtable<Double, Double> pi) {
//        this.pi = pi;
//    }
//	
//    public Hashtable<Double, Double> getInitialProbabilities() {
//        return this.pi;
//    }
//
//    public Double getInitialProbability(double state) {
//        return this.pi.get(state);
//    }

	@SuppressWarnings("unchecked")
	public static void main(String[] args) throws Exception {
		HMM2 hmm2 = new HMM2(2,3);
		hmm2.a = new double[][]{{0.7,0.3},{0.4,0.6}};
		hmm2.b = new double[][]{{0.1,0.4,0.5},{0.7,0.2,0.1}};
		hmm2.pi = new double[]{0.0, 1.0};
		double[] states = new double[]{0.0,0.0,0.0,0.0};
		double[] observations = new double[]{1.0,0.0,0.0,2.0};
	    ICombinatoricsVector<Double> observationVector = Factory.createVector(new Double[] {0.0,1.0,2.0});
	    ICombinatoricsVector<Double> stateVector = Factory.createVector(new Double[]{0.0,1.0});
	    Generator<Double> observationGen = Factory.createPermutationWithRepetitionGenerator(observationVector, 4);
	    Generator<Double> stateGen = Factory.createPermutationWithRepetitionGenerator(stateVector, 4);
	    for (ICombinatoricsVector<Double> ob : observationGen){
	    	for (ICombinatoricsVector<Double> st : stateGen){
	    		double result = hmm2.evaluateUsingBruteForce(st.getVector(), ob.getVector());
	    		System.out.println("results from using observation sequence " + ob.getVector() +
	    							" and state sequence " + st.getVector() +": " + result);
	    	}
	    }
//		hmm2.evaluateUsingBruteForce(states, observations);
	}

}
