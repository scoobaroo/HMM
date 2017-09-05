package HMM;
import java.util.Hashtable;
import java.util.Vector;

public class HMM {
	@supresswarnings 'unused';
	double[][] directComputationResult = new double[81][16]; 
	
//	public int numStates;
//	public int sigmaSize;
//	public double pi[];
//	public double A[][];
//	public double B[][];

	public Hashtable<Double, Double> pi;
    private Vector<Double> states;
    private Vector<Double> observations;

    private Hashtable<Pair<Double, Double>, Double> transitionMatrix;
    private Hashtable<Pair<Double, Double>, Double> emissionMatrix;


//	public HMM(int numStates, int sigmaSize) {
//		this.numStates = numStates;
//		this.sigmaSize = sigmaSize;
//		double[] pi = new double[numStates];
//		double[][] A = new double[numStates][numStates];
//		double[][] B = new double[numStates][sigmaSize];
//	}
//	
    public void setTransitionMatrix(Hashtable<Pair<Double, Double>, Double> transitionMatrix) {
        this.transitionMatrix = transitionMatrix;
    }

    public Hashtable<Pair<Double, Double>, Double> getEmissionMatrix() {
        return emissionMatrix;
    }

    public void setEmissionMatrix(Hashtable<Pair<Double, Double>, Double> emissionMatrix) {
        this.emissionMatrix = emissionMatrix;
    }

    public Double getTransitionValue(Double firstState, Double secondState) {
        return this.transitionMatrix.get(new Pair<Double, Double>(firstState, secondState));
    }

    public Double getEmissionValue(Double state, Double observation) {
        return this.emissionMatrix.get(new Pair<Double, Double>(state, observation));
    }

    public void setInitialProbabilities(Hashtable<Double, Double> pi) {
        this.pi = pi;
    }
	
    public Hashtable<Double, Double> getInitialProbabilities() {
        return this.pi;
    }

    public Double getInitialProbability(double state) {
        return this.pi.get(state);
    }
    
    @SuppressWarnings("null")
	public double evaluateUsingBruteForce(Vector<Double> states, Vector<Double> observations) throws Exception {
        if (states.size() != observations.size())
            throw new Exception("States and Observations must be at a same size!");

        double previousState = (Double) null;
        double probability = 0.0;
        double result = 0.0;

        for (int i = 0; i < states.size(); i++) {
            probability = this.getInitialProbability(states.get(i));
            previousState = (Double) null;
            for (int j = 0; j < observations.size(); j++) {
                double emissionValue = this.getEmissionValue(states.get(j), observations.get(j));
                double transitionValue = 0.0;
                if (j != 0) {
                    transitionValue += this.getTransitionValue(previousState, states.get(j));
                    probability *= transitionValue * emissionValue;
                }
                previousState = states.get(j);
            }
            result += probability;
        }

        return result;
    }

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		HMM hmm = new HMM(2,3);
		hmm.A = new double[][]{{0.7,0.3},{0.4,0.6}};
		hmm.B = new double[][]{{0.1,0.4,0.5},{0.7,0.2,0.1}};
		hmm.pi = new double[]{0.0, 1.0};
	}

}
