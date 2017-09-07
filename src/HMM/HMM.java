package HMM;
import java.io.BufferedReader;
import java.io.FileReader;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;
import org.paukov.combinatorics.Factory;
import org.paukov.combinatorics.Generator;
import org.paukov.combinatorics.ICombinatoricsVector;
import org.paukov.combinatorics.util.Util;

@SuppressWarnings("unused")
public class HMM { 
	  public int N;
	  public int M;
	  public double pi[];
	  public double a[][];
	  public double b[][];

	public HMM(int N, int M) {
	    this.N = N;
	    this.M = M;
	    pi = new double[N];
	    a = new double[N][N];
	    b = new double[N][M];
	  }
	
	  double[] c2 = new double[50000];
	  public double[][] alphaPass(double[] o){
		  int T = o.length;
		  double[] c = new double[T];
		  c[0]=0.0;
		  double[][] alpha= new double[N][T];
		  //initialize
		  for (int i = 0; i < N; i++){
		    alpha[i][0] = pi[i] * b[i][(int) o[0]];
		    c[0] += alpha[i][0];
		  }
		  c[0] = 1/c[0];
		  for(int i=0; i< N-1; i++){
		  	alpha[i][0] = c[0] * alpha[i][0];
		  }
		  for(int t = 1; t<T; t++){
			  c[t] = 0;
			  for(int i = 0; i< N; i++){
				  alpha[i][t] = 0;
				  for(int j=0; j<N; j++){
					  alpha[i][t] += alpha[j][t-1] * a[j][i];
				  }
				  alpha[i][t] *= b[i][(int) o[t]];
				  c[t] += alpha[i][t];
			  }
			  c[t] = 1/c[t];
			  for(int i = 0; i<N; i++){
				  alpha[i][t] *= c[t];
			  }
		  }
		  c2=c;
		  return alpha;
	  }
	  
	  public double[][] betaPass(double[] o) {
	    int T = o.length;
	    double[][] beta = new double[N][T];  
	    double[] c = new double[T];
	    /* initialization (time 0) */
	    for (int i = 0; i < N; i++){
	    		beta[i][T-1] = c2[T-1];
//	    		beta[i][T-1] = 1;
	    }
	    for(int t=T-2; t>=0; t--){
		    	for(int i = 0; i<N; i++){
		        	beta[i][t] = 0;
		        	for(int j =0; j<N; j++){
		        		beta[i][t] += a[i][j]*beta[j][t+1]*b[j][(int) o[t+1]];
		        	}
		        	beta[i][t] *= c2[t];
		    	}
	    }
	    return beta;
	  }
// Baum-Welch Algorithm
	  public void train(double[] o, int steps) {
	    int T = o.length;
	    double[][] alpha;
	    double[][] beta;
	    double pi1[] = new double[N];
	    double a1[][] = new double[N][N];
	    double b1[][] = new double[N][M];

	    for (int s = 0; s < steps; s++) {
	      /* calculation of alpha and beta from current model */
	      alpha = alphaPass(o);
	      beta = betaPass(o);
	      /* re-estimation of pi */
	      for (int i = 0; i < N; i++)
	    	  pi1[i] = gamma(i, 0, o, alpha, beta);
	      /* re-estimation of A */ 
	      for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
			  double num = 0;
			  double denom = 0;
			  for (int t = 0; t <= T - 1; t++) {
			    num += p(t, i, j, o, alpha, beta);
			    denom += gamma(o, alpha, beta);
			  }
			  a1[i][j] =  divide(num,denom);
			}
	      }
	      /* re-estimation of B  */
	      for (int i = 0; i < N; i++) {
			for (int k = 0; k < M; k++) {
			  double num = 0;
			  double denom = 0;
			  for (int t = 0; t <= T - 1; t++) {
			    double[][] g = gamma(o, alpha, beta);
			    num += g * (k == o[t] ? 1 : 0);
			    denom += g;
			  }
			  b1[i][k] =  divide(num,denom);;
			}
	      }
	      pi = pi1;
	      a = a1;
	      b = b1;
		  }
	  	}

	  public double[][] gamma(double[] o, double[][] alpha, double[][] beta) {
		  int T = o.length;
		  double denom;
		  double[][] gamma = new double[N][T];
		  double[][][] digamma = new double[N][N][T];
		  for( int t = 0 ; t<T-1; t++) {
			denom = 0.0;
			for( int i = 0; i<N ; i++) {
				for( int j = 0; j< N; j++) {
					denom += alpha[i][t] * a[i][j] * b[j][(int) o[t+1]] * beta[j][t+1];
				}
			}
			for( int i = 0; i<N; i++) {
				gamma[i][t] = 0;
				for ( int j = 0; j<N ; j++) {
					digamma[i][j][t] = (alpha[i][t] * a[i][j] * b[j][(int) o[t+1]]) / denom;
					gamma[i][t] += digamma[i][j][t];
				}
			}
		  }
		  denom = 0;
		  for(int i = 0; i<N ; i++) {
			  denom += alpha[i][T-1];
		  }
		  for(int i = 0; i<N ; i++) {
			  gamma[i][T-1] = alpha[i][T-1] / denom;
		  }
		  return gamma;
	  }
	  
	  public void print() {
		NumberFormat formatter = new DecimalFormat("#.#####"); 
		System.out.println();
	    for (int i = 0; i < N; i++)
	    	System.out.println("pi[" + i + "]: " + formatter.format(pi[i]));
	    	System.out.println();

	    for (int i = 0; i < N; i++) {
		    	for (int j = 0; j < N; j++) {
		    		System.out.print("a[" + i + "," + j + "]: " + formatter.format(a[i][j]) + "  ");
		    	}
	    		System.out.println();
	    }

	    System.out.println();
	    for (int i = 0; i < N; i++) {
		    	for (int k = 0; k < M; k++)
		    		System.out.print("b[" + i + "," + k + "]: " + formatter.format(b[i][k]) + "  ");
		    		System.out.println();
		    }
	  	}
	  
	  public double divide(double n, double d) {
	    if (n == 0)
	      return 0;
	    else
	      return n / d;
	  }
	  
	public double evaluateUsingBruteForce(List<Double> states, List<Double> observations) throws Exception{
      if (states.size() != observations.size())
          throw new Exception("States and Observations must be at a same size!");
      double previousState = 0.0;
      double probability = 0.0;
      double result = 0.0;
 	  for (int i = 0; i<states.size(); i++){
 		  previousState=0.0;
 		  probability = pi[states.get(i).intValue()];
 		  for(int j = 0; j<observations.size(); j++){
 			  double emissionValue = b[states.get(j).intValue()][observations.get(j).intValue()];
 			  double transitionValue = 0.0;
 			  if(j!=0){
 				  transitionValue += a[(int) previousState][states.get(j).intValue()];
 			      probability *= transitionValue * emissionValue;
 			  }
 			  previousState = states.get(j).intValue();
 		  }
 		  result+=probability;
 	  }
 	  return result;
  }  
	public static void printMatrix(double[][] matrix) {
		System.out.println("Matrix size: " +matrix.length + " x " + matrix[0].length);
	    for (double[] row : matrix) 
	        System.out.println(Arrays.toString(row));       
	}
	public static void printArray(double[] array) {
		System.out.println("Array Length: "+array.length);
	    for (double element : array) {
	        System.out.print(element + " ");
	    }
	}
	@SuppressWarnings("unchecked")
	public static void main(String[] args) throws Exception {
		NumberFormat nf = new DecimalFormat("0.######");
		HMM hmm = new HMM(2,3);
		hmm.a = new double[][]{{0.7,0.3},{0.4,0.6}};
		hmm.b = new double[][]{{0.1,0.4,0.5},{0.7,0.2,0.1}};
		hmm.pi = new double[]{0.0, 1.0};
		double[] states = new double[]{0.0,0.0,0.0,0.0};
		double[] observations = new double[]{1.0,0.0,0.0,2.0};
	    ICombinatoricsVector<Double> observationVector = Factory.createVector(new Double[] {0.0,1.0,2.0});
	    ICombinatoricsVector<Double> stateVector = Factory.createVector(new Double[]{0.0,1.0});
	    Generator<Double> observationGen = Factory.createPermutationWithRepetitionGenerator(observationVector, 4);
	    Generator<Double> stateGen = Factory.createPermutationWithRepetitionGenerator(stateVector, 4);
	    // Brute Force
	    for (ICombinatoricsVector<Double> ob : observationGen){
		    	for (ICombinatoricsVector<Double> st : stateGen){
		    		double result = hmm.evaluateUsingBruteForce(st.getVector(), ob.getVector());
		    		System.out.println("Brute Force results from using observation sequence " + ob.getVector() +
		    							" and state sequence " + st.getVector() +": " + result);
		    		System.out.println("");
		    	}
	    }
	    //Alpha Pass
	    for(ICombinatoricsVector<Double> ob: observationGen){
		    	double[] observationSequence = ob.getVector().stream().mapToDouble(d -> d).toArray(); 
	    		double[][] alphaPass = hmm.alphaPass(observationSequence);
	    		System.out.println("Alpha Pass Results from observationSequence "+ ob.getVector() +
	    							" with results below:");
	    		printMatrix(alphaPass);
	    		System.out.println("");
	    }
	    //Problem 10
	    HMM hmmtext = new HMM(2, 27);
		hmmtext.a = new double[][]{{0.45,0.55},{0.25,0.75}};
		hmmtext.b = new double[][]{
			{0.039,0.035,0.035,0.036,0.038,0.039,0.039,0.035,0.035,0.036,0.038,0.039,0.039,0.035,0.035,0.036,0.038,0.039,0.039,0.035,0.035,0.036,0.038,0.039,0.037,0.037,0.037},
			{0.039,0.035,0.035,0.036,0.038,0.039,0.039,0.035,0.035,0.036,0.038,0.039,0.039,0.035,0.035,0.036,0.038,0.039,0.039,0.035,0.035,0.036,0.038,0.039,0.037,0.037,0.037}};
		hmmtext.pi = new double[]{0.35,0.65};
		System.out.println("");
		double[] doubleArr = {};
	    try(BufferedReader br = new BufferedReader(new FileReader("resources/BrownCorpus"))) {
	        StringBuilder sb = new StringBuilder();
	        String line = br.readLine();
	        while (line != null) {
	        	line = line.substring(14);
	        	line = line.replaceAll("[^a-zA-Z ]", "").toLowerCase();
	            sb.append(line);
	            line = br.readLine();
	        }
	        String everything = sb.toString();
	        everything = everything.trim().replaceAll(" +", " ");
	        everything = everything.substring(0,50000);
	        System.out.println(everything);
	        String[] corpusArr = everything.split("");
		    for(String s:corpusArr){
		    		int num ="abcdefghijklmnopqrstuvwxyz ".indexOf(s);
		    		doubleArr = addElementDouble(doubleArr, (double) num);
		    }
		    System.out.println("double array: ");
		    printArray(doubleArr);
	    }
	    hmmtext.print();
	    hmmtext.train(doubleArr, 100);
//	    printMatrix(hmmtext.alphaPass(doubleArr));
//	    printMatrix(hmmtext.betaPass(doubleArr));
	    hmmtext.print();
//		printMatrix(hmmtext.a);
//		printMatrix(hmmtext.b);
//		printArray(hmmtext.pi);
	}

	static int[] addElement(int[] a, int e) {
	    a  = Arrays.copyOf(a, a.length + 1);
	    a[a.length - 1] = e;
	    return a;
	}
	static double[] addElementDouble(double[] a, double e) {
	    a  = Arrays.copyOf(a, a.length + 1);
	    a[a.length - 1] = e;
	    return a;
	}
}
