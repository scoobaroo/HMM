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

	  /** implementation of the Baum-Welch Algorithm for HMMs.
	      @param o the training set
	      @param steps the number of steps
	  */
	  public void train(double[] o, int steps) {
	    int T = o.length;
	    double[][] alpha;
	    double[][] beta;
	    double pi1[] = new double[N];
	    double a1[][] = new double[N][N];
	    double b1[][] = new double[N][M];

    	printMatrix(a);
		printMatrix(b);
	    for (int s = 0; s < steps; s++) {
	      /* calculation of Forward- and Backward Variables from the
		 current model */
	      alpha = alphaPass(o);
	      beta = betaPass(o);
	      /* re-estimation of initial state probabilities */
	      for (int i = 0; i < N; i++)
	    	  pi1[i] = gamma(i, 0, o, alpha, beta);
	      /* re-estimation of transition probabilities */ 
	      for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
			  double num = 0;
			  double denom = 0;
			  for (int t = 0; t <= T - 1; t++) {
			    num += p(t, i, j, o, alpha, beta);
			    denom += gamma(i, t, o, alpha, beta);
			  }
			  a1[i][j] = num / denom;
			}
	      }
	      /* re-estimation of emission probabilities */
	      for (int i = 0; i < N; i++) {
			for (int k = 0; k < M; k++) {
			  double num = 0;
			  double denom = 0;
			  for (int t = 0; t <= T - 1; t++) {
			    double g = gamma(i, t, o, alpha, beta);
			    num += g * (k == o[t] ? 1 : 0);
			    denom += g;
			  }
			  b1[i][k] = num / denom;
			}
	      }
	      pi = pi1;
	      a = a1;
	      b = b1;
		  }
	  	}
	  
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
		  System.out.println("alpha");
		  printMatrix(alpha);
		  System.out.println("b");
		  printMatrix(b);
		  System.out.println("a");
		  printMatrix(a);
		  
		  /* induction */
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
				  alpha[i][t]= c[t] * alpha[i][t];
			  }
		  }
		  return alpha;
	  }
	  
	  public double[][] betaPass(double[] o) {
	    int T = o.length;
	    double[][] beta = new double[N][T];  
	    double[] c = new double[T];
	    /* initialization (time 0) */
	    for (int i = 0; i < N; i++){
	    	beta[i][T-1] = c[T-1];
	    }
	    for(int t=T-2; t>=0; t--){
	    	for(int i = 0; i<N; i++){
	        	beta[i][t] = 0;
	        	for(int j =0; j<N; j++){
	        		beta[i][t] += a[i][j]*beta[j][t+1]*b[j][(int) o[t+1]];
	        	}
	        	beta[i][t] *= c[t];
	    	}
	    }
	    return beta;
	  }

	  public double p(int t, int i, int j, double[] o, double[][] alpha, double[][] beta) {
	    double num;
	    if (t == o.length - 1)
	      num = alpha[i][t] * a[i][j];
	    else
	      num = alpha[i][t] * a[i][j] * b[j][(int) o[t+1]] * beta[j][t+1];
	    double denom = 0;
	    for (int k = 0; k < N; k++)
	      denom += (alpha[k][t] * beta[k][t]);
	    return  num / denom;
	  }

	  public double gamma(int i, int t, double[] o, double[][] alpha, double[][] beta) {
	    double num = alpha[i][t] * beta[i][t];
	    double denom = 0;
	    for (int j = 0; j < N; j++)
	      denom += alpha[j][t] * beta[j][t];
	    return  num / denom;
	  }

	  public void print() {
		NumberFormat formatter = new DecimalFormat("#.#####"); 
	    
	    for (int i = 0; i < N; i++)
	    	System.out.println("pi(" + i + ") = " + formatter.format(pi[i]));
	    	System.out.println();

	    for (int i = 0; i < N; i++) {
	    	for (int j = 0; j < N; j++)
	    		System.out.print("a(" + i + "," + j + ") = " + 
	    							formatter.format(a[i][j]) + "  ");
	    	System.out.println();
	    }

	    System.out.println();
	    for (int i = 0; i < N; i++) {
	    	for (int k = 0; k < M; k++)
	    		System.out.print("b(" + i + "," + k + ") = " + 
	    							formatter.format(b[i][k]) + "  ");
	    	System.out.println();
	    }
	  }
	  
	  @SuppressWarnings("null")
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
	public static void printMatrix(double matrix[][]) {
	    for (double[] row : matrix) 
	        System.out.println(Arrays.toString(row));       
	}
	public static void printArray(double array[]) {
	    for (double element : array) 
	        System.out.print(element + " ");       
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
		    double total = 0.0;
	    	for (ICombinatoricsVector<Double> st : stateGen){
//	    		double[] stateSequence = st.getVector().stream().mapToDouble(d -> d).toArray();
	    		double result = hmm.evaluateUsingBruteForce(st.getVector(), ob.getVector());
	    		System.out.println("Brute Force results from using observation sequence " + ob.getVector() +
	    							" and state sequence " + st.getVector() +": " + result);
	    		System.out.println("");
	    		total+=result;
	    		System.out.println("total for this sequence:"+total);
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
	    String doi = "We hold these truths to be self-evident, that all men are created equal,"+ 
						"that they are endowed by their Creator with certain unalienable rights, that"+ 
						"mong these are life, liberty and the pursuit of happiness. That to secure "+
						"these rights, governments are instituted among men, deriving their just" +
						"powers from the consent of the governed. That whenever any form of "+
						"government becomes destructive of these ends, it is the right of the people"+ 
						"to alter or to abolish it, and to institute new government, laying its ";
	    String words = doi.replaceAll("[^a-zA-Z ]", "").toLowerCase();
	    System.out.println(words);
	    String[] charArr = words.split("");
	    double[] doubleArr = {};
	    for(String s:charArr){
	    	int num ="abcdefghijklmnopqrstuvwxyz ".indexOf(s);
	    	doubleArr = addElementDouble(doubleArr, (double) num);
	    	System.out.print(s);
	    }
	    HMM hmmtext = new HMM(2, 27);
		hmmtext.a = new double[][]{{0.7,0.3},{0.4,0.6}};
		hmmtext.b = new double[][]{{0.039,0.035,0.035,0.036,0.038,0.039,0.039,0.035,0.035,0.036,0.038,0.039,0.039,0.035,0.035,0.036,0.038,0.039,0.039,0.035,0.035,0.036,0.038,0.039,0.037,0.037,0.037},
			{0.039,0.035,0.035,0.036,0.038,0.039,0.039,0.035,0.035,0.036,0.038,0.039,0.039,0.035,0.035,0.036,0.038,0.039,0.039,0.035,0.035,0.036,0.038,0.039,0.037,0.037,0.037}};
		hmmtext.pi = new double[]{0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
		System.out.println("");
		printMatrix(hmmtext.a);
		printMatrix(hmmtext.b);
		System.out.println("training now");
	    hmmtext.train(doubleArr, 100);
	    hmmtext.print();
	    double[][] alpha = hmmtext.alphaPass(doubleArr);
	    printMatrix(alpha);
	    double[][] beta = hmmtext.betaPass(doubleArr);
	    printMatrix(beta);
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
	        everything = everything.substring(0,10000);
	        System.out.println(everything);
	        String[] corpusArr = everything.split("");
		    double[] doubleArr2 = {};
		    for(String s:corpusArr){
		    	int num ="abcdefghijklmnopqrstuvwxyz ".indexOf(s);
		    	doubleArr2 = addElementDouble(doubleArr2, (double) num);
		    }
		    System.out.println("double array 2");
		    printArray(doubleArr2);
	    }
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
