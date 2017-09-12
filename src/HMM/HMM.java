package HMM;
import java.io.BufferedReader;
import java.io.FileReader;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
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
	  double[] c = new double[50000];
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
	  public void trainWithoutA(double[] o, int steps) {
		    int T = o.length;
		    double[][] alpha;
		    double[][] beta;
		    double[][] gamma;
		    double[][][] digamma;
		    double delta = 0.0;
		    double epsilon = 0.0005;
		    double logProb = 0.0;
		    double oldLogProb = Double.NEGATIVE_INFINITY;
		    for ( int s = 0; s< steps; s++) {
			    logProb = computeLogProb(o);
		    		delta = Math.abs(logProb - oldLogProb);
		    		if (s<steps || delta> epsilon) {
				    oldLogProb = logProb;
		    			reestimateWithoutA(o);
		    		}
		    }
		  }
	  public void train(double[] o, int steps) {
	    int T = o.length;
	    double[][] alpha;
	    double[][] beta;
	    double[][] gamma;
	    double[][][] digamma;
	    double delta = 0.0;
	    double epsilon = 0.0005;
	    double logProb = 0.0;
	    double oldLogProb = Double.NEGATIVE_INFINITY;
	    for ( int s = 0; s< steps; s++) {
		    logProb = computeLogProb(o);
	    		delta = Math.abs(logProb - oldLogProb);
	    		if (s<steps || delta> epsilon) {
			    oldLogProb = logProb;
	    			reestimate(o);
	    		}
	    }
	  }
	  public void reestimateWithoutA(double[] o) {
		  	int T = o.length;
		    double[][] alpha;
		    double[][] beta;
		    double[][] gamma;
		    double[][][] digamma;
		    alpha = alphaPass(o);
		    beta = betaPass(o);
		    gamma = gamma(o,alpha,beta);
		    digamma = digamma(o,alpha,beta);
		    double numer;
		    double denom;
		    //reestimate pi
		    for ( int i =0; i< N ; i++) {
		    		pi[i] = gamma[i][0];
		    }
		    //reestimate b
		    for( int i = 0; i < N; i++) {
		    		for(int j = 0; j < M; j++) {
		    			numer = 0;
		    			denom = 0;
		    			for( int t = 0 ; t < T ; t++) {
		    				if( (int) o[t] == j ) {
		    					numer += gamma[i][t];
		    				}
		    				denom += gamma[i][t];
		    			}
		    			b[i][j] = divide(numer,denom);
		    		}
		    }
	  }
	  public void reestimate(double[] o) {
		  	int T = o.length;
		    double[][] alpha;
		    double[][] beta;
		    double[][] gamma;
		    double[][][] digamma;
		    alpha = alphaPass(o);
		    beta = betaPass(o);
		    gamma = gamma(o,alpha,beta);
		    digamma = digamma(o,alpha,beta);
		    double numer;
		    double denom;
		    //reestimate pi
		    for ( int i =0; i< N ; i++) {
		    		pi[i] = gamma[i][0];
		    }
		    //reestimate a
		    for ( int i = 0 ; i< N; i++) {
		    		for( int j=0; j< N; j++) {
		    			numer= 0;
		    			denom= 0;
		    			for( int t = 0; t< T-1; t++) {
						numer+= digamma[i][j][t];
						denom+= gamma[i][t];
		    			}
		    			a[i][j] = divide(numer,denom);
		    		}
		    }
		    //reestimate b
		    for( int i = 0; i < N; i++) {
		    		for(int j = 0; j < M; j++) {
		    			numer = 0;
		    			denom = 0;
		    			for( int t = 0 ; t < T ; t++) {
		    				if( (int) o[t] == j ) {
		    					numer += gamma[i][t];
		    				}
		    				denom += gamma[i][t];
		    			}
		    			b[i][j] = divide(numer,denom);
		    		}
		    }
	  }
	  public double computeLogProb(double[] o) {
		  int T = o.length;
		  double logProb = 0.0;
		  for (int i = 0; i<T-1; i++) {
			  logProb = logProb + Math.log(c[i]);
		  }
		  logProb = -logProb;
		  return logProb;
	  }
	  public double[][][] digamma(double[] o, double[][] alpha, double[][] beta){
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
		  return digamma;
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

	    	double[] atotal = new double[N];
	    	double[] btotal = new double[N];
	    	double pitotal = 0.0;
		for(int z = 0; z< pi.length; z++) {
	    		pitotal += pi[z];
	    }
		System.out.println("total for pi: "+ pitotal);
    		System.out.println();
	    for (int i = 0; i < N; i++) {
	    		atotal[i] = 0.0;
		    	for (int j = 0; j < N; j++) {
		    		System.out.print("a[" + i + "," + j + "]: " + formatter.format(a[i][j]) + "  ");
		    		atotal[i] += a[i][j];
		    	}
		    	System.out.println();
	    		System.out.println("total for a for row "+i+": "+atotal[i]);
	    }

	    System.out.println();
	    for (int i = 0; i < N; i++) {
	    		btotal[i] = 0.0;
		    	for (int k = 0; k < M; k++) {
		    		System.out.print("b[" + i + "," + k + "]: " + formatter.format(b[i][k]) + "  ");
		    		btotal[i] += b[i][k];
		    }
		    	System.out.println();
	    		System.out.println("total for b for row " + i + ": " + btotal[i]);
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
	public static void printMatrixInteger(Integer[][] matrix) {
		System.out.println("Matrix size: " +matrix.length + " x " + matrix[0].length);
	    for (Integer[] row : matrix) 
	        System.out.println(Arrays.toString(row));       
	}
	public static void printMatrixString(String[][] matrix) {
		System.out.println("Matrix size: " +matrix.length + " x " + matrix[0].length);
	    for (String[] row : matrix) 
	        System.out.println(Arrays.toString(row));       
	}
	public static void printArray(double[] array) {
		System.out.println("Array Length: "+array.length);
	    for (double element : array) {
	        System.out.print(element + " ");
	    }
	}
	public static void printArrayString(String[] array) {
		System.out.println("Array Length: "+array.length);
	    for (String s : array) {
	        System.out.print(s + " ");
	    }
	}
	static String cipher(String msg, int shift){
	    String s = "";
	    int len = msg.length();
	    for(int x = 0; x < len; x++){
	        char c = (char)(msg.charAt(x) + shift);
	        if (c > 'z')
	            s += (char)(msg.charAt(x) - (26-shift));
	        else
	            s += (char)(msg.charAt(x) + shift);
	    }
	    return s;
	}
	public static String[] split(String string) {
		String[] stringArray = new String[string.length() *2];
		String[] stringArray2 = new String[string.length() *2];
		stringArray = string.split("");
		for(int i = 0; i < string.length(); i++) {
			if(i < string.length() - 1) {
				String firstChar = stringArray[i];
				String secondChar = stringArray[i+1];
				String finalStr = firstChar + secondChar;
				stringArray2[i] = finalStr;
			}
		}
		return stringArray2;
	}
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
	    SeedGenerator sg  = new SeedGenerator();
		hmmtext.a = new double[][]{sg.generate(2), sg.generate(2)};
		hmmtext.b = new double[][]{sg.generate(27), sg.generate(27)};
		hmmtext.pi = sg.generate(2);
		System.out.println("");
		String everything;
		double[] doubleArr = {};
		String[] everythingArr = {};
	    try(BufferedReader br = new BufferedReader(new FileReader("resources/BrownCorpus"))) {
	        StringBuilder sb = new StringBuilder();
	        String line = br.readLine();
	        while (line != null) {
	        	line = line.substring(14);
	        	line = line.replaceAll("[^a-zA-Z ]", "").toLowerCase();
	            sb.append(line);
	            line = br.readLine();
	        }
	        everything = sb.toString();
	        everything = everything.trim().replaceAll(" +", " ");
	        // remove below line if you want spaces, keep it if you want to remove them
	        everything = everything.trim().replaceAll(" ", "");
	        
	        everything = everything.substring(0,50000);
	        System.out.println(everything);
	        //generate cipher text with alphabet shift of 3
	        everything = cipher(everything,3);
	        everything = everything.replace("#", " ");
	        everythingArr = split(everything);
	        String[] corpusArr = everything.split("");
		    for(String s:corpusArr){
		    		int num ="abcdefghijklmnopqrstuvwxyz ".indexOf(s);
		    		doubleArr = addElementDouble(doubleArr, (double) num);
		    }
		    System.out.println("double array: ");
		    printArray(doubleArr);
	    }
	    System.out.println();
	    System.out.println(everything);
	    for( double num : doubleArr) {
    			char character = "abcdefghijklmnopqrstuvwxyz ".charAt((int) num);
    			System.out.print(character);
	    }
//	    hmmtext.train(doubleArr, 100);
//	    hmmtext.print();
	    ArrayList<String> stringArr = new ArrayList<>(); 
	    String[][] digraph = new String[26][26];
	    for (int i = 0; i<26; i++) {
	    		for(int j = 0; j<26; j++) {
	    			char first = "abcdefghijklmnopqrstuvwxyz ".charAt(i);
	    			char second = "abcdefghijklmnopqrstuvwxyz ".charAt(j);
	    			digraph[i][j]= ""+first+second;
	    			stringArr.add("" +first+second);
	    		}
	    }
//	    printMatrixString(digraph);
	    for(String s : stringArr) {
	    		System.out.print(s+" ");
	    }
	    Map<String, Integer> hashmap = new HashMap<String, Integer>();
	    for(int i = 0; i< stringArr.size(); i++) {
	    		hashmap.put(stringArr.get(i),0);
	    }
//	    System.out.println();
//	    System.out.println(Collections.singletonList(hashmap));
	    Map<Pair<Character,Character>,Double> map = new LinkedHashMap<Pair<Character,Character>, Double>();
	    for (int i = 0; i<26; i++) {
	    		for(int j = 0; j<26; j++) {
	    			char first = "abcdefghijklmnopqrstuvwxyz ".charAt(i);
	    			char second = "abcdefghijklmnopqrstuvwxyz ".charAt(j);
	    			Pair<Character,Character> p = new Pair<Character, Character>(first,second);
	    			map.put(p, (double) 0);
	    		}
	    }
	    printArrayString(everythingArr);	   
	    System.out.println();
	    for(int i = 0; i<everythingArr.length; i++) {
	    		String elem = everythingArr[i];
	    		if(elem != null) {
	    			char[] charArr = elem.toCharArray();
		    		Pair<Character,Character> pair = new Pair<Character,Character>(charArr[0], charArr[1]);
		    		map.put(pair, map.get(pair)+1);
	    		}
	    }
	    //add 5 to each element in everythingArr
	    for(Entry<Pair<Character, Character>, Double> entry: map.entrySet()) {
	    		double value = entry.getValue();
	    		value +=5;
	    		map.put(entry.getKey(), value);
	    }
	    System.out.println(Collections.singletonList(map));
	    int rowSum = 0;
	    double[][] charArr = new double[26][26];
	    for( Entry<Pair<Character, Character>, Double> entry: map.entrySet()) {
	    		char first = entry.getKey().getFirst();
	    		char second = entry.getKey().getSecond();
	    		double firstInt = Double.valueOf("abcdefghijklmnopqrstuvwxyz".indexOf(String.valueOf(first)));
	    		double secondInt = Double.valueOf("abcdefghijklmnopqrstuvwxyz".indexOf(String.valueOf(second)));
	    		charArr[(int) firstInt][(int) secondInt]= entry.getValue();
	    }
	    double[][] newCharArr = new double[26][26];
	    for(int i = 0; i< charArr.length; i++) {
	    		int total = 0;
	    		for( int j = 0; j< charArr[0].length; j++) {
	    			total+=charArr[i][j];
	    		}
	    		for(int k = 0 ; k<charArr[0].length; k++) {
	    			double newValue =  charArr[i][k]/Double.valueOf(total);
	    			newCharArr[i][k]= newValue;
	    		}
	    }
	    printMatrix(newCharArr);
//	    for(double[] row : newCharArr) {
//	    		double total = 0.0;
//	    		for(double elem: row) {
//	    			total+=elem;
//	    		}
//	    		System.out.println("total for row :"+ total);
//	    }
	    HMM hmmDigraph = new HMM(26,26);
	    hmmDigraph.a = newCharArr;
	    hmmDigraph.b = new double[][] { sg.generate(26), sg.generate(26), sg.generate(26), sg.generate(26), sg.generate(26),
	    									sg.generate(26), sg.generate(26), sg.generate(26), sg.generate(26), sg.generate(26),
	    									sg.generate(26), sg.generate(26), sg.generate(26), sg.generate(26), sg.generate(26),
	    									sg.generate(26), sg.generate(26), sg.generate(26), sg.generate(26), sg.generate(26),
	    									sg.generate(26), sg.generate(26), sg.generate(26), sg.generate(26), sg.generate(26),
	    									sg.generate(26)};
	    hmmDigraph.pi = sg.generate(26);
	    hmmDigraph.trainWithoutA(doubleArr, 100);
	    hmmDigraph.print();
//	    HMM hmmtext2 = new HMM(3, 27);
//		hmmtext2.a = new double[][]{sg.generate(3),sg.generate(3),sg.generate(3)};
//		hmmtext2.b = new double[][]{sg.generate(27), sg.generate(27), sg.generate(27)};
//		hmmtext2.pi = sg.generate(3);
//		hmmtext2.train(doubleArr, 100);
//		hmmtext2.print();
//		HMM hmmtext3 = new HMM(4, 27);
//		hmmtext3.a = new double[][]{sg.generate(4),sg.generate(4),sg.generate(4),sg.generate(4)};
//		hmmtext3.b = new double[][]{sg.generate(27), sg.generate(27), sg.generate(27), sg.generate(27)};
//		hmmtext3.pi = sg.generate(4);
//		hmmtext3.train(doubleArr, 100);
//		hmmtext3.print();
//		HMM hmmtext4 = new HMM(26, 27);
//		hmmtext4.a = new double[][]{ sg.generate(26), sg.generate(26), sg.generate(26), sg.generate(26), sg.generate(26),
//		 							sg.generate(26), sg.generate(26), sg.generate(26), sg.generate(26), sg.generate(26),
//		 							 sg.generate(26), sg.generate(26), sg.generate(26), sg.generate(26), sg.generate(26),
//		 							 sg.generate(26), sg.generate(26), sg.generate(26), sg.generate(26), sg.generate(26),
//		 							 sg.generate(26), sg.generate(26), sg.generate(26), sg.generate(26), sg.generate(26),
//		 							 sg.generate(26)};
//		hmmtext4.b = new double[][]{ sg.generate(27), sg.generate(27), sg.generate(27), sg.generate(27), sg.generate(27),
//			 						sg.generate(27), sg.generate(27), sg.generate(27), sg.generate(27), sg.generate(27),
//			 						 sg.generate(27), sg.generate(27), sg.generate(27), sg.generate(27), sg.generate(27),
//			 						 sg.generate(27), sg.generate(27), sg.generate(27), sg.generate(27), sg.generate(27),
//			 						 sg.generate(27), sg.generate(27), sg.generate(27), sg.generate(27), sg.generate(27),
//			 						 sg.generate(27), sg.generate(27)};
//		hmmtext4.pi = sg.generate(26);
//		hmmtext4.train(doubleArr, 100);
//		hmmtext4.print();
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
