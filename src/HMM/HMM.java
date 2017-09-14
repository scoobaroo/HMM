package HMM;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
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
	public double[] c = new double[1000000];
	public double[] c2 = new double[1000000];
	  
	public HMM(int N, int M) {
	    this.N = N;
	    this.M = M;
	    pi = new double[N];
	    a = new double[N][N];
	    b = new double[N][M];
	  }
	
	  public double[][] alphaPass(double[] o){
		  printMatrix(b);
		  System.out.println(o.length);
		  System.out.println(N);
		  System.out.println(M);
		  System.out.println(pi.length);
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
		  for(int i=0; i< N; i++){
		  	alpha[i][0] *= c[0];
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
	  public void randomRestartWithoutA(int numberOfTimes, int iterations, int n, int m, double[] o) {
		  for(int k=0; k<numberOfTimes ; k++) {
			  HMM hmm = new HMM(n,m);
			  SeedGenerator sg = new SeedGenerator();
			  hmm.pi = sg.generate(n);
			  hmm.a = new double[n][n];
			  hmm.b = new double[n][m];
			  for(int i = 0; i<N; i++) {
				  hmm.a[i] = sg.generate(n);
				  hmm.b[i] = sg.generate(m);
			  }
			  hmm.trainWithoutA(o, iterations);
			  hmm.printB();
		  }
	  }
	  public void randomRestart(int numberOfTimes, int iterations, int N,int M,double[] o) {
		  for(int n=0; n<numberOfTimes ; n++) {
			  HMM hmm = new HMM(N,M);
			  SeedGenerator sg = new SeedGenerator();
			  hmm.pi = sg.generate(N);
			  hmm.a = new double[N][N];
			  hmm.b = new double[N][M];
			  for(int i = 0; i<N; i++) {
				  hmm.a[i] = sg.generate(N);
				  hmm.b[i] = sg.generate(M);
			  }
			  hmm.train(o, iterations);
			  hmm.printB();
		  }
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
					digamma[i][j][t] = (alpha[i][t] * a[i][j] * b[j][(int) o[t+1]] * beta[j][t+1]) / denom;
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

	  public void printB() {
		NumberFormat formatter = new DecimalFormat("#.#####"); 
		System.out.println();
	    	double[] btotal = new double[N];
	    System.out.println();
	    for (int i = 0; i < N; i++) {
	    		btotal[i] = 0.0;
		    	for (int k = 0; k < M; k++) {
		    		System.out.print("B[" + i + "," + k + "]: " + formatter.format(b[i][k]) + "  ");
		    		btotal[i] += b[i][k];
		    }
		    	System.out.println();
	  	}
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
		System.out.println("total for pi: "+ formatter.format(pitotal));
    		System.out.println();
	    for (int i = 0; i < N; i++) {
	    		atotal[i] = 0.0;
		    	for (int j = 0; j < N; j++) {
		    		System.out.print("A[" + i + "," + j + "]: " + formatter.format(a[i][j]) + "  ");
		    		atotal[i] += a[i][j];
		    	}
		    	System.out.println();
	    		System.out.println("total for A for row "+i+": "+formatter.format(atotal[i]));
	    }

	    System.out.println();
	    for (int i = 0; i < N; i++) {
	    		btotal[i] = 0.0;
		    	for (int k = 0; k < M; k++) {
		    		System.out.print("B[" + i + "," + k + "]: " + formatter.format(b[i][k]) + "  ");
		    		btotal[i] += b[i][k];
		    }
		    	System.out.println();
	    		System.out.println("total for B for row " + i + ": " + formatter.format(btotal[i]));
	  	}
	  }
	  
	  public double divide(double n, double d) {
	    if (n == 0)
	      return 0;
	    else
	      return n / d;
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
		System.out.println("String before shift of "+ shift);
		System.out.println(msg);
	    String s = "";
	    int len = msg.length();
	    for(int x = 0; x < len; x++){
	        char c = (char)(msg.charAt(x) + shift);
	        if (c > 'z')
	            s += (char)(msg.charAt(x) - (26-shift));
	        else
	            s += (char)(msg.charAt(x) + shift);
	    }
	    System.out.println("String after shift of "+ shift);
	    System.out.println(s);
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

	        //we are interested in 1,000,000 characters for the A matrix
//	        everything = everythingA.substring(0,1000000);
	        // we are interested only in the first 1000 or however many characters for training HMM N=26 M=26
//	        everything = everything.substring(0,1000);
	        //obtain 1000 plaintext characters for problem 11a
	        everything = everything.substring(3000000,3001000);
	        //generate cipher text with alphabet shift of 3
//	        everything = cipher(everything,3);
	        //generate cipher text with alphabet shift of 1
	        everything = cipher(everything,1);
	        everything = everything.replace("#", " ");
	        everythingArr = split(everything);
	        String[] corpusArr = everything.split("");
		    for(String s:corpusArr){
		    		int num ="abcdefghijklmnopqrstuvwxyz ".indexOf(s);
		    		doubleArr = addElementDouble(doubleArr, (double) num);
		    }
		    System.out.println("doubleArr: ");
		    printArray(doubleArr);
	    }
	    System.out.println();
	    for( double num : doubleArr) {
    			char character = "abcdefghijklmnopqrstuvwxyz ".charAt((int) num);
    			System.out.print(character);
	    }
	    System.out.println();
	    SeedGenerator sg = new SeedGenerator();
	    //training HMM with N=2, M=26
//	    HMM hmmCipher = new HMM(2,26);
//	    hmmCipher.pi = sg.generate(2);
//	    hmmCipher.a = new double[][] {sg.generate(2), sg.generate(2)};
//	    hmmCipher.b = new double[][] {sg.generate(26), sg.generate(26)};
//	    hmmCipher.print();
//	    hmmCipher.train(doubleArr, 100);
//	    hmmCipher.print();
//   
//	    ArrayList<String> stringArr = new ArrayList<>(); 
//	    String[][] digraph = new String[26][26];
//	    for (int i = 0; i<26; i++) {
//	    		for(int j = 0; j<26; j++) {
//	    			char first = "abcdefghijklmnopqrstuvwxyz ".charAt(i);
//	    			char second = "abcdefghijklmnopqrstuvwxyz ".charAt(j);
//	    			digraph[i][j]= ""+first+second;
//	    			stringArr.add("" +first+second);
//	    		}
//	    }
//	    printMatrixString(digraph);
//	    for(String s : stringArr) {
//	    		System.out.print(s+" ");
//	    }
//	    Map<String, Integer> hashmap = new HashMap<String, Integer>();
//	    for(int i = 0; i< stringArr.size(); i++) {
//	    		hashmap.put(stringArr.get(i),0);
//	    }
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
//	    printArrayString(everythingArr);	   
	    System.out.println();
	    //here we are going through everythingArr(elements of 2 characters), and incrementing map's pair(char1,char2) + 1 for each element
	    // element = pair of 2 characters
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
//	    System.out.println(Collections.singletonList(map));
	    int rowSum = 0;
	    double[][] charArr = new double[26][26];
	    //we are converting characters to integer values
	    for( Entry<Pair<Character, Character>, Double> entry: map.entrySet()) {
	    		char first = entry.getKey().getFirst();
	    		char second = entry.getKey().getSecond();
	    		double firstInt = Double.valueOf("abcdefghijklmnopqrstuvwxyz".indexOf(String.valueOf(first)));
	    		double secondInt = Double.valueOf("abcdefghijklmnopqrstuvwxyz".indexOf(String.valueOf(second)));
	    		charArr[(int) firstInt][(int) secondInt]= entry.getValue();
	    }
//	    making a new A matrix from the old one
	    double[][] newCharArr = new double[26][26];
	    for(int i = 0; i< charArr.length; i++) {
	    		//we are calculating row sums and storing it in total
	    		int total = 0;
	    		for( int j = 0; j< charArr[0].length; j++) {
	    			total+=charArr[i][j];
	    		}
	    		//we are dividing each element of a row by the row sum to obtain row stochasticity
	    		for(int k = 0 ; k<charArr[0].length; k++) {
	    			double newValue =  charArr[i][k]/Double.valueOf(total);
	    			newCharArr[i][k]= newValue;
	    		}
	    }
	    //show A matrix
	    printMatrix(newCharArr);
	    // print Totals in A Matrix rows
//	    for(double[] row : newCharArr) {
//	    		double total = 0.0;
//	    		for(double elem: row) {
//	    			total+=elem;
//	    		}
//	    		System.out.println("total for row :"+ total);
//	    }
	    //training HMM with N=M=26
	    HMM hmmDigraph = new HMM(26,26);
	    hmmDigraph.a = newCharArr;
//	    hmmDigraph.a = new double[][] 
//	    		{{0.0029371617998342515, 0.022717301223614293, 0.046543655243016624, 0.04077901818359089, 0.0012796763028323501, 0.012687076488080729, 0.023728854872519866, 0.0043508994296299905, 0.035282503778091945, 0.001620923316920977, 0.011102715351240676, 0.1041412762638327, 0.035940623019548576, 0.19393799054258276, 0.0025959147857456248, 0.020742943499244382, 6.21557061375713E-4, 0.11400087749232195, 0.099485692000195, 0.1451152927411885, 0.012504265587676108, 0.021449812314142253, 0.009969287768732023, 0.002656851752547165, 0.03215034368449276, 0.0016574854970019012},
//	    		{0.09310193596892426, 0.0070171041914667, 0.0026314140718000124, 0.001441012467890483, 0.3028632291209824, 0.0010650961719190527, 8.771380239333375E-4, 0.0016916233318714365, 0.06108639809535743, 0.006327924315519077, 5.638744439571456E-4, 0.11985464569889105, 0.002694066787795251, 7.518325919428607E-4, 0.11402794311133388, 0.0010024434559238143, 3.1326357997619195E-4, 0.059332122047490755, 0.018482551218595326, 0.008144853079380992, 0.10832654595576718, 0.003320593947747635, 0.0022554977758285823, 3.1326357997619195E-4, 0.08213771066975753, 3.7591629597143037E-4},
//	    		{0.1340026189436927, 0.0018706740662218619, 0.01923676498098148, 0.0025254099893995134, 0.15239758059487435, 0.0016212508573922803, 0.0010600486375257217, 0.16215626364033173, 0.06703248737295005, 4.053127143480701E-4, 0.036104009478081935, 0.03647814429132631, 0.001777140362910769, 0.0010912265386294195, 0.19879029743717652, 0.0030554343081623746, 0.0010600486375257217, 0.03703934651119287, 0.009509259836627797, 0.09347134750888571, 0.028247178399950116, 6.547359231776517E-4, 0.001777140362910769, 1.558895055184885E-4, 0.00826214379247989, 2.1824530772588389E-4},
//	    		{0.10246402209123254, 0.037836546120758326, 0.022303648239604908, 0.022250544315224895, 0.16855185598215708, 0.023471934575965164, 0.016621528330943658, 0.02718920928256598, 0.12529870957463757, 0.004912113005151081, 0.001858637353300409, 0.019913971642504382, 0.028623015240826298, 0.016674632255323667, 0.07840794434708726, 0.018612925495194096, 0.0016993255801603737, 0.03337581647283734, 0.061202272847963465, 0.1009240082842122, 0.03300408900217726, 0.00695661409378153, 0.030056821199086613, 1.327598109500292E-4, 0.017285327385693803, 3.717274706600818E-4},
//		    	{0.08105897406825507, 0.017533487746805097, 0.049725457166458806, 0.0827677804953109, 0.037018740180437, 0.024579277280163268, 0.014715171933461831, 0.015565525842660232, 0.03165746124815757, 0.002607751988208426, 0.0048024749348062005, 0.04306840084873419, 0.038662757738220574, 0.1016780317141515, 0.02710604318178137, 0.02897682178201785, 0.0031827532029997246, 0.14204959587942792, 0.11036783880529324, 0.06223780754466383, 0.0074669171836278525, 0.019558139911563195, 0.02925217447642495, 0.011467629861189847, 0.01232608237904728, 5.669026061322664E-4},
//		    	{0.10180327146069397, 0.010140646355980778, 0.022618050350513647, 0.01053745425686698, 0.08258013315109564, 0.05916846699880958, 0.009920197522155108, 0.017724086239583794, 0.10832855694193377, 0.0036153608747409726, 0.00185177020413562, 0.026542039592610554, 0.015916405802213308, 0.008685684052731362, 0.18394250694413827, 0.014240994665138222, 6.613465014770072E-4, 0.08002292667871787, 0.022353511749922842, 0.15995767382390547, 0.037079493849477535, 0.0023808474053172257, 0.013359199329835545, 2.2044883382566906E-4, 0.005775759446232529, 5.731669679467395E-4},
//		    	{0.10664407475867735, 0.012733620866707743, 0.013298418566440748, 0.008420620250564798, 0.16050523721503387, 0.014684740193058123, 0.01494146642020949, 0.11414048059149723, 0.0918052988293284, 0.0027726432532347504, 0.0012836311357568289, 0.030499075785582256, 0.015095502156500308, 0.03008831382214007, 0.098120764017252, 0.012939001848428836, 6.674881905935511E-4, 0.09581022797288971, 0.03835489833641405, 0.07943109468063257, 0.03224481413021154, 0.0024132265352228384, 0.01591702608338468, 2.567262271513658E-4, 0.0066235366605052375, 3.0807147258163895E-4},
//		    	{0.16200396634383965, 0.0045439667289215785, 0.007759400812522865, 0.00294587673527543, 0.4656410651366078, 0.0036390242023990604, 0.002040934208752912, 0.005968769855786819, 0.12619134720911873, 8.279261412865587E-4, 6.353851781966613E-4, 0.003735294683944009, 0.0065849009376744904, 0.007451335271579028, 0.09315131794289234, 0.003927835647033907, 2.1179505939888712E-4, 0.01971619462040549, 0.010628261162562335, 0.04126152839016501, 0.015480293432427749, 0.001097483489612415, 0.00812522864239367, 9.627048154494869E-5, 0.006161310818876716, 1.7328686678090764E-4},
//		    	{0.031066686721452696, 0.01039202548780988, 0.07716078924698837, 0.03724720714315015, 0.04119891157206733, 0.018008286272954754, 0.028413985478511754, 0.001449414081194536, 9.981813955396333E-4, 3.418429436779566E-4, 0.005127644155169349, 0.05433935432704798, 0.030218916221131365, 0.24771307070679446, 0.07694200976303447, 0.01015957228610887, 8.887916535626871E-4, 0.03619433087662204, 0.1279449769597856, 0.12213364691726034, 0.001449414081194536, 0.029986463019430355, 0.0015724775409186004, 0.002023710226573503, 1.0938974197694612E-4, 0.0069189011800418415},
//		    	{0.1486810551558753, 0.0057553956834532375, 0.005275779376498801, 0.004316546762589928, 0.15587529976019185, 0.0028776978417266188, 0.0038369304556354917, 0.0057553956834532375, 0.028297362110311752, 0.0038369304556354917, 0.002398081534772182, 0.003357314148681055, 0.004316546762589928, 0.0028776978417266188, 0.294484412470024, 0.003357314148681055, 0.002398081534772182, 0.026378896882494004, 0.004796163069544364, 0.003357314148681055, 0.2729016786570743, 0.002398081534772182, 0.005275779376498801, 0.002398081534772182, 0.002398081534772182, 0.002398081534772182},
//		    	{0.06977818853974121, 0.0118607516943931, 0.016019716574245224, 0.009242144177449169, 0.3285582255083179, 0.015403573629081947, 0.006469500924214418, 0.03126925446703635, 0.14525569932224275, 0.005083179297597043, 0.002618607516943931, 0.024953789279112754, 0.012322858903265557, 0.05730129390018484, 0.05560690080098583, 0.010166358595194085, 0.0012322858903265558, 0.010320394331484904, 0.08872458410351201, 0.04882932840418977, 0.006315465187923599, 0.0018484288354898336, 0.024029574861367836, 7.701786814540973E-4, 0.015249537892791128, 7.701786814540973E-4},
//		    	{0.11746211852587286, 0.015034394723778455, 0.013568777628064204, 0.06306881308654233, 0.1658274826844432, 0.017634683119400515, 0.00520057679124412, 0.0076117533035482115, 0.1276032432687989, 0.001040115358248824, 0.005389688674562088, 0.13597144410561898, 0.014372503132165567, 0.005011464907926152, 0.07807956882490603, 0.014041557336359124, 6.61891591612888E-4, 0.00898281445760348, 0.04089544476751058, 0.03680590029075952, 0.02352079048767227, 0.0075171973618892276, 0.01016476372834078, 1.1819492707373E-4, 0.084131149091081, 2.8366782497695197E-4},
//		    	{0.18915635203887918, 0.04073961576429493, 0.007517655099096363, 0.003037436403675298, 0.24098261067658897, 0.005695193256891184, 0.001822461842205179, 0.006606424177993773, 0.11276482648644544, 0.001328878426607943, 7.593591009188245E-4, 0.003872731414686005, 0.04165084668539752, 0.004290378920191358, 0.12070012909104716, 0.06348242083681373, 3.037436403675298E-4, 0.027678639228491154, 0.033183992710152634, 0.02729895967803174, 0.04267598147163794, 7.973270559647657E-4, 0.008542789885336776, 1.8983977522970614E-4, 0.014655630647733313, 2.657756853215886E-4},
//		    	{0.0740925740925741, 0.011474636474636474, 0.053058053058053056, 0.1426906426906427, 0.0892995892995893, 0.01567876567876568, 0.1108058608058608, 0.013875013875013876, 0.060786435786435784, 0.003357753357753358, 0.007506382506382506, 0.012737262737262738, 0.013167388167388168, 0.018689643689643688, 0.06746031746031746, 0.010142635142635142, 0.0011932511932511933, 0.007006882006882007, 0.06998556998556998, 0.16697191697191696, 0.012057387057387058, 0.006438006438006438, 0.016330891330891332, 2.4975024975024975E-4, 0.014221889221889222, 7.215007215007215E-4},
//		    	{0.01898014750588173, 0.018477887335113272, 0.02155753522429882, 0.02400274921356631, 0.008194771207274841, 0.11340505961035184, 0.012411113167146898, 0.009820508075814851, 0.013640328848238124, 0.0014935631393904148, 0.008617727140553543, 0.04451611197758334, 0.06455364931666183, 0.18047000978085595, 0.028998916175420972, 0.028959264056676096, 5.154775436834175E-4, 0.13779111263845198, 0.04142324671548283, 0.06069417642549367, 0.09178143752147823, 0.022231621242961748, 0.04022046578022152, 0.0011895635623463481, 0.0055248618784530384, 5.286949165983769E-4},
//		    	{0.1317783971443939, 0.0032761234169478267, 0.0023959708571707986, 0.0020536893061463986, 0.17304777272505012, 0.004498557527749254, 0.0012713314752334848, 0.032370055254021805, 0.06004596352256614, 7.334604664808567E-4, 9.290499242090852E-4, 0.09671898684660897, 0.009437191335387023, 8.80152559777028E-4, 0.13515231529020585, 0.05202679575570877, 2.9338418659234265E-4, 0.17783971443939173, 0.0230306586474989, 0.04317637279350643, 0.040046941469854776, 3.911789154564569E-4, 0.004498557527749254, 2.444868221602856E-4, 0.0036184049679722262, 2.444868221602856E-4},
//		    	{0.0044964028776978415, 0.0044964028776978415, 0.00539568345323741, 0.0044964028776978415, 0.00539568345323741, 0.0044964028776978415, 0.0044964028776978415, 0.00539568345323741, 0.0044964028776978415, 0.0044964028776978415, 0.0044964028776978415, 0.0044964028776978415, 0.0044964028776978415, 0.0044964028776978415, 0.0044964028776978415, 0.00539568345323741, 0.0044964028776978415, 0.0044964028776978415, 0.0044964028776978415, 0.00539568345323741, 0.8812949640287769, 0.00539568345323741, 0.0044964028776978415, 0.00539568345323741, 0.0044964028776978415, 0.0044964028776978415},
//		    	{0.10308531089470864, 0.010258107213765718, 0.024865273707100313, 0.033390060193501624, 0.22561532885821436, 0.011896883174182976, 0.015694431313226814, 0.012101730169235133, 0.10138350508965996, 0.0018751378777851313, 0.015127162711543915, 0.01974409883079638, 0.02945069490403706, 0.02650404966751757, 0.10138350508965996, 0.013992625508178123, 6.460559074721881E-4, 0.01878289370016703, 0.07889760801739623, 0.07989032807034131, 0.02018530774321641, 0.009706596073240679, 0.013094450222180202, 1.733320727364407E-4, 0.031861586461189376, 3.9393652894645614E-4},
//		    	{0.09845728964052919, 0.020638762268185125, 0.035932233589214393, 0.013185051002984454, 0.10870243062257791, 0.019287591500987394, 0.007023118383346449, 0.05388350235341282, 0.09723975114701035, 0.00267264547357793, 0.007438863234791905, 0.01717917118294258, 0.02348958410666825, 0.015219231168985435, 0.083921067870347, 0.038560334971566024, 0.00203418016600098, 0.010438165377362692, 0.07069147277613624, 0.19777576504476682, 0.03536800843368127, 0.002687493503986696, 0.03082451112859879, 2.524165169490267E-4, 0.006874638079258786, 2.2272045613149416E-4},
//		    	{0.06887022487130859, 0.010295312923327011, 0.011454890273638581, 0.0055377946356001085, 0.10546735302086155, 0.008615551341099974, 0.003500406393931184, 0.31280411812516934, 0.11806014630181523, 0.0015713898672446492, 0.0015605526957464101, 0.014966133839068004, 0.011238146843673802, 0.0056136548360877815, 0.11214305066377675, 0.007867786507721485, 3.793010024383636E-4, 0.038081820644811706, 0.043034408019506906, 0.05342725548631807, 0.0214467623950149, 0.001332972094283392, 0.020872392305608237, 6.502302898943376E-5, 0.021230018965050123, 5.635329179084259E-4},
//		    	{0.035906167234203555, 0.03185773741959894, 0.04903518728717367, 0.028717366628830874, 0.04309496783957624, 0.0061293984108967085, 0.03658721150208097, 0.0014755959137343927, 0.026598562239878925, 9.837306091562618E-4, 0.0018917896329928112, 0.0968974650018918, 0.03658721150208097, 0.1388951948543322, 0.0035944003026863415, 0.04074914869466515, 2.2701475595913735E-4, 0.1453651153991676, 0.13984108967082862, 0.12750662126371548, 4.161937192584185E-4, 0.00170261066969353, 0.0016269390843738176, 0.0012864169504351116, 0.002497162315550511, 5.297010972379871E-4},
//		    	{0.08454198473282443, 0.00200381679389313, 0.001049618320610687, 0.0018129770992366412, 0.6177480916030534, 0.0011450381679389313, 7.633587786259542E-4, 0.0014312977099236641, 0.21125954198473282, 8.587786259541985E-4, 0.0011450381679389313, 9.541984732824427E-4, 0.002480916030534351, 6.679389312977099E-4, 0.049904580152671754, 0.0012404580152671756, 5.725190839694657E-4, 0.0016221374045801526, 0.005057251908396947, 0.0022900763358778627, 0.002099236641221374, 9.541984732824427E-4, 0.0026717557251908397, 4.7709923664122136E-4, 0.004770992366412214, 4.7709923664122136E-4},
//		    	{0.18162500698987866, 0.00682212156796958, 0.0051445506906000115, 0.005759660012302187, 0.1719510149303808, 0.004361684281160879, 0.0018453279651065259, 0.17804618911815692, 0.1863781244757591, 0.001677570877369569, 0.0018453279651065259, 0.006933959626460885, 0.009058882737795673, 0.035732259687971814, 0.12044958899513504, 0.006262931275513057, 6.151093217021753E-4, 0.01627243751048482, 0.022535368785997874, 0.017390818095397866, 0.0017334899066152212, 0.0013979757311413073, 0.00631885030475871, 2.795951462282615E-4, 0.009170720796286976, 3.9143320471956607E-4},
//		    	{0.11284807034684904, 0.017586712261846604, 0.0986809965803615, 0.006350757205666829, 0.07278944797264289, 0.01172447484123107, 0.0039081582804103565, 0.023937469467513434, 0.09330727894479726, 0.0029311187103077674, 0.004885197850512946, 0.010258915486077186, 0.01074743527112848, 0.0039081582804103565, 0.032730825598436736, 0.2598925256472887, 0.003419638495359062, 0.008793356130923302, 0.01074743527112848, 0.1621885686370298, 0.008304836345872008, 0.005862237420615535, 0.01709819247679531, 0.003419638495359062, 0.011235955056179775, 0.002442598925256473},
//		    	{0.10702360391479562, 0.04277489925158319, 0.04214162348877375, 0.028151986183074264, 0.08370754173862982, 0.033045480713874496, 0.014219919401266552, 0.038860103626943004, 0.06810592976396085, 0.006390328151986183, 0.0037996545768566492, 0.023949337938975246, 0.04179620034542314, 0.020782959124928037, 0.10391479562464019, 0.03172135866436385, 0.0016119746689694876, 0.02498560736902706, 0.09884858952216465, 0.1158894645941278, 0.008175014392630972, 0.004778353483016695, 0.05043177892918826, 3.454231433506045E-4, 0.0034542314335060447, 0.001093839953943581},
//		    	{0.17531305903398928, 0.009838998211091235, 0.011627906976744186, 0.008050089445438283, 0.40518783542039355, 0.008944543828264758, 0.008050089445438283, 0.017889087656529516, 0.08586762075134168, 0.005366726296958855, 0.011627906976744186, 0.025044722719141325, 0.012522361359570662, 0.004472271914132379, 0.03667262969588551, 0.007155635062611807, 0.004472271914132379, 0.008944543828264758, 0.014311270125223614, 0.016100178890876567, 0.01699463327370304, 0.006261180679785331, 0.009838998211091235, 0.004472271914132379, 0.01520572450805009, 0.06976744186046512}};
	    hmmDigraph.b = new double[][] { sg.generate(26), sg.generate(26), sg.generate(26), sg.generate(26), sg.generate(26),
	    									sg.generate(26), sg.generate(26), sg.generate(26), sg.generate(26), sg.generate(26),
	    									sg.generate(26), sg.generate(26), sg.generate(26), sg.generate(26), sg.generate(26),
	    									sg.generate(26), sg.generate(26), sg.generate(26), sg.generate(26), sg.generate(26),
	    									sg.generate(26), sg.generate(26), sg.generate(26), sg.generate(26), sg.generate(26),
	    									sg.generate(26)};
	    hmmDigraph.pi = sg.generate(26);
	    // random restart 10 times
//	    hmmDigraph.randomRestartWithoutA(10, 100, 26, 26, doubleArr);
	    
//	    System.out.println("initial values for hmmDigraph");
//	    hmmDigraph.print();
//	    hmmDigraph.trainWithoutA(doubleArr, 100);
//	    System.out.println("values for hmmDigraph after training");
//	    hmmDigraph.print();
	    
	    // Zodiac408 Cipher
	    String Zstring;
		String[] ZArr;
		double[] ZdoubleArr = new double[408];
		try(BufferedReader br = new BufferedReader(new FileReader("resources/Z408.txt"))) {
	        StringBuilder sb = new StringBuilder();
	        String line = br.readLine();
	        while (line != null) {
	            sb.append(line + " ");
	            line = br.readLine();
	        }
	        Zstring = sb.toString();
	        Zstring = Zstring.replaceAll(" +", " ");
	        ZArr = Zstring.split(" ");
	        // remove below line if you want spaces, keep it if you want to remove them
	        System.out.println("Zstring: ");
	        System.out.println(Zstring);
		    System.out.println("ZArr: ");
		    printArrayString(ZArr);
		    for(int i=0; i<ZArr.length; i++) {
		    		double d = Double.valueOf(ZArr[i]);
		    		ZdoubleArr[i]= d;
		    }
		    printArray(ZdoubleArr);
	    }
		HMM hmmZodiac = new HMM(26,53);
		hmmZodiac.pi = sg.generate(26);
	    hmmZodiac.a = new double[][] 
			{{0.0029371617998342515, 0.022717301223614293, 0.046543655243016624, 0.04077901818359089, 0.0012796763028323501, 0.012687076488080729, 0.023728854872519866, 0.0043508994296299905, 0.035282503778091945, 0.001620923316920977, 0.011102715351240676, 0.1041412762638327, 0.035940623019548576, 0.19393799054258276, 0.0025959147857456248, 0.020742943499244382, 6.21557061375713E-4, 0.11400087749232195, 0.099485692000195, 0.1451152927411885, 0.012504265587676108, 0.021449812314142253, 0.009969287768732023, 0.002656851752547165, 0.03215034368449276, 0.0016574854970019012},
			{0.09310193596892426, 0.0070171041914667, 0.0026314140718000124, 0.001441012467890483, 0.3028632291209824, 0.0010650961719190527, 8.771380239333375E-4, 0.0016916233318714365, 0.06108639809535743, 0.006327924315519077, 5.638744439571456E-4, 0.11985464569889105, 0.002694066787795251, 7.518325919428607E-4, 0.11402794311133388, 0.0010024434559238143, 3.1326357997619195E-4, 0.059332122047490755, 0.018482551218595326, 0.008144853079380992, 0.10832654595576718, 0.003320593947747635, 0.0022554977758285823, 3.1326357997619195E-4, 0.08213771066975753, 3.7591629597143037E-4},
			{0.1340026189436927, 0.0018706740662218619, 0.01923676498098148, 0.0025254099893995134, 0.15239758059487435, 0.0016212508573922803, 0.0010600486375257217, 0.16215626364033173, 0.06703248737295005, 4.053127143480701E-4, 0.036104009478081935, 0.03647814429132631, 0.001777140362910769, 0.0010912265386294195, 0.19879029743717652, 0.0030554343081623746, 0.0010600486375257217, 0.03703934651119287, 0.009509259836627797, 0.09347134750888571, 0.028247178399950116, 6.547359231776517E-4, 0.001777140362910769, 1.558895055184885E-4, 0.00826214379247989, 2.1824530772588389E-4},
			{0.10246402209123254, 0.037836546120758326, 0.022303648239604908, 0.022250544315224895, 0.16855185598215708, 0.023471934575965164, 0.016621528330943658, 0.02718920928256598, 0.12529870957463757, 0.004912113005151081, 0.001858637353300409, 0.019913971642504382, 0.028623015240826298, 0.016674632255323667, 0.07840794434708726, 0.018612925495194096, 0.0016993255801603737, 0.03337581647283734, 0.061202272847963465, 0.1009240082842122, 0.03300408900217726, 0.00695661409378153, 0.030056821199086613, 1.327598109500292E-4, 0.017285327385693803, 3.717274706600818E-4},
		    	{0.08105897406825507, 0.017533487746805097, 0.049725457166458806, 0.0827677804953109, 0.037018740180437, 0.024579277280163268, 0.014715171933461831, 0.015565525842660232, 0.03165746124815757, 0.002607751988208426, 0.0048024749348062005, 0.04306840084873419, 0.038662757738220574, 0.1016780317141515, 0.02710604318178137, 0.02897682178201785, 0.0031827532029997246, 0.14204959587942792, 0.11036783880529324, 0.06223780754466383, 0.0074669171836278525, 0.019558139911563195, 0.02925217447642495, 0.011467629861189847, 0.01232608237904728, 5.669026061322664E-4},
		    	{0.10180327146069397, 0.010140646355980778, 0.022618050350513647, 0.01053745425686698, 0.08258013315109564, 0.05916846699880958, 0.009920197522155108, 0.017724086239583794, 0.10832855694193377, 0.0036153608747409726, 0.00185177020413562, 0.026542039592610554, 0.015916405802213308, 0.008685684052731362, 0.18394250694413827, 0.014240994665138222, 6.613465014770072E-4, 0.08002292667871787, 0.022353511749922842, 0.15995767382390547, 0.037079493849477535, 0.0023808474053172257, 0.013359199329835545, 2.2044883382566906E-4, 0.005775759446232529, 5.731669679467395E-4},
		    	{0.10664407475867735, 0.012733620866707743, 0.013298418566440748, 0.008420620250564798, 0.16050523721503387, 0.014684740193058123, 0.01494146642020949, 0.11414048059149723, 0.0918052988293284, 0.0027726432532347504, 0.0012836311357568289, 0.030499075785582256, 0.015095502156500308, 0.03008831382214007, 0.098120764017252, 0.012939001848428836, 6.674881905935511E-4, 0.09581022797288971, 0.03835489833641405, 0.07943109468063257, 0.03224481413021154, 0.0024132265352228384, 0.01591702608338468, 2.567262271513658E-4, 0.0066235366605052375, 3.0807147258163895E-4},
		    	{0.16200396634383965, 0.0045439667289215785, 0.007759400812522865, 0.00294587673527543, 0.4656410651366078, 0.0036390242023990604, 0.002040934208752912, 0.005968769855786819, 0.12619134720911873, 8.279261412865587E-4, 6.353851781966613E-4, 0.003735294683944009, 0.0065849009376744904, 0.007451335271579028, 0.09315131794289234, 0.003927835647033907, 2.1179505939888712E-4, 0.01971619462040549, 0.010628261162562335, 0.04126152839016501, 0.015480293432427749, 0.001097483489612415, 0.00812522864239367, 9.627048154494869E-5, 0.006161310818876716, 1.7328686678090764E-4},
		    	{0.031066686721452696, 0.01039202548780988, 0.07716078924698837, 0.03724720714315015, 0.04119891157206733, 0.018008286272954754, 0.028413985478511754, 0.001449414081194536, 9.981813955396333E-4, 3.418429436779566E-4, 0.005127644155169349, 0.05433935432704798, 0.030218916221131365, 0.24771307070679446, 0.07694200976303447, 0.01015957228610887, 8.887916535626871E-4, 0.03619433087662204, 0.1279449769597856, 0.12213364691726034, 0.001449414081194536, 0.029986463019430355, 0.0015724775409186004, 0.002023710226573503, 1.0938974197694612E-4, 0.0069189011800418415},
		    	{0.1486810551558753, 0.0057553956834532375, 0.005275779376498801, 0.004316546762589928, 0.15587529976019185, 0.0028776978417266188, 0.0038369304556354917, 0.0057553956834532375, 0.028297362110311752, 0.0038369304556354917, 0.002398081534772182, 0.003357314148681055, 0.004316546762589928, 0.0028776978417266188, 0.294484412470024, 0.003357314148681055, 0.002398081534772182, 0.026378896882494004, 0.004796163069544364, 0.003357314148681055, 0.2729016786570743, 0.002398081534772182, 0.005275779376498801, 0.002398081534772182, 0.002398081534772182, 0.002398081534772182},
		    	{0.06977818853974121, 0.0118607516943931, 0.016019716574245224, 0.009242144177449169, 0.3285582255083179, 0.015403573629081947, 0.006469500924214418, 0.03126925446703635, 0.14525569932224275, 0.005083179297597043, 0.002618607516943931, 0.024953789279112754, 0.012322858903265557, 0.05730129390018484, 0.05560690080098583, 0.010166358595194085, 0.0012322858903265558, 0.010320394331484904, 0.08872458410351201, 0.04882932840418977, 0.006315465187923599, 0.0018484288354898336, 0.024029574861367836, 7.701786814540973E-4, 0.015249537892791128, 7.701786814540973E-4},
		    	{0.11746211852587286, 0.015034394723778455, 0.013568777628064204, 0.06306881308654233, 0.1658274826844432, 0.017634683119400515, 0.00520057679124412, 0.0076117533035482115, 0.1276032432687989, 0.001040115358248824, 0.005389688674562088, 0.13597144410561898, 0.014372503132165567, 0.005011464907926152, 0.07807956882490603, 0.014041557336359124, 6.61891591612888E-4, 0.00898281445760348, 0.04089544476751058, 0.03680590029075952, 0.02352079048767227, 0.0075171973618892276, 0.01016476372834078, 1.1819492707373E-4, 0.084131149091081, 2.8366782497695197E-4},
		    	{0.18915635203887918, 0.04073961576429493, 0.007517655099096363, 0.003037436403675298, 0.24098261067658897, 0.005695193256891184, 0.001822461842205179, 0.006606424177993773, 0.11276482648644544, 0.001328878426607943, 7.593591009188245E-4, 0.003872731414686005, 0.04165084668539752, 0.004290378920191358, 0.12070012909104716, 0.06348242083681373, 3.037436403675298E-4, 0.027678639228491154, 0.033183992710152634, 0.02729895967803174, 0.04267598147163794, 7.973270559647657E-4, 0.008542789885336776, 1.8983977522970614E-4, 0.014655630647733313, 2.657756853215886E-4},
		    	{0.0740925740925741, 0.011474636474636474, 0.053058053058053056, 0.1426906426906427, 0.0892995892995893, 0.01567876567876568, 0.1108058608058608, 0.013875013875013876, 0.060786435786435784, 0.003357753357753358, 0.007506382506382506, 0.012737262737262738, 0.013167388167388168, 0.018689643689643688, 0.06746031746031746, 0.010142635142635142, 0.0011932511932511933, 0.007006882006882007, 0.06998556998556998, 0.16697191697191696, 0.012057387057387058, 0.006438006438006438, 0.016330891330891332, 2.4975024975024975E-4, 0.014221889221889222, 7.215007215007215E-4},
		    	{0.01898014750588173, 0.018477887335113272, 0.02155753522429882, 0.02400274921356631, 0.008194771207274841, 0.11340505961035184, 0.012411113167146898, 0.009820508075814851, 0.013640328848238124, 0.0014935631393904148, 0.008617727140553543, 0.04451611197758334, 0.06455364931666183, 0.18047000978085595, 0.028998916175420972, 0.028959264056676096, 5.154775436834175E-4, 0.13779111263845198, 0.04142324671548283, 0.06069417642549367, 0.09178143752147823, 0.022231621242961748, 0.04022046578022152, 0.0011895635623463481, 0.0055248618784530384, 5.286949165983769E-4},
		    	{0.1317783971443939, 0.0032761234169478267, 0.0023959708571707986, 0.0020536893061463986, 0.17304777272505012, 0.004498557527749254, 0.0012713314752334848, 0.032370055254021805, 0.06004596352256614, 7.334604664808567E-4, 9.290499242090852E-4, 0.09671898684660897, 0.009437191335387023, 8.80152559777028E-4, 0.13515231529020585, 0.05202679575570877, 2.9338418659234265E-4, 0.17783971443939173, 0.0230306586474989, 0.04317637279350643, 0.040046941469854776, 3.911789154564569E-4, 0.004498557527749254, 2.444868221602856E-4, 0.0036184049679722262, 2.444868221602856E-4},
		    	{0.0044964028776978415, 0.0044964028776978415, 0.00539568345323741, 0.0044964028776978415, 0.00539568345323741, 0.0044964028776978415, 0.0044964028776978415, 0.00539568345323741, 0.0044964028776978415, 0.0044964028776978415, 0.0044964028776978415, 0.0044964028776978415, 0.0044964028776978415, 0.0044964028776978415, 0.0044964028776978415, 0.00539568345323741, 0.0044964028776978415, 0.0044964028776978415, 0.0044964028776978415, 0.00539568345323741, 0.8812949640287769, 0.00539568345323741, 0.0044964028776978415, 0.00539568345323741, 0.0044964028776978415, 0.0044964028776978415},
		    	{0.10308531089470864, 0.010258107213765718, 0.024865273707100313, 0.033390060193501624, 0.22561532885821436, 0.011896883174182976, 0.015694431313226814, 0.012101730169235133, 0.10138350508965996, 0.0018751378777851313, 0.015127162711543915, 0.01974409883079638, 0.02945069490403706, 0.02650404966751757, 0.10138350508965996, 0.013992625508178123, 6.460559074721881E-4, 0.01878289370016703, 0.07889760801739623, 0.07989032807034131, 0.02018530774321641, 0.009706596073240679, 0.013094450222180202, 1.733320727364407E-4, 0.031861586461189376, 3.9393652894645614E-4},
		    	{0.09845728964052919, 0.020638762268185125, 0.035932233589214393, 0.013185051002984454, 0.10870243062257791, 0.019287591500987394, 0.007023118383346449, 0.05388350235341282, 0.09723975114701035, 0.00267264547357793, 0.007438863234791905, 0.01717917118294258, 0.02348958410666825, 0.015219231168985435, 0.083921067870347, 0.038560334971566024, 0.00203418016600098, 0.010438165377362692, 0.07069147277613624, 0.19777576504476682, 0.03536800843368127, 0.002687493503986696, 0.03082451112859879, 2.524165169490267E-4, 0.006874638079258786, 2.2272045613149416E-4},
		    	{0.06887022487130859, 0.010295312923327011, 0.011454890273638581, 0.0055377946356001085, 0.10546735302086155, 0.008615551341099974, 0.003500406393931184, 0.31280411812516934, 0.11806014630181523, 0.0015713898672446492, 0.0015605526957464101, 0.014966133839068004, 0.011238146843673802, 0.0056136548360877815, 0.11214305066377675, 0.007867786507721485, 3.793010024383636E-4, 0.038081820644811706, 0.043034408019506906, 0.05342725548631807, 0.0214467623950149, 0.001332972094283392, 0.020872392305608237, 6.502302898943376E-5, 0.021230018965050123, 5.635329179084259E-4},
		    	{0.035906167234203555, 0.03185773741959894, 0.04903518728717367, 0.028717366628830874, 0.04309496783957624, 0.0061293984108967085, 0.03658721150208097, 0.0014755959137343927, 0.026598562239878925, 9.837306091562618E-4, 0.0018917896329928112, 0.0968974650018918, 0.03658721150208097, 0.1388951948543322, 0.0035944003026863415, 0.04074914869466515, 2.2701475595913735E-4, 0.1453651153991676, 0.13984108967082862, 0.12750662126371548, 4.161937192584185E-4, 0.00170261066969353, 0.0016269390843738176, 0.0012864169504351116, 0.002497162315550511, 5.297010972379871E-4},
		    	{0.08454198473282443, 0.00200381679389313, 0.001049618320610687, 0.0018129770992366412, 0.6177480916030534, 0.0011450381679389313, 7.633587786259542E-4, 0.0014312977099236641, 0.21125954198473282, 8.587786259541985E-4, 0.0011450381679389313, 9.541984732824427E-4, 0.002480916030534351, 6.679389312977099E-4, 0.049904580152671754, 0.0012404580152671756, 5.725190839694657E-4, 0.0016221374045801526, 0.005057251908396947, 0.0022900763358778627, 0.002099236641221374, 9.541984732824427E-4, 0.0026717557251908397, 4.7709923664122136E-4, 0.004770992366412214, 4.7709923664122136E-4},
		    	{0.18162500698987866, 0.00682212156796958, 0.0051445506906000115, 0.005759660012302187, 0.1719510149303808, 0.004361684281160879, 0.0018453279651065259, 0.17804618911815692, 0.1863781244757591, 0.001677570877369569, 0.0018453279651065259, 0.006933959626460885, 0.009058882737795673, 0.035732259687971814, 0.12044958899513504, 0.006262931275513057, 6.151093217021753E-4, 0.01627243751048482, 0.022535368785997874, 0.017390818095397866, 0.0017334899066152212, 0.0013979757311413073, 0.00631885030475871, 2.795951462282615E-4, 0.009170720796286976, 3.9143320471956607E-4},
		    	{0.11284807034684904, 0.017586712261846604, 0.0986809965803615, 0.006350757205666829, 0.07278944797264289, 0.01172447484123107, 0.0039081582804103565, 0.023937469467513434, 0.09330727894479726, 0.0029311187103077674, 0.004885197850512946, 0.010258915486077186, 0.01074743527112848, 0.0039081582804103565, 0.032730825598436736, 0.2598925256472887, 0.003419638495359062, 0.008793356130923302, 0.01074743527112848, 0.1621885686370298, 0.008304836345872008, 0.005862237420615535, 0.01709819247679531, 0.003419638495359062, 0.011235955056179775, 0.002442598925256473},
		    	{0.10702360391479562, 0.04277489925158319, 0.04214162348877375, 0.028151986183074264, 0.08370754173862982, 0.033045480713874496, 0.014219919401266552, 0.038860103626943004, 0.06810592976396085, 0.006390328151986183, 0.0037996545768566492, 0.023949337938975246, 0.04179620034542314, 0.020782959124928037, 0.10391479562464019, 0.03172135866436385, 0.0016119746689694876, 0.02498560736902706, 0.09884858952216465, 0.1158894645941278, 0.008175014392630972, 0.004778353483016695, 0.05043177892918826, 3.454231433506045E-4, 0.0034542314335060447, 0.001093839953943581},
		    	{0.17531305903398928, 0.009838998211091235, 0.011627906976744186, 0.008050089445438283, 0.40518783542039355, 0.008944543828264758, 0.008050089445438283, 0.017889087656529516, 0.08586762075134168, 0.005366726296958855, 0.011627906976744186, 0.025044722719141325, 0.012522361359570662, 0.004472271914132379, 0.03667262969588551, 0.007155635062611807, 0.004472271914132379, 0.008944543828264758, 0.014311270125223614, 0.016100178890876567, 0.01699463327370304, 0.006261180679785331, 0.009838998211091235, 0.004472271914132379, 0.01520572450805009, 0.06976744186046512}};
		hmmZodiac.b = new double[][] {   sg.generate(53), sg.generate(53), sg.generate(53), sg.generate(53), sg.generate(53),
										sg.generate(53), sg.generate(53), sg.generate(53), sg.generate(53), sg.generate(53),
										sg.generate(53), sg.generate(53), sg.generate(53), sg.generate(53), sg.generate(53),
										sg.generate(53), sg.generate(53), sg.generate(53), sg.generate(53), sg.generate(53),
										sg.generate(53), sg.generate(53), sg.generate(53), sg.generate(53), sg.generate(53), 
										sg.generate(53) };
		hmmZodiac.print();
		// 1000 random restarts, 200 iterations per restart
		hmmZodiac.randomRestartWithoutA(1000, 200, 26, 53, ZdoubleArr);
		// 10000 random restarts, 200 iterations per restart
//		hmmZodiac.randomRestartWithoutA(10000, 200, 26, 53, ZdoubleArr);
		// 100000 random restarts, 200 iterations per restart
//		hmmZodiac.randomRestartWithoutA(100000, 200, 26, 53, ZdoubleArr);
		// 1000000 random restarts, 200 iterations per restart
//		hmmZodiac.randomRestartWithoutA(1000000, 200, 26, 53, ZdoubleArr);
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
//