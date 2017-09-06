package HMM;
import java.text.*;
import java.util.List;
/** This class implements a Hidden Markov Model, as well as
    the Baum-Welch Algorithm for training HMMs.
    @author Holger Wunsch (wunsch@sfs.nphil.uni-tuebingen.de) 
*/
public class HMM2 {
  /** number of states */
  public int numStates;

  /** size of output vocabulary */
  public int sigmaSize;

  /** initial state probabilities */
  public double pi[];

  /** transition probabilities */
  public double a[][];

  /** emission probabilities */
  public double b[][];

  /** initializes an HMM.
      @param numStates number of states
      @param sigmaSize size of output vocabulary 
  */
  public HMM2(int numStates, int sigmaSize) {
    this.numStates = numStates;
    this.sigmaSize = sigmaSize;

    pi = new double[numStates];
    a = new double[numStates][numStates];
    b = new double[numStates][sigmaSize];
  }

  /** implementation of the Baum-Welch Algorithm for HMMs.
      @param o the training set
      @param steps the number of steps
  */
  public void train(double[] o, int steps) {
    int T = o.length;
    double[][] alpha;
    double[][] beta;
    double pi1[] = new double[numStates];
    double a1[][] = new double[numStates][numStates];
    double b1[][] = new double[numStates][sigmaSize];

    for (int s = 0; s < steps; s++) {
      /* calculation of Forward- und Backward Variables from the
	 current model */
      alpha = alphaPass(o);
      beta = betaPass(o);
      /* re-estimation of initial state probabilities */
      for (int i = 0; i < numStates; i++)
    	  pi1[i] = gamma(i, 0, o, alpha, beta);
      /* re-estimation of transition probabilities */ 
      for (int i = 0; i < numStates; i++) {
		for (int j = 0; j < numStates; j++) {
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
      for (int i = 0; i < numStates; i++) {
		for (int k = 0; k < sigmaSize; k++) {
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
  
  public double[][] alphaPass(double[] o) {
    int T = o.length;
    double[][] alpha = new double[numStates][T];
    /* initialization (time 0) */
    for (int i = 0; i < numStates; i++)
      alpha[i][0] = pi[i] * b[i][(int) o[0]];
    /* induction */
    for (int t = 0; t <= T-2; t++) {
      for (int j = 0; j < numStates; j++) {
		alpha[j][t+1] = 0;
		for (int i = 0; i < numStates; i++)
		  alpha[j][t+1] += (alpha[i][t] * a[i][j]);
		  alpha[j][t+1] *= b[j][(int) o[t+1]];
	    }
    }
    return alpha;
  }

  public double[][] betaPass(double[] o) {
    int T = o.length;
    double[][] beta = new double[numStates][T];  
    /* initialization (time 0) */
    for (int i = 0; i < numStates; i++)
      beta[i][T-1] = 1;
    /* induction */
    for (int t = T - 2; t >= 0; t--) {
      for (int i = 0; i < numStates; i++) {
	beta[i][t] = 0;
	for (int j = 0; j < numStates; j++)
	  beta[i][t] += (beta[j][t+1] * a[i][j] * b[j][(int) o[t+1]]);
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
    for (int k = 0; k < numStates; k++)
      denom += (alpha[k][t] * beta[k][t]);
    return  num / denom;
  }

  public double gamma(int i, int t, double[] o, double[][] alpha, double[][] beta) {
    double num = alpha[i][t] * beta[i][t];
    double denom = 0;
    for (int j = 0; j < numStates; j++)
      denom += alpha[j][t] * beta[j][t];
    return  num / denom;
  }

  public void print() {
    DecimalFormat fmt = new DecimalFormat();
    fmt.setMinimumFractionDigits(5);
    fmt.setMaximumFractionDigits(5);
    
    for (int i = 0; i < numStates; i++)
    	System.out.println("pi(" + i + ") = " + fmt.format(pi[i]));
    	System.out.println();

    for (int i = 0; i < numStates; i++) {
    	for (int j = 0; j < numStates; j++)
    		System.out.print("a(" + i + "," + j + ") = " + 
    							fmt.format(a[i][j]) + "  ");
    	System.out.println();
    }

    System.out.println();
    for (int i = 0; i < numStates; i++) {
    	for (int k = 0; k < sigmaSize; k++)
    		System.out.print("b(" + i + "," + k + ") = " + 
    							fmt.format(b[i][k]) + "  ");
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
}
