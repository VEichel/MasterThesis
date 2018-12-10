package notSorted;

import net.finmath.montecarlo.BrownianMotion;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretization;

import java.util.function.DoubleUnaryOperator;

public class LogNormalSumTest {

	public static void main(String[] args) {
		int numberOfPaths = 1000000;		
		
		
		
		double variance1   = 0.12;
		double variance2   = 0.13;
		double covariance  = 0; //Math.sqrt(variance1)*Math.sqrt(variance2); //same seed
		double a 		   = 13;
		double b           = -13;
		
		TimeDiscretization timeDiscretization1 = new TimeDiscretization(0, 10, variance1);
		BrownianMotion bm1 = new BrownianMotion(timeDiscretization1, 1, numberOfPaths, 12222);
		RandomVariableInterface nv1 = bm1.getBrownianIncrement(0, 0).exp().mult(a);
		
		TimeDiscretization timeDiscretization2 = new TimeDiscretization(0, 10, variance2);
		BrownianMotion bm2 = new BrownianMotion(timeDiscretization2, 1, numberOfPaths, 122222);
		RandomVariableInterface nv2 = bm2.getBrownianIncrement(0, 0).exp().mult(b);
	
		RandomVariableInterface buildSum = nv1.add(nv2);
		/*System.out.println(buildSum.getMin());
		System.out.println(buildSum.getAverage() - 10*log(buildSum.getVariance()));
	    System.out.println(buildSum.getAverage() +"\t"  + buildSum.floor(0.0).getAverage());*/
		/*
		double calculatedVariance       = exp(variance1)*(exp(variance1) - 1)+exp(variance2)*(exp(variance2) - 1)+2*exp((variance1 + variance2) / 2)*
				(exp(Math.sqrt(variance1)*Math.sqrt(variance2)) - 1);
		double calculatedNewSigmaSqaure = log((calculatedVariance / Math.pow(calculatedAverage,2)) + 1);
		double calculatedNewMu			= log(calculatedAverage) - calculatedNewSigmaSqaure / 2;
		*/
		double calculatedAverage        = a*exp(variance1/2) + b*exp(variance2/2);
		double calculated2Moment        = a*a* exp(2 * variance1) + b*b* exp(2 * variance2) + 2*a*b* exp((variance1 + variance2 + 2 * covariance) / 2);
		double calculatedVariance       = a*a* exp(variance1)*(exp(variance1) - 1) + b*b* exp(variance2)*(exp(variance2) - 1) + 2*a*b* exp((variance1 + variance2) / 2)*
				(exp(covariance) - 1);
		System.out.println("Ave + Var: " + calculatedAverage +"\t " + calculatedVariance);
		
		double calculated3Moment        = Math.pow(a, 3)*exp(9.0/2.0 * variance1) + Math.pow(b, 3)*exp(9.0/2.0 * variance2) 
				+ 3.0* Math.pow(a, 2)*b*exp((4.0*variance1+variance2+4.0*covariance)/2.0) + 3.0*a*Math.pow(b, 2)*exp((variance1+4.0*variance2+4.0*covariance)/2.0);
		
		double skewness = (calculated3Moment - 3.0*calculatedVariance*calculatedAverage - Math.pow(calculatedAverage, 3)) / Math.pow(calculatedVariance, 1.5);
		
		//New Outputs
		System.out.println("Discriminate: " + solveThirdDegreePolynominal(-2*pow(calculatedAverage,3) + 3*calculatedAverage*calculated2Moment - calculated3Moment,
				  3*pow(calculatedAverage,2)*calculated2Moment + 3*pow(calculated2Moment,2) - 3*pow(calculatedAverage, 4) + 3*calculatedAverage*calculated2Moment,
				  4*calculatedAverage*pow(calculated2Moment,2) - 3*pow(calculatedAverage,3)*calculated2Moment - 3*pow(calculatedAverage,2)*calculated3Moment,
				  pow(calculated2Moment,3) - pow(calculatedAverage,3)*calculated3Moment)[0]);
		
		
		/** change Z to -Z, if Skeweness<0 **/
		boolean negativeZUsed = false;
		if(skewness<0) {
			negativeZUsed = true;
			System.out.println("-Z used");
			double tempVar1 = variance1;
			double tempA   = a;
			a=-b;
			b=-tempA;
			
			variance1=variance2;
			variance2=tempVar1;
			buildSum = buildSum.mult(-1.0);
			calculatedAverage = -calculatedAverage;
		}
		
		
		
		

		System.out.println("Max/Min: " + buildSum.getMax() +"\t " + buildSum.getMin());
		System.out.println("Quantil. " + buildSum.getQuantile(0.99) +"\t "+ buildSum.getQuantile(0.01));
		/** change Z to Z + |min(Z)|, if min(Z)<0 **/
		double addedMinPre = 0.0;
		boolean addedMinUsed = false;
		if( true ) {
			addedMinPre = 17;
			addedMinUsed = true;
			System.out.println("min(Z) used.");
			calculatedAverage += addedMinPre;
			buildSum = buildSum.add(addedMinPre);
			System.out.println("AddedMin: " + addedMinPre);
		}
		final double addedMin = addedMinPre;
		
		
		
		
		double calculatedNewSigmaSqaure = log(calculatedVariance/Math.pow(calculatedAverage,2) + 1.0);
		System.out.println("calculatedNewSigmaSqaure: " + calculatedNewSigmaSqaure);
				/*log(a*a* exp(2*variance1) + b*b* exp(2*variance2) + 2*a*b* exp((variance1+variance2+covariance)/2)) - 
		 										log(Math.pow(calculatedAverage,2));*/
		double calculatedNewMu			= log(calculatedAverage) - calculatedNewSigmaSqaure / 2.0;
		
		
		System.out.println("Linke Seite: " + (log(-b/a) / Math.sqrt(variance1 + variance2 - covariance)) );
		System.out.println("Rechte Seite:  " + (log(addedMin) - calculatedNewMu) / Math.sqrt(calculatedNewSigmaSqaure) );
		
		
		
		
		TimeDiscretization timeDiscretizationSum = new TimeDiscretization(0, 10, calculatedNewSigmaSqaure);
		BrownianMotion bmSum = new BrownianMotion(timeDiscretizationSum, 1, numberOfPaths, 12222);
		RandomVariableInterface nvSum = (bmSum.getBrownianIncrement(0, 0).add(calculatedNewMu)).exp();
		
		/*System.out.println(buildSum.getAverage());
        System.out.println(buildSum.getVariance());
		System.out.println(buildSum);
		System.out.println(calculatedAverage);*/
		DoubleUnaryOperator operator = null;
		if(negativeZUsed) {
			operator = x -> x < 0 ? - x : 0;
			if(addedMinUsed) {
				operator = x -> x < addedMin ? - ( x - addedMin) : 0;
			}
		} else {
			operator = x -> x > 0 ? x : 0;
			if(addedMinUsed) {
				operator = x -> x > addedMin ? x - addedMin : 0;
			}
		}
		
		
		//if(addedMinUsed)  nvSum = nvSum.sub(addedMin);
		//if(negativeZUsed) nvSum = nvSum.mult(-1.0); 
		
		if(numberOfPaths < 10001) System.out.println("nvSum: " + nvSum);
		if(numberOfPaths < 10001) System.out.println("buidlSum: " + buildSum);
		
		System.out.println("nvSum Ave and Var:    " + nvSum.sub(addedMin).getAverage() +    "\t" + nvSum.getVariance());
		System.out.println("buildSum Ave and VAr: " + buildSum.sub(addedMin).getAverage() + "\t" + buildSum.getVariance());
		
		//System.out.println(nvSum.apply(operator).getAverage());
		System.out.println(buildSum.sub(addedMin).floor(0.0).getAverage());
		System.out.println(nvSum.sub(addedMin).floor(0.0).getAverage());
		
		
		//temp for m-test:
		System.out.println(nvSum.apply(x -> x>0 ? 1 : 0).getAverage());
		System.out.println(buildSum.apply(x -> x>0 ? 1 : 0).getAverage());
	}

	
	
	
	
	private static double log(double a) {
		return Math.log(a);
	}

	private static double pow(double a, double b) {
		return Math.pow(a, b);
	}



	private static double exp(double a) {
		return Math.exp(a);
	}

	private static double[] solveThirdDegreePolynominal(double a, double b, double c, double d) {
		/*double[] roots;
		double discriminate = 18*a*b*c*d - 4*b*b*b*d + b*b*c*c - 4*a*c*c*c - 27*a*a*d*d;
		if(discriminate > 0) {
			roots = new double[3];
			double delta0 = b*b - 3*a*c;
			double delta1 = 2*b*b*b - 9*a*b*c + 27*a*a*d;
			double C1     = Math.pow((delta1 + Math.sqrt(-27*a*a*discriminate)) / 2, 1/3)
			double[0] = (b + )
		}
		return roots;*/
		// Verify preconditions. 
		double TWO_PI = 2.0 * Math.PI; 
		double FOUR_PI = 4.0 * Math.PI; 
		int nRoots; 
		 
		/**
		 * The first real root. 
		 */ 
		double x1; 
		 
		/**
		 * The second real root. 
		 */ 
		double x2; 
		 
		/**
		 * The third real root. 
		 */ 
		double x3; 
		
		 if (a == 0.0) 
		  { 
		  throw new RuntimeException ("Cubic.solve(): a = 0"); 
		  } 
		 
		 // Normalize coefficients. 
		 double denom = a; 
		 a = b/denom; 
		 b = c/denom; 
		 c = d/denom; 
		 
		 // Commence solution. 
		 double a_over_3 = a / 3.0; 
		 double Q = (3*b - a*a) / 9.0; 
		 double Q_CUBE = Q*Q*Q; 
		 double R = (9*a*b - 27*c - 2*a*a*a) / 54.0; 
		 double R_SQR = R*R; 
		 double D = Q_CUBE + R_SQR; 
		 
		 if (D < 0.0) 
		  { 
		  // Three unequal real roots. 
		  nRoots = 3; 
		  double theta = Math.acos (R / Math.sqrt (-Q_CUBE)); 
		  double SQRT_Q = Math.sqrt (-Q); 
		  x1 = 2.0 * SQRT_Q * Math.cos (theta/3.0) - a_over_3; 
		  x2 = 2.0 * SQRT_Q * Math.cos ((theta+TWO_PI)/3.0) - a_over_3; 
		  x3 = 2.0 * SQRT_Q * Math.cos ((theta+FOUR_PI)/3.0) - a_over_3;  
		  } 
		 else if (D > 0.0) 
		  { 
		  // One real root. 
		  nRoots = 1; 
		  double SQRT_D = Math.sqrt (D); 
		  double S = Math.cbrt (R + SQRT_D); 
		  double T = Math.cbrt (R - SQRT_D); 
		  x1 = (S + T) - a_over_3; 
		  x2 = Double.NaN; 
		  x3 = Double.NaN; 
		  } 
		 else 
		  { 
		  // Three real roots, at least two equal. 
		  nRoots = 3; 
		  double CBRT_R = Math.cbrt (R); 
		  x1 = 2*CBRT_R - a_over_3; 
		  x2 = x3 = CBRT_R - a_over_3;  
		  } 
		 double[] roots = new double[nRoots];
		 if(nRoots>0) roots[0] = x1;
		 if(nRoots>1) roots[1] = x2;
		 if(nRoots>2) roots[2] = x3;
		 
		 return roots;
	}

	
	

}
