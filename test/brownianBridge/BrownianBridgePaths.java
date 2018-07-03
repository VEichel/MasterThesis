package brownianBridge;

import java.security.GeneralSecurityException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.TrapezoidIntegrator;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.special.Erf;

import montecarlo.interestrates.BrownianBridgeWithVariance;
import net.finmath.montecarlo.BrownianBridge;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.BrownianMotionInterface;
import net.finmath.montecarlo.RandomVariable;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationInterface;

public class BrownianBridgePaths {

	static int seed = 13412;
	
	private static DecimalFormat formatterValue		= new DecimalFormat(" ##0.000;-##0.000", new DecimalFormatSymbols(Locale.ENGLISH));
	private static int path = 0;
	
	public static void main(String[] args) {
		
		long startTime = System.currentTimeMillis();
		
		RandomVariableInterface start = new RandomVariable(0.0);
		RandomVariableInterface end   = new RandomVariable(0.0);
		int numberOfPaths = 100000;
		TimeDiscretizationInterface timeDiscretization = new TimeDiscretization(10, 100, 0.1);
		BrownianMotion standardMotion = new BrownianMotion(timeDiscretization, 1, numberOfPaths, seed);
		BrownianBridge standardBridge = new BrownianBridge(timeDiscretization, numberOfPaths, seed, start, end);
		
		RandomVariableInterface[] variances = new RandomVariableInterface[timeDiscretization.getNumberOfTimeSteps()];
		for (int i = 0; i < variances.length; i++) {
			variances[i] = new RandomVariable(0.02);
		}
		BrownianBridgeWithVariance varianceBridge = new BrownianBridgeWithVariance(timeDiscretization, standardMotion, variances);
		
		RandomVariableInterface[] bridgeValue = new RandomVariableInterface[timeDiscretization.getNumberOfTimes()];
		RandomVariableInterface[] bridgeValueWithVariance = new RandomVariableInterface[timeDiscretization.getNumberOfTimes()];
		bridgeValue[0] = new RandomVariable(0.0);
		bridgeValueWithVariance[0] = new RandomVariable(0.0);
		for (int timeIndex = 0; timeIndex < bridgeValue.length - 1; timeIndex++) {
			bridgeValue[timeIndex+1] = bridgeValue[timeIndex].add(standardBridge.getIncrement(timeIndex,0));
			bridgeValueWithVariance[timeIndex+1] = bridgeValueWithVariance[timeIndex].add(varianceBridge.getIncrement(timeIndex,0));
		}
		
		int n = 20;
		int m = 30;
		System.out.println("Variance: ");
		System.out.println(bridgeValue[2].getVariance());
		System.out.println("with extra: " + bridgeValueWithVariance[n].getVariance());
		System.out.println();
		System.out.println("Covariance: " + bridgeValueWithVariance[m].mult(bridgeValueWithVariance[n]).getAverage());
		System.out.println();
		System.out.println();
		System.out.println("Increment Variance: ");
		System.out.println(standardBridge.getBrownianIncrement(6, 0).getVariance());
		System.out.println(varianceBridge.getBrownianIncrement(6, 0).getVariance());
		System.out.println();
		System.out.println();
		System.out.println(bridgeValue[6].getAverage());
		System.out.println(bridgeValueWithVariance[6].getAverage());
		System.out.println();
		System.out.println(standardBridge.getBrownianIncrement(5, 0).getAverage());
		System.out.println(varianceBridge.getBrownianIncrement(5, 0).getAverage());
		System.out.println();
		
		System.out.println("Test 2:");
		System.out.println("---------");
		//System.out.println(cov(timeDiscretization.getTime(6),timeDiscretization.getTime(5),timeDiscretization.getTime(0),timeDiscretization.getTime(timeDiscretization.getNumberOfTimeSteps())));
		//System.out.println(bridgeValue[6].mult(bridgeValue[5]).getAverage());
		
		long startTime9 = System.currentTimeMillis();
		double varianceOfBB = 0.0;
		for(int temp=0;temp<10000;temp++) {
		varianceOfBB = 0.0;
		for (int i = 1; i < n+1; i++) {
			for (int j = 1; j < n+1; j++) {
				double ti = timeDiscretization.getTime(i-1);
				double ti2 = timeDiscretization.getTime(i);
				double tj= timeDiscretization.getTime(j-1);
				double tj2 = timeDiscretization.getTime(j);
				double T1 = timeDiscretization.getTime(0);
				double T2 = timeDiscretization.getTime(timeDiscretization.getNumberOfTimeSteps());
				varianceOfBB += Math.sqrt(variances[i-1].getAverage())*Math.sqrt(variances[j-1].getAverage())*covOfIncrement(ti, ti2, tj, tj2, T1, T2);
			}
		}}
		long varianceOfBBOldTime = System.currentTimeMillis() - startTime9;
		System.out.println("old calculated: " + varianceOfBB +  "  time needed: " + varianceOfBBOldTime);
		
		double varianceOfBBNew = 0.0;
		long startTime10 = System.currentTimeMillis();
		for(int temp=0;temp<10000;temp++) {
		varianceOfBBNew = 0.0;
		double periodLenght    = timeDiscretization.getTime(timeDiscretization.getNumberOfTimeSteps()) - timeDiscretization.getTime(0);
		for (int i = 1; i < n+1; i++) {
			double ti = timeDiscretization.getTime(i-1);
			double ti2 = timeDiscretization.getTime(i);
			varianceOfBBNew += variances[i-1].getAverage() * ((ti2-ti)-(ti2-ti)*(ti2-ti)/periodLenght);
			for (int j = 1; j < i; j++) {
				double tj= timeDiscretization.getTime(j-1);
				double tj2 = timeDiscretization.getTime(j);
				varianceOfBBNew -= 2* Math.sqrt(variances[i-1].getAverage())*Math.sqrt(variances[j-1].getAverage())*(ti2-ti)*(tj2-tj)/periodLenght;
			}
		}
		}
		long varianceOfBBNewTime = System.currentTimeMillis()-startTime10;
		System.out.println("new calculated: " + varianceOfBBNew + "  time needed: " + varianceOfBBNewTime);
		
		
		
		double covOfBBNew = 0.0;
		long startTime11 = System.currentTimeMillis();
		for(int temp=0;temp<10000;temp++) {
		covOfBBNew = 0.0;
		double periodLenght    = timeDiscretization.getTime(timeDiscretization.getNumberOfTimeSteps()) - timeDiscretization.getTime(0);
		for (int i = 1; i < n+1; i++) {
			double ti = timeDiscretization.getTime(i-1);
			double ti2 = timeDiscretization.getTime(i);
			if(i<m+1) {
				covOfBBNew += variances[i-1].getAverage() * ((ti2-ti)-(ti2-ti)*(ti2-ti)/periodLenght);
			}
			for (int j = 1; j<m+1; j++) {
				if(j != i) {
					double tj= timeDiscretization.getTime(j-1);
					double tj2 = timeDiscretization.getTime(j);
					covOfBBNew -= Math.sqrt(variances[i-1].getAverage())*Math.sqrt(variances[j-1].getAverage())*(ti2-ti)*(tj2-tj)/periodLenght;
				}
			}
		}
		}
		long covOfBBNewTime = System.currentTimeMillis()-startTime11;
		System.out.println("covariance calculated: " + covOfBBNew + "  time needed: " + covOfBBNewTime);
		
		System.out.println();
		
		/*
		long startTime12 = System.currentTimeMillis();
		double result = 0;
		for(int temp =0; temp<100000; temp++) {
			org.apache.commons.math3.distribution.NormalDistribution nv = new org.apache.commons.math3.distribution.NormalDistribution(9, 0.21);
			result = nv.density(10.2);
		}
		System.out.println(result + "      " + (System.currentTimeMillis() - startTime12));
		*/
		/*
		getArray(standardBridge);
		
		long time1 = System.currentTimeMillis()-startTime;
		System.out.println("time1 "  + (time1));
		
		RandomVariableInterface startV = new RandomVariable(5.0);
		RandomVariableInterface endV   = new RandomVariable(21.0);
		getModifiedArray(standardBridge, startV, endV);
		long time2 = System.currentTimeMillis()-startTime-time1;
		System.out.println("time2 "  + time2);
		*/
		
		
		//printPaths(getArray(standardBridge));
		//printPaths(getArray(varianceBridge));
		
		/*
		System.out.println();
		System.out.println("Test für Normalverteilungshypothese");
		
		for (int i = 0; i < numberOfPaths && i<10000; i++) {
			System.out.println(bridgeValueWithVariance[4].get(i));
		}*/
		double residuelVariable =  -0.006473509461259366;
		double varianceVariable = 0.01873444531509433;
		System.out.println();
		System.out.println();
		System.out.println("Test6");
		NormalDistribution nd = new NormalDistribution(0.0, Math.sqrt(varianceVariable));
		double cumulative = nd.cumulativeProbability(residuelVariable);
		
		TrapezoidIntegrator integrator = new TrapezoidIntegrator();
		
		UnivariateFunction f = new UnivariateFunction() {
			
			public double value(double x)
			{
				return nd.cumulativeProbability(x);
			}
		};
		/*double firstPart =  integrator.integrate(1000000000, f, -500000, -residuelVariable);
		System.out.println(firstPart  + residuelVariable );*/
		
		
		double erf = Erf.erf(residuelVariable / Math.sqrt(2 * varianceVariable)); 
		double lt  = 0.5 * residuelVariable * (erf - 1) 
					+ Math.sqrt(varianceVariable / (2 * Math.PI)) * Math.exp(- residuelVariable * residuelVariable / (2 * varianceVariable)); 
		
		double erf2 = 2*cumulative - 1;
		double e   = Math.exp((residuelVariable*residuelVariable)/ (2*varianceVariable));
		//double lt  = (Math.PI * residuelVariable * e * erf - Math.PI * residuelVariable * e + Math.sqrt(2 * Math.PI * varianceVariable)) / (e * 2 * Math.PI);
		
		System.out.println(lt +  residuelVariable);
		TimeDiscretization timeDiscretization2 = new TimeDiscretization(0, 10, varianceVariable);
		BrownianMotion bm = new BrownianMotion(timeDiscretization2, 1, 10000000, 1222143);
		RandomVariableInterface nv = bm.getBrownianIncrement(0, 0);
		
		System.out.println(nv.apply(x -> x>-residuelVariable ? residuelVariable + x : 0.0).getAverage());
	}

	 
	public static double[][] getArray(BrownianMotionInterface bridge) {
		
			double[][] values = new double[bridge.getTimeDiscretization().getNumberOfTimeSteps()+1][bridge.getNumberOfPaths()];
			for (int pathIndex = 0; pathIndex < bridge.getNumberOfPaths(); pathIndex++) {
				values[0][pathIndex] = 0.0;
				for (int timeIndex = 1; timeIndex < bridge.getTimeDiscretization().getNumberOfTimeSteps()+1; timeIndex++) {
					values[timeIndex][pathIndex] = values[timeIndex-1][pathIndex]+bridge.getIncrement(timeIndex-1, 0).get(pathIndex);
				}
			}
			
		return values;	
	}
	
	public static double[][] getModifiedArray(BrownianBridge bridge, RandomVariableInterface start, RandomVariableInterface end) {
		
		double[][] values = getArray(bridge);
		TimeDiscretizationInterface timeDiscretization = bridge.getTimeDiscretization();
		double endTime = timeDiscretization.getTime(timeDiscretization.getNumberOfTimeSteps());
		double timeLenght = endTime-0.0;
		for (int pathIndex = 0; pathIndex < bridge.getNumberOfPaths(); pathIndex++) {
			for (int timeIndex = 0; timeIndex < values.length; timeIndex++) {
				double time = timeDiscretization.getTime(timeIndex);
				values[timeIndex][pathIndex] = (endTime-time)/timeLenght*start.get(pathIndex)+time/timeLenght*end.get(pathIndex)+values[timeIndex][pathIndex];
			}
		}
	return values;	
	}
	
	
	public static void printPaths(double[] values) {
		
		for (int timeIndex = 0; timeIndex < values.length; timeIndex++) {
			System.out.println(formatterValue.format(values[timeIndex]));
		}
		System.out.println();
	}
	
	public static void printPaths(double[][] values) {
		for (int pathIndex = 0; pathIndex < values[0].length; pathIndex++) {
			for (int timeIndex = 0; timeIndex < values.length; timeIndex++) {
				System.out.println(formatterValue.format(values[timeIndex][pathIndex]));
			}
			System.out.println();
		}
	}
	
	public static double cov(double t1, double t2, double T1, double T2) {
		if(t2<t1) {
			double temp = t2;
			t2=t1;
			t1=temp;
		}
		return ((T2-t2)*(t1-T1))/(T2-T1);
	}
	
	public static double covOfIncrement (double ti, double ti2, double tj, double tj2, double T1, double T2) {
		return cov(ti2,tj2,T1,T2)-cov(ti,tj2,T1,T2)-cov(ti2,tj,T1,T2)+cov(ti,tj,T1,T2);
	}

}
