package brownianBridge;

import java.security.GeneralSecurityException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;

import net.finmath.montecarlo.BrownianBridge;
import net.finmath.montecarlo.RandomVariable;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationInterface;

public class BrownianBridgePaths {

	private static DecimalFormat formatterValue		= new DecimalFormat(" ##0.000;-##0.000", new DecimalFormatSymbols(Locale.ENGLISH));
	private static int path = 0;
	
	public static void main(String[] args) {
		
		long startTime = System.currentTimeMillis();
		
		RandomVariableInterface start = new RandomVariable(0.0);
		RandomVariableInterface end   = new RandomVariable(0.0);
		int numberOfPaths = 1;
		TimeDiscretizationInterface timeDiscretization = new TimeDiscretization(0.0, 100, 0.1);
		
		BrownianBridge standardBridge = new BrownianBridge(timeDiscretization, numberOfPaths, /* seed */ 123, start, end);
		getArray(standardBridge);
		
		long time1 = System.currentTimeMillis()-startTime;
		System.out.println("time1 "  + (time1));
		
		RandomVariableInterface startV = new RandomVariable(5.0);
		RandomVariableInterface endV   = new RandomVariable(21.0);
		getModifiedArray(standardBridge, startV, endV);
		long time2 = System.currentTimeMillis()-startTime-time1;
		System.out.println("time2 "  + time2);
		
		printPaths(getArray(standardBridge));
	}

	 
	public static double[][] getArray(BrownianBridge bridge) {
		
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

}
