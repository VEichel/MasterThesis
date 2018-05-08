package montecarlo.interestrate;


import java.io.IOException;
import java.util.Arrays;

import net.finmath.montecarlo.AbstractRandomVariableFactory;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.BrownianMotionInterface;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretizationInterface;

public class BrownianBridgeWithVariance implements BrownianMotionInterface {

	private final TimeDiscretizationInterface						timeDiscretization;

	private final int	numberOfFactors;
	private final int	numberOfPaths;
	private final int	seed;

	private RandomVariableInterface[] start;
	private RandomVariableInterface[] end;

	private AbstractRandomVariableFactory randomVariableFactory;

	private transient RandomVariableInterface[][]	    brownianIncrements;
	private transient Object							brownianIncrementsLazyInitLock = new Object();

	/**
	 * Construct a Brownian bridge, bridging from a given start to a given end.
	 * 
	 * @param timeDiscretization The time discretization used for the Brownian increments.
	 * @param numberOfPaths Number of paths to simulate.
	 * @param seed The seed of the random number generator.
	 * @param start Start value of the Brownian bridge.
	 * @param end End value of the Brownian bridge.
	 */
	public BrownianBridgeWithVariance(TimeDiscretizationInterface timeDiscretization, int numberOfPaths, int seed, RandomVariableInterface[] start, RandomVariableInterface[] end) {
		super();
		this.timeDiscretization = timeDiscretization;
		this.numberOfFactors = start.length;
		this.numberOfPaths = numberOfPaths;
		this.seed = seed;
		this.start = start;
		this.end = end;
	}

	/**
	 * Construct a Brownian bridge, bridging from a given start to a given end.
	 * 
	 * @param timeDiscretization The time discretization used for the Brownian increments.
	 * @param numberOfPaths Number of paths to simulate.
	 * @param seed The seed of the random number generator.
	 * @param start Start value of the Brownian bridge.
	 * @param end End value of the Brownian bridge.
	 */
	public BrownianBridgeWithVariance(TimeDiscretizationInterface timeDiscretization, int numberOfPaths, int seed, RandomVariableInterface start, RandomVariableInterface end) {
		this(timeDiscretization, numberOfPaths, seed, new RandomVariableInterface[] {start}, new RandomVariableInterface[] {end});
	}

	/* (non-Javadoc)
	 * @see net.finmath.montecarlo.BrownianMotionInterface#getBrownianIncrement(int, int)
	 */
	@Override
	public RandomVariableInterface getBrownianIncrement(int timeIndex, int factor) {
		// Thread safe lazy initialization
		synchronized(brownianIncrementsLazyInitLock) {
			if(brownianIncrements == null) doGenerateBrownianMotion();
		}

		/*
		 *  For performance reasons we return directly the stored data (no defensive copy).
		 *  We return an immutable object to ensure that the receiver does not alter the data.
		 */
		return brownianIncrements[timeIndex][factor];
	}

	/**
	 * Lazy initialization of brownianIncrement. Synchronized to ensure thread safety of lazy init.
	 */
	private void doGenerateBrownianMotion() {
		if(brownianIncrements != null) return;	// Nothing to do
		
		BrownianMotion generator = new BrownianMotion(timeDiscretization, numberOfFactors, numberOfPaths, seed);

		// Allocate memory
		brownianIncrements = new RandomVariableInterface[generator.getTimeDiscretization().getNumberOfTimeSteps()][generator.getNumberOfFactors()];

		double endTime 		= getTimeDiscretization().getTime(getTimeDiscretization().getNumberOfTimeSteps());
		for(int factor=0; factor<generator.getNumberOfFactors(); factor++) {
			// The end point
			RandomVariableInterface endOfFactor		= end[factor];
			// Initialized the bridge to the start point
			RandomVariableInterface brownianBridge	= start[factor];
			for(int timeIndex=0; timeIndex<getTimeDiscretization().getNumberOfTimeSteps(); timeIndex++) {
				double currentTime	= getTimeDiscretization().getTime(timeIndex);
				double nextTime		= getTimeDiscretization().getTime(timeIndex+1);
				double alpha		= (nextTime-currentTime)/(endTime-currentTime);

				// Calculate the next point using the "scheme" of the Brownian bridge
				RandomVariableInterface nextRealization = brownianBridge.mult(1.0-alpha).add(endOfFactor.mult(alpha)).add(generator.getBrownianIncrement(timeIndex, factor).mult(0.00).mult(Math.sqrt(1-alpha)));
				
				// Store the increment
				brownianIncrements[timeIndex][factor] = nextRealization.sub(brownianBridge);
				
				// Update the bridge to the current point
				brownianBridge = nextRealization;
			}
		}
	}

	/* (non-Javadoc)
	 * @see net.finmath.montecarlo.BrownianMotionInterface#getTimeDiscretization()
	 */
	@Override
	public TimeDiscretizationInterface getTimeDiscretization() {
		return timeDiscretization;
	}

	/* (non-Javadoc)
	 * @see net.finmath.montecarlo.BrownianMotionInterface#getNumberOfFactors()
	 */
	@Override
	public int getNumberOfFactors() {
		return numberOfFactors;
	}

	/* (non-Javadoc)
	 * @see net.finmath.montecarlo.BrownianMotionInterface#getNumberOfPaths()
	 */
	@Override
	public int getNumberOfPaths() {
		return numberOfPaths;
	}

	@Override
	public RandomVariableInterface getRandomVariableForConstant(double value) {
		return randomVariableFactory.createRandomVariable(value);
	}

	/* (non-Javadoc)
	 * @see net.finmath.montecarlo.BrownianMotionInterface#getCloneWithModifiedSeed(int)
	 */
	@Override
	public BrownianMotionInterface getCloneWithModifiedSeed(int seed) {
		return new BrownianBridgeWithVariance(timeDiscretization, numberOfPaths, seed, start, end);
	}

	/* (non-Javadoc)
	 * @see net.finmath.montecarlo.BrownianMotionInterface#getCloneWithModifiedTimeDiscretization(net.finmath.time.TimeDiscretizationInterface)
	 */
	@Override
	public BrownianMotionInterface getCloneWithModifiedTimeDiscretization(TimeDiscretizationInterface newTimeDiscretization) {
		return new BrownianBridgeWithVariance(newTimeDiscretization, getNumberOfFactors(), seed, start, end);
	}

	@Override
	public RandomVariableInterface getIncrement(int timeIndex, int factor) {
		return getBrownianIncrement(timeIndex, factor);
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "BrownianBridge [timeDiscretization=" + timeDiscretization
				+ ", numberOfFactors=" + numberOfFactors + ", numberOfPaths="
				+ numberOfPaths + ", seed=" + seed + ", start="
				+ Arrays.toString(start) + ", end=" + Arrays.toString(end)
				+ "]";
	}

	private void readObject(java.io.ObjectInputStream in) throws ClassNotFoundException, IOException {
		in.defaultReadObject();
		// initialization of transients
		brownianIncrementsLazyInitLock = new Object();
	}
}



