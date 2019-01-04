package montecarlo.interestrates;


import net.finmath.montecarlo.AbstractRandomVariableFactory;
import net.finmath.montecarlo.BrownianMotionInterface;
import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretizationInterface;

import java.io.IOException;
import java.util.Arrays;

public class VariatingBrownianBridge implements BrownianMotionInterface {

	private final TimeDiscretizationInterface						timeDiscretization;


	private AbstractRandomVariableFactory randomVariableFactory ;

	private transient RandomVariableInterface[][]	    brownianIncrements;
	private transient Object							brownianIncrementsLazyInitLock = new Object();

	private final RandomVariableInterface[] variances;
	
	private final BrownianMotionInterface generator;

	public VariatingBrownianBridge(TimeDiscretizationInterface timeDiscretization,
								   BrownianMotionInterface generator, RandomVariableInterface[] variances, AbstractRandomVariableFactory randomVariableFactory) {
		super();
		this.timeDiscretization = timeDiscretization;	
		this.variances = variances;
		this.generator = generator;
		this.randomVariableFactory = randomVariableFactory;
	}
	
	/**
	 * 
	 * @param timeDiscretization
	 * @param generator
	 * @param variances
	 */
	public VariatingBrownianBridge(TimeDiscretizationInterface timeDiscretization,
								   BrownianMotionInterface generator, RandomVariableInterface[] variances) {
		this(timeDiscretization, generator, variances, new RandomVariableFactory());
	}
	
	public VariatingBrownianBridge(TimeDiscretizationInterface timeDiscretization,
								   BrownianMotionInterface generator, RandomVariableInterface[] variances,
								   RandomVariableFactory randomVariableFactory) {
		this.timeDiscretization = timeDiscretization;
		this.variances = variances;
		this.generator = generator;
		this.randomVariableFactory = randomVariableFactory;
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
		

		// Allocate memory
		brownianIncrements = new RandomVariableInterface[getTimeDiscretization().getNumberOfTimeSteps()][getNumberOfFactors()];

		double endTime 		= getTimeDiscretization().getTime(getTimeDiscretization().getNumberOfTimeSteps());
		for(int factor=0; factor<generator.getNumberOfFactors(); factor++) {
			
			// Initialized the bridge to the start point
			RandomVariableInterface brownianBridge	= getRandomVariableForConstant(0.0);
			for(int timeIndex=0; timeIndex<getTimeDiscretization().getNumberOfTimeSteps(); timeIndex++) {
				double currentTime	= getTimeDiscretization().getTime(timeIndex);
				double nextTime		= getTimeDiscretization().getTime(timeIndex+1);
				double alpha		= (nextTime-currentTime)/(endTime-currentTime);
				
				int    generatorTimeIndex = generator.getTimeDiscretization().getTimeIndex(currentTime);

				// Calculate the next point using the "scheme" of the Brownian bridge
				RandomVariableInterface nextRealization = brownianBridge.mult(1.0-alpha).add(generator.getBrownianIncrement(generatorTimeIndex, factor)
						.mult(variances[timeIndex].sqrt()).mult(Math.sqrt(1.0-alpha)));
				
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


	@Override
	public RandomVariableInterface getRandomVariableForConstant(double value) {
		return randomVariableFactory.createRandomVariable(value);
	}

	/* (non-Javadoc)
	 * @see net.finmath.montecarlo.BrownianMotionInterface#getCloneWithModifiedSeed(int)
	 */
	@Override
	public BrownianMotionInterface getCloneWithModifiedSeed(int seed) {
		return new VariatingBrownianBridge(getTimeDiscretization(), this.generator.getCloneWithModifiedSeed(seed), variances);
	}

	/* (non-Javadoc)
	 * @see net.finmath.montecarlo.BrownianMotionInterface#getCloneWithModifiedTimeDiscretization(net.finmath.time.TimeDiscretizationInterface)
	 */
	@Override
	public BrownianMotionInterface getCloneWithModifiedTimeDiscretization(TimeDiscretizationInterface newTimeDiscretization) {
		return new VariatingBrownianBridge(newTimeDiscretization, this.generator, this.variances);
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
		return "BrownianBridge [timeDiscretization=" + timeDiscretization + "]";
	}

	private void readObject(java.io.ObjectInputStream in) throws ClassNotFoundException, IOException {
		in.defaultReadObject();
		// initialization of transients
		brownianIncrementsLazyInitLock = new Object();
	}


	@Override
	public int getNumberOfFactors() {
		return generator.getNumberOfFactors();
	}


	@Override
	public int getNumberOfPaths() {
		return generator.getNumberOfPaths();
	}
}



