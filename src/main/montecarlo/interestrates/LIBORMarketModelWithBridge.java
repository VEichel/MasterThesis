package montecarlo.interestrates;

import montecarlo.interestrates.modelplugins.AbstractLiborCovarianceModelWithInterpolation;
import net.finmath.exception.CalculationException;
import net.finmath.functions.AnalyticFormulas;
import net.finmath.marketdata.model.AnalyticModelInterface;
import net.finmath.marketdata.model.curves.DiscountCurveInterface;
import net.finmath.marketdata.model.curves.ForwardCurveInterface;
import net.finmath.marketdata.model.volatilities.AbstractSwaptionMarketData;
import net.finmath.marketdata.products.Swap;
import net.finmath.marketdata.products.SwapAnnuity;
import net.finmath.montecarlo.AbstractRandomVariableFactory;
import net.finmath.montecarlo.BrownianMotionInterface;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.interestrate.LIBORMarketModel.CalibrationItem;
import net.finmath.montecarlo.interestrate.LIBORMarketModelInterface;
import net.finmath.montecarlo.interestrate.modelplugins.AbstractLIBORCovarianceModel;
import net.finmath.montecarlo.interestrate.modelplugins.AbstractLIBORCovarianceModelParametric;
import net.finmath.montecarlo.interestrate.products.AbstractLIBORMonteCarloProduct;
import net.finmath.montecarlo.interestrate.products.SwaptionAnalyticApproximation;
import net.finmath.montecarlo.interestrate.products.SwaptionSimple;
import net.finmath.montecarlo.model.AbstractModel;
import net.finmath.montecarlo.process.AbstractProcessInterface;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.RegularSchedule;
import net.finmath.time.ScheduleInterface;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationInterface;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

public class LIBORMarketModelWithBridge extends AbstractModel implements LIBORMarketModelInterface {

	public enum Driftapproximation	{ EULER, LINE_INTEGRAL, PREDICTOR_CORRECTOR }
	public enum Measure				{ SPOT, TERMINAL }
	public enum StateSpace			{ NORMAL, LOGNORMAL }

	private final TimeDiscretizationInterface		liborPeriodDiscretization;

	private String							forwardCurveName;
	private AnalyticModelInterface			curveModel;

	private ForwardCurveInterface			forwardRateCurve;
	private DiscountCurveInterface			discountCurve;

	private final AbstractRandomVariableFactory	randomVariableFactory;
	private 	  AbstractLiborCovarianceModelWithInterpolation	covarianceModel;

	private AbstractSwaptionMarketData		swaptionMarketData;

	private Driftapproximation	driftApproximationMethod	= Driftapproximation.EULER;
	private Measure				measure						= Measure.SPOT;
	private StateSpace			stateSpace					= StateSpace.LOGNORMAL;
	private double				liborCap					= 1E5;

	// This is a cache of the integrated covariance.
	private double[][][]	integratedLIBORCovariance;
	private final Object	integratedLIBORCovarianceLazyInitLock = new Object();

	// Cache for the numeraires, needs to be invalidated if process changes
	private final ConcurrentHashMap<Integer, RandomVariableInterface>	numeraires;
	private AbstractProcessInterface									numerairesProcess = null;

	
	public enum InterpolationScheme { LINEAR, LOGLINEAR }
	private InterpolationScheme interpolationScheme = InterpolationScheme.LOGLINEAR;
	
	private RandomVariableInterface[][] brownianBridgeValues;

	private final BrownianMotionInterface interpolationDriver;
	
	/**
	 * 
	 * @param liborPeriodDiscretization
	 * @param analyticModel
	 * @param forwardRateCurve
	 * @param discountCurve
	 * @param randomVariableFactory
	 * @param covarianceModel
	 * @param calibrationItems
	 * @param properties
	 * @throws CalculationException
	 */
	public LIBORMarketModelWithBridge(
			TimeDiscretizationInterface							liborPeriodDiscretization,
			AnalyticModelInterface								analyticModel,
			ForwardCurveInterface								forwardRateCurve,
			DiscountCurveInterface								discountCurve,
			AbstractRandomVariableFactory						randomVariableFactory,
			AbstractLiborCovarianceModelWithInterpolation		covarianceModel,
			CalibrationItem[]									calibrationItems,
			Map<String, ?>									 	properties,

			BrownianMotionInterface 							interpolationDriver
			) throws CalculationException {

		// Set some properties
		if(properties != null && properties.containsKey("measure"))					measure		= Measure.valueOf(((String)properties.get("measure")).toUpperCase());
		if(properties != null && properties.containsKey("stateSpace"))				stateSpace	= StateSpace.valueOf(((String)properties.get("stateSpace")).toUpperCase());
		if(properties != null && properties.containsKey("liborCap"))				liborCap	= (Double)properties.get("liborCap");

		Map<String,Object> calibrationParameters = null;
		if(properties != null && properties.containsKey("calibrationParameters"))	calibrationParameters	= (Map<String,Object>)properties.get("calibrationParameters");

		this.liborPeriodDiscretization	= liborPeriodDiscretization;
		this.curveModel					= analyticModel;
		this.forwardRateCurve		= forwardRateCurve;
		this.discountCurve			= discountCurve;
		this.randomVariableFactory	= randomVariableFactory;
		this.covarianceModel		= covarianceModel;

		this.interpolationDriver = interpolationDriver;

		double[] times = new double[liborPeriodDiscretization.getNumberOfTimeSteps()];
		for(int i=0; i<times.length; i++) times[i] = liborPeriodDiscretization.getTime(i);

		// Perform calibration, if data is given
		if(calibrationItems != null && calibrationItems.length > 0) {

			AbstractLiborCovarianceModelWithInterpolation covarianceModelParametric = null;
			try {
				covarianceModelParametric = (AbstractLiborCovarianceModelWithInterpolation)covarianceModel;
			}
			catch(Exception e) {
				throw new ClassCastException("Calibration is currently restricted to parametric covariance models (AbstractLIBORCovarianceModelParametric).");
			}

			// @TODO Should be more elegant. Convert array for constructor
			AbstractLIBORMonteCarloProduct[]	calibrationProducts		= new AbstractLIBORMonteCarloProduct[calibrationItems.length];
			double[]							calibrationTargetValues	= new double[calibrationItems.length];
			double[]							calibrationWeights		= new double[calibrationItems.length];
			for(int i=0; i<calibrationTargetValues.length; i++) {
				calibrationProducts[i]		= calibrationItems[i].calibrationProduct;
				calibrationTargetValues[i]	= calibrationItems[i].calibrationTargetValue;
				calibrationWeights[i]		= calibrationItems[i].calibrationWeight;
			}
			this.covarianceModel    = (AbstractLiborCovarianceModelWithInterpolation) covarianceModelParametric.getCloneCalibrated(this, calibrationProducts, calibrationTargetValues, calibrationWeights, calibrationParameters);
		}

		numeraires = new ConcurrentHashMap<Integer, RandomVariableInterface>();
		//if(interpolationDriver.getTimeDiscretization().equals(covarianceModel.getTimeDiscretization()) ) throw new IllegalArgumentException("Time discretization for the Brownian Motion used in the Bridge does not match the one in Covariance Model!");
	}

	
	
	public LIBORMarketModelWithBridge(
			TimeDiscretizationInterface						liborPeriodDiscretization,
			ForwardCurveInterface 							forwardRateCurve, 
			DiscountCurveInterface 							discountCurve,
			AbstractLiborCovarianceModelWithInterpolation 	covarianceModel,
			AbstractSwaptionMarketData 						swaptionMarketData, 
			Map<String, Object> 							properties,
			BrownianMotionInterface 						interpolationDriver
			) throws CalculationException {
		this(
				liborPeriodDiscretization,
				null, 
				forwardRateCurve,
				discountCurve,
				new RandomVariableFactory(), 
				covarianceModel,
				getCalibrationItems(
						liborPeriodDiscretization,
						forwardRateCurve,
						swaptionMarketData,
						// Condition under which we use analytic approximation
						(properties == null || properties.get("stateSpace") == null || ((String)properties.get("stateSpace")).toUpperCase().equals(StateSpace.LOGNORMAL.name()))
						&& AbstractLIBORCovarianceModelParametric.class.isAssignableFrom(covarianceModel.getClass())
						),
				properties,
				interpolationDriver
				);	
		}



	private static CalibrationItem[] getCalibrationItems(TimeDiscretizationInterface liborPeriodDiscretization, ForwardCurveInterface forwardCurve, AbstractSwaptionMarketData swaptionMarketData, boolean isUseAnalyticApproximation) {
		if(swaptionMarketData == null) return null;

		TimeDiscretizationInterface	optionMaturities		= swaptionMarketData.getOptionMaturities();
		TimeDiscretizationInterface	tenor					= swaptionMarketData.getTenor();
		double						swapPeriodLength		= swaptionMarketData.getSwapPeriodLength();

		ArrayList<CalibrationItem> calibrationItems = new ArrayList<CalibrationItem>();
		for(int exerciseIndex=0; exerciseIndex<=optionMaturities.getNumberOfTimeSteps(); exerciseIndex++) {
			for(int tenorIndex=0; tenorIndex<=tenor.getNumberOfTimeSteps()-exerciseIndex; tenorIndex++) {

				// Create a swaption
				double exerciseDate	= optionMaturities.getTime(exerciseIndex);
				double swapLength	= tenor.getTime(tenorIndex);

				if(liborPeriodDiscretization.getTimeIndex(exerciseDate) < 0) continue;
				if(liborPeriodDiscretization.getTimeIndex(exerciseDate+swapLength) <= liborPeriodDiscretization.getTimeIndex(exerciseDate)) continue;

				int numberOfPeriods = (int)(swapLength / swapPeriodLength);

				double[] fixingDates      = new double[numberOfPeriods];
				double[] paymentDates     = new double[numberOfPeriods];
				double[] swapTenorTimes   = new double[numberOfPeriods+1];

				for(int periodStartIndex=0; periodStartIndex<numberOfPeriods; periodStartIndex++) {
					fixingDates[periodStartIndex] = exerciseDate + periodStartIndex * swapPeriodLength;
					paymentDates[periodStartIndex] = exerciseDate + (periodStartIndex+1) * swapPeriodLength;
					swapTenorTimes[periodStartIndex] = exerciseDate + periodStartIndex * swapPeriodLength;
				}
				swapTenorTimes[numberOfPeriods] = exerciseDate + numberOfPeriods * swapPeriodLength;


				// Swaptions swap rate
				ScheduleInterface swapTenor = new RegularSchedule(new TimeDiscretization(swapTenorTimes));
				double swaprate = Swap.getForwardSwapRate(swapTenor, swapTenor, forwardCurve, null);

				// Set swap rates for each period
				double[] swaprates        = new double[numberOfPeriods];
				for(int periodStartIndex=0; periodStartIndex<numberOfPeriods; periodStartIndex++) {
					swaprates[periodStartIndex] = swaprate;
				}

				if(isUseAnalyticApproximation) {
					AbstractLIBORMonteCarloProduct swaption = new SwaptionAnalyticApproximation(swaprate, swapTenorTimes, SwaptionAnalyticApproximation.ValueUnit.VOLATILITY);
					double impliedVolatility = swaptionMarketData.getVolatility(exerciseDate, swapLength, swaptionMarketData.getSwapPeriodLength(), swaprate);

					calibrationItems.add(new CalibrationItem(swaption, impliedVolatility, 1.0));
				}
				else {
					AbstractLIBORMonteCarloProduct swaption = new SwaptionSimple(swaprate, swapTenorTimes, SwaptionSimple.ValueUnit.VALUE);

					double forwardSwaprate		= Swap.getForwardSwapRate(swapTenor, swapTenor, forwardCurve);
					double swapAnnuity 			= SwapAnnuity.getSwapAnnuity(swapTenor, forwardCurve);
					double impliedVolatility	= swaptionMarketData.getVolatility(exerciseDate, swapLength, swaptionMarketData.getSwapPeriodLength(), swaprate);

					double targetValue = AnalyticFormulas.blackModelSwaptionValue(forwardSwaprate, impliedVolatility, exerciseDate, swaprate, swapAnnuity);

					calibrationItems.add(new CalibrationItem(swaption, targetValue, 1.0));
				}
			}
		}

		return calibrationItems.toArray(new CalibrationItem[calibrationItems.size()]);
	}
	
	
	
	

	@Override
	public RandomVariableInterface getNumeraire(double time) throws CalculationException {
		int liborTimeIndex = getLiborPeriodIndex(time);
		int timeIndex	   = Math.max(getTimeIndex(time), - getTimeIndex(time) - 1);
		
		if(liborTimeIndex < 0) {
			// Interpolation of Numeraire: log linear interpolation.
			int upperIndex = -liborTimeIndex-1;
			int lowerIndex = upperIndex-1;
			if(lowerIndex < 0) throw new IllegalArgumentException("Numeraire requested for time " + time + ". Unsupported");

			/*RandomVariableInterface numeraire = getNumeraire(getLiborPeriod(upperIndex)).div( 
					getInterpolatedLibor(timeIndex , timeIndex).mult(getLiborPeriod(upperIndex) - time).add(1.0) );*/
			/*
			 * Adjust for discounting, i.e. funding or collateralization
			 */
			
			
			
			double upperTime    				= getLiborPeriod(upperIndex);   //T_i+1
			
			RandomVariableInterface numeraire   = getNumeraire(upperTime).div(getLIBOR(time, time, upperTime).mult(upperTime - time).add(1.0));

			
			/*
			RandomVariableInterface bridge		= getBrownianBridge(lowerIndex, time);
			RandomVariableInterface evaluationTimeScalingFactor = covarianceModel.getEvaluationTimeScalingFactor(timeIndex);
			double lowerTime = getLiborPeriod(lowerIndex);
			double alpha = (upperTime - time) / (upperTime - lowerTime);
			
			double analyticLiborShortPeriod				= getForwardRateCurve().getForward(getAnalyticModel(), time, upperTime -time);
			double analyticLibor					 	= getForwardRateCurve().getForward(getAnalyticModel(), lowerTime, upperTime -lowerTime);
			double analyticInterpolatedOnePlusLiborDt		= Math.exp(Math.log(1 + analyticLibor * (upperTime - lowerTime)) * alpha);
			double analyticOnePlusLiborDt					= (1 + analyticLiborShortPeriod * (upperTime -time));
			double adjustment = analyticOnePlusLiborDt / analyticInterpolatedOnePlusLiborDt;
			RandomVariableInterface numeraire2 = (getUnAdjustedNumeraire(upperTime).log().mult(1 - alpha).add(getUnAdjustedNumeraire(lowerTime).log().mult(alpha))
					.sub(bridge.mult(evaluationTimeScalingFactor))).exp();
		    numeraire = numeraire2;*/
			
			
			if(discountCurve != null) {
				// This includes a control for zero bonds
				double deterministicNumeraireAdjustment = numeraire.invert().getAverage() / discountCurve.getDiscountFactor(curveModel, time);
				numeraire = numeraire.mult(deterministicNumeraireAdjustment);
				//System.out.println("inLMM numeraireAdj interpol " + time + "\t" + deterministicNumeraireAdjustment);
				
				/*
				RandomVariableInterface bridge		= getBrownianBridge(lowerIndex, time);
				RandomVariableInterface evaluationTimeScalingFactor = covarianceModel.getEvaluationTimeScalingFactor(timeIndex);
				double lowerTime = getLiborPeriod(lowerIndex);
				double alpha = (upperTime - time) / (upperTime - lowerTime);
				
				double analyticLiborShortPeriod				= getForwardRateCurve().getForward(getAnalyticModel(), time, upperTime -time);
				double analyticLibor					 	= getForwardRateCurve().getForward(getAnalyticModel(), lowerTime, upperTime -lowerTime);
				double analyticInterpolatedOnePlusLiborDt		= Math.exp(Math.log(1 + analyticLibor * (upperTime - lowerTime)) * alpha);
				double analyticOnePlusLiborDt					= (1 + analyticLiborShortPeriod * (upperTime -time));
				double adjustment = analyticOnePlusLiborDt / analyticInterpolatedOnePlusLiborDt;
				RandomVariableInterface numeraire2 = (getUnAdjustedNumeraire(upperTime).log().mult(1 - alpha).add(getUnAdjustedNumeraire(lowerTime).log().mult(alpha))
						.sub(bridge.mult(evaluationTimeScalingFactor))).exp().div(adjustment);
				System.out.println(numeraire);
				System.out.println(numeraire2);
				*/
			}

			return numeraire;
		}

		/*
		 * Calculate the numeraire, when time is part of liborPeriodDiscretization
		 */

		/*
		 * Check if numeraire cache is values (i.e. process did not change)
		 */
		if(getProcess() != numerairesProcess) {
			numeraires.clear();
			numerairesProcess = getProcess();
		}

		/*
		 * Check if numeraire is part of the cache
		 */
		RandomVariableInterface numeraire = numeraires.get(liborTimeIndex);
		if(numeraire == null) {
			/*
			 * Calculate the numeraire for timeIndex
			 */

			// Initialize to 1.0
			numeraire = getRandomVariableForConstant(1.0);


			// Get the start and end of the product
			int firstLiborIndex, lastLiborIndex;

			if(measure == Measure.TERMINAL) {
				firstLiborIndex	= getLiborPeriodIndex(time);
				if(firstLiborIndex < 0) {
					throw new CalculationException("Simulation time discretization not part of forward rate tenor discretization.");
				}

				lastLiborIndex 	= liborPeriodDiscretization.getNumberOfTimeSteps()-1;
			}
			else if(measure == Measure.SPOT) {
				// Spot measure
				firstLiborIndex	= 0;
				lastLiborIndex	= getLiborPeriodIndex(time)-1;
			}
			else {
				throw new CalculationException("Numeraire not implemented for specified measure.");
			}

			// The product 
			for(int liborIndex = firstLiborIndex; liborIndex<=lastLiborIndex; liborIndex++) {
				RandomVariableInterface libor = getLIBOR(getTimeIndex(Math.min(time,liborPeriodDiscretization.getTime(liborIndex))), liborIndex);

				double periodLength = liborPeriodDiscretization.getTimeStep(liborIndex);

				if(measure == Measure.SPOT) {
					numeraire = numeraire.accrue(libor, periodLength);
				}
				else {
					numeraire = numeraire.discount(libor, periodLength);
				}
			}
			numeraires.put(liborTimeIndex, numeraire);
		}

		/*
		 * Adjust for discounting, i.e. funding or collateralization
		 */
		if(discountCurve != null) {
			// This includes a control for zero bonds
			double deterministicNumeraireAdjustment = numeraire.invert().getAverage() / discountCurve.getDiscountFactor(curveModel, time);
			numeraire = numeraire.mult(deterministicNumeraireAdjustment);
		}
		return numeraire;
	}
	
	
	@Override
	public RandomVariableInterface[] getInitialState() {
		double[] liborInitialStates = new double[liborPeriodDiscretization.getNumberOfTimeSteps()];
		for(int timeIndex=0; timeIndex<liborPeriodDiscretization.getNumberOfTimeSteps(); timeIndex++) {
			double rate = forwardRateCurve.getForward(curveModel, liborPeriodDiscretization.getTime(timeIndex), liborPeriodDiscretization.getTimeStep(timeIndex));
			liborInitialStates[timeIndex] = (stateSpace == StateSpace.LOGNORMAL) ? Math.log(Math.max(rate,0)) : rate;
		}

		RandomVariableInterface[] initialStateRandomVariable = new RandomVariableInterface[getNumberOfComponents()];
		for(int componentIndex=0; componentIndex<getNumberOfComponents(); componentIndex++) {
			initialStateRandomVariable[componentIndex] = getRandomVariableForConstant(liborInitialStates[componentIndex]);
		}
		return initialStateRandomVariable;
	}
	
	
	
	@Override
	public RandomVariableInterface[] getDrift(int timeIndex, RandomVariableInterface[] realizationAtTimeIndex, RandomVariableInterface[] realizationPredictor) {
		double	time				= getTime(timeIndex);
		int		firstLiborIndex		= this.getLiborPeriodIndex(time)+1;
		if(firstLiborIndex<0) firstLiborIndex = -firstLiborIndex-1 + 1;

		RandomVariableInterface		zero	= getRandomVariableForConstant(0.0);

		// Allocate drift vector and initialize to zero (will be used to sum up drift components)
		RandomVariableInterface[]	drift = new RandomVariableInterface[getNumberOfComponents()];
		for(int componentIndex=firstLiborIndex; componentIndex<getNumberOfComponents(); componentIndex++) {
			drift[componentIndex] = zero;
		}

		RandomVariableInterface[]	covarianceFactorSums	= new RandomVariableInterface[getNumberOfFactors()];
		for(int factorIndex=0; factorIndex<getNumberOfFactors(); factorIndex++) {
			covarianceFactorSums[factorIndex] = zero;
		}

		if(measure == Measure.SPOT) {
			// Calculate drift for the component componentIndex (starting at firstLiborIndex, others are zero)
			for(int componentIndex=firstLiborIndex; componentIndex<getNumberOfComponents(); componentIndex++) {
				double						periodLength	= liborPeriodDiscretization.getTimeStep(componentIndex);
				RandomVariableInterface		libor			= realizationAtTimeIndex[componentIndex];
				RandomVariableInterface		oneStepMeasureTransform = getRandomVariableForConstant(periodLength).discount(libor, periodLength);

				if(stateSpace == StateSpace.LOGNORMAL) oneStepMeasureTransform = oneStepMeasureTransform.mult(libor);

				RandomVariableInterface[]	factorLoading   	= getFactorLoading(timeIndex, componentIndex, realizationAtTimeIndex);
				for(int factorIndex=0; factorIndex<getNumberOfFactors(); factorIndex++) {
					covarianceFactorSums[factorIndex] = covarianceFactorSums[factorIndex].add(oneStepMeasureTransform.mult(factorLoading[factorIndex]));
					drift[componentIndex] = drift[componentIndex].addProduct(covarianceFactorSums[factorIndex], factorLoading[factorIndex]);
				}
			}
		}
		else if(measure == Measure.TERMINAL) {
			// Calculate drift for the component componentIndex (starting at firstLiborIndex, others are zero)
			for(int componentIndex=getNumberOfComponents()-1; componentIndex>=firstLiborIndex; componentIndex--) {
				double					periodLength	= liborPeriodDiscretization.getTimeStep(componentIndex);
				RandomVariableInterface libor			= realizationAtTimeIndex[componentIndex];
				RandomVariableInterface oneStepMeasureTransform = getRandomVariableForConstant(periodLength).discount(libor, periodLength);

				if(stateSpace == StateSpace.LOGNORMAL) oneStepMeasureTransform = oneStepMeasureTransform.mult(libor);

				RandomVariableInterface[]	factorLoading   	= getFactorLoading(timeIndex, componentIndex, realizationAtTimeIndex);
				for(int factorIndex=0; factorIndex<getNumberOfFactors(); factorIndex++) {
					drift[componentIndex] = drift[componentIndex].addProduct(covarianceFactorSums[factorIndex], factorLoading[factorIndex]);
					covarianceFactorSums[factorIndex] = covarianceFactorSums[factorIndex].sub(oneStepMeasureTransform.mult(factorLoading[factorIndex]));
				}
			}
		}

		if(stateSpace == StateSpace.LOGNORMAL) {
			// Drift adjustment for log-coordinate in each component
			for(int componentIndex=firstLiborIndex; componentIndex<getNumberOfComponents(); componentIndex++) {
				RandomVariableInterface		variance		= covarianceModel.getCovariance(getTime(timeIndex), componentIndex, componentIndex, realizationAtTimeIndex);
				drift[componentIndex] = drift[componentIndex].addProduct(variance, -0.5);
			}
		}

		return drift;
	}

	@Override
	public	RandomVariableInterface[]	getFactorLoading(int timeIndex, int componentIndex, RandomVariableInterface[] realizationAtTimeIndex)
	{
		return covarianceModel.getFactorLoading(getTime(timeIndex), getLiborPeriod(componentIndex), realizationAtTimeIndex);
	}

	@Override
	public RandomVariableInterface applyStateSpaceTransform(int componentIndex, RandomVariableInterface randomVariable) {
		RandomVariableInterface value = randomVariable;

		if(stateSpace == StateSpace.LOGNORMAL)	value = value.exp();

		if(!Double.isInfinite(liborCap)) value = value.cap(liborCap);

		return value;
	}

	@Override
	public RandomVariableInterface applyStateSpaceTransformInverse(int componentIndex, RandomVariableInterface randomVariable) {
		RandomVariableInterface value = randomVariable;

		if(stateSpace == StateSpace.LOGNORMAL)	value = value.log();

		return value;
	}

	/* (non-Javadoc)
	 * @see net.finmath.montecarlo.model.AbstractModelInterface#getRandomVariableForConstant(double)
	 */
	@Override
	public RandomVariableInterface getRandomVariableForConstant(double value) {
		return randomVariableFactory.createRandomVariable(value);
	}

	/**
	 * @return Returns the driftApproximationMethod.
	 */
	public Driftapproximation getDriftApproximationMethod() {
		return driftApproximationMethod;
	}
	
	
	

	@Override
	public RandomVariableInterface getLIBOR(double time, double periodStart, double periodEnd)
			throws CalculationException {
		int periodStartIndex    = getLiborPeriodIndex(periodStart);
		int periodEndIndex      = getLiborPeriodIndex(periodEnd);

		// The forward rates are provided on fractional tenor discretization points using linear interpolation. See ISBN 0470047224.

		if(periodEndIndex < 0) {
			int nextPeriodEndIndex = -periodEndIndex - 1;
			double nextPeriodEnd   = getLiborPeriod(nextPeriodEndIndex);
			
			RandomVariableInterface longOnePlusLibordt = getRandomVariableForConstant(1.0);
			if(nextPeriodEnd != periodStart) {
				longOnePlusLibordt         = getLIBOR(time, periodStart, nextPeriodEnd).mult(nextPeriodEnd - periodStart).add(1.0);
			}
			RandomVariableInterface interpolatedOnePlusLibordt = getInterpolatedLibor(getTimeIndex(time), getTimeIndex(periodEnd)).mult(nextPeriodEnd - periodEnd).add(1.0);
			
			return longOnePlusLibordt.div(interpolatedOnePlusLibordt).sub(1.0).div(periodEnd - periodStart);
		}

		if(periodStartIndex < 0) {
			int nextPeriodStartIndex = -periodStartIndex - 1;
			double nextPeriodStart	 = getLiborPeriod(nextPeriodStartIndex);
			
			RandomVariableInterface longOnePlusLibordt = getRandomVariableForConstant(1.0);
			if(nextPeriodStart != periodEnd) {
				longOnePlusLibordt         = getLIBOR(time, nextPeriodStart, periodEnd).mult(periodEnd - nextPeriodStart).add(1.0);
			}
			
			RandomVariableInterface interpolatedOnePlusLibordt = getInterpolatedLibor(getTimeIndex(time), getTimeIndex(periodStart)).mult(nextPeriodStart - periodStart).add(1.0);
			return longOnePlusLibordt.mult(interpolatedOnePlusLibordt).sub(1.0).div(periodEnd - periodStart);
		}

		if(periodStartIndex < 0 || periodEndIndex < 0) throw new AssertionError("LIBOR requested outside libor discretization points and interpolation was not performed.");

		// If time is beyond fixing, use the fixing time.
		time = Math.min(time, periodStart);
		int timeIndex           = getTimeIndex(time);

		// If time is not part of the discretization, use the latest available point.
		if(timeIndex < 0) {
			timeIndex = -timeIndex-2;
			//			double timeStep = getTimeDiscretization().getTimeStep(timeIndex);
			//			return getLIBOR(getTime(timeIndex), periodStart, periodEnd).mult((getTime(timeIndex+1)-time)/timeStep).add(getLIBOR(getTime(timeIndex+1), periodStart, periodEnd).mult((time-getTime(timeIndex))/timeStep));
		}

		// If this is a model primitive then return it
		if(periodStartIndex+1==periodEndIndex) return getLIBOR(timeIndex, periodStartIndex);

		// The requested LIBOR is not a model primitive. We need to calculate it (slow!)
		RandomVariableInterface accrualAccount = null; //=randomVariableFactory.createRandomVariable(1.0);

		// Calculate the value of the forward bond
		for(int periodIndex = periodStartIndex; periodIndex<periodEndIndex; periodIndex++)
		{
			double subPeriodLength = getLiborPeriod(periodIndex+1) - getLiborPeriod(periodIndex);
			RandomVariableInterface liborOverSubPeriod = getLIBOR(timeIndex, periodIndex);

			accrualAccount = accrualAccount == null ? liborOverSubPeriod.mult(subPeriodLength).add(1.0) : accrualAccount.accrue(liborOverSubPeriod, subPeriodLength);
		}

		RandomVariableInterface libor = accrualAccount.sub(1.0).div(periodEnd - periodStart);

		return libor;
	}
	
	@Override
	public RandomVariableInterface getLIBOR(int timeIndex, int liborIndex) throws CalculationException {
		timeIndex = Math.min(timeIndex,getTimeIndex(getLiborPeriod(liborIndex)));
		return getProcessValue(timeIndex, liborIndex);
	}
	/**
	 * 
	 * @param evaluationTimeIndex 
	 * @param processTimeIndex    t
	 * @return L(t,m(t)+1;evaluationTimeIndex), where m(t) is the smallest LiborIndex greater than t.
	 * @throws CalculationException
	 */
	public RandomVariableInterface getInterpolatedLibor(int evaluationTimeIndex, int processTimeIndex) throws CalculationException {
		//insure Index is on TimeDiscretization
		evaluationTimeIndex = Math.max(evaluationTimeIndex, -evaluationTimeIndex - 1);
		processTimeIndex    = Math.max(processTimeIndex, -processTimeIndex - 1);
		double processTime        = getTime(processTimeIndex);
		int    previousLiborIndex = getLiborPeriodIndex(processTime); 
		
		//if(previousLiborIndex>=0) throw new UnsupportedOperationException("This method is only for inner period LIBORs!");
		RandomVariableInterface evaluationTimeScalingFactor = covarianceModel.getEvaluationTimeScalingFactor(evaluationTimeIndex);
		
		if(previousLiborIndex<0)	previousLiborIndex = (-previousLiborIndex-1)-1;	//i
		double previousLiborTime 			= getLiborPeriod(previousLiborIndex);   //T_i
		double nextLiborTime     			= getLiborPeriod(previousLiborIndex+1); //T_i+1
		double periodLenght     			= nextLiborTime - previousLiborTime;
		double shortPeriodLenght			= nextLiborTime	- processTime;
		double alpha            			= (nextLiborTime - processTime) / periodLenght;
		RandomVariableInterface startLibor	= getLIBOR(evaluationTimeIndex, previousLiborIndex); //L(T_i,T_i+1)
		RandomVariableInterface bridge		= getBrownianBridge(previousLiborIndex, processTime);
		RandomVariableInterface libor;
		/*if(interpolationScheme == InterpolationScheme.LINEAR) {
			libor = startLibor.mult(periodLenght).add(1.0).mult(alpha).add(1-alpha).mult(bridge.mult(evaluationTimeScalingFactor)).sub(1.0).div(shortPeriodLenght);
		}*/
		
		if(interpolationScheme == InterpolationScheme.LOGLINEAR) {
			libor = (startLibor.mult(periodLenght).add(1.0).log().mult(alpha).add(bridge.mult(evaluationTimeScalingFactor))).exp().sub(1.0).div(shortPeriodLenght);
		}
		else { throw new UnsupportedOperationException("InterpolationScheme not supported!"); }
		
		if(getForwardRateCurve() != null) {
			double analyticLiborShortPeriod				= getForwardRateCurve().getForward(getAnalyticModel(), processTime, shortPeriodLenght);
			double analyticLibor					 	= getForwardRateCurve().getForward(getAnalyticModel(), previousLiborTime, periodLenght);
			double analyticInterpolatedOnePlusLiborDt		= Math.exp(Math.log(1 + analyticLibor * periodLenght) * alpha);
			double analyticOnePlusLiborDt					= (1 + analyticLiborShortPeriod * (shortPeriodLenght));
			double adjustment = analyticOnePlusLiborDt / analyticInterpolatedOnePlusLiborDt;
			libor = libor.mult(shortPeriodLenght).add(1.0).mult(adjustment).sub(1.0).div(shortPeriodLenght);
		}
	
		return libor;
	}
	
	/**
	 * 
	 * @param liborIndex the previous liborIndex of processTime.
	 * @param time processTime.
	 * @return BB(time)
	 */
	public RandomVariableInterface getBrownianBridge(int liborIndex, double time) {
		
		if(brownianBridgeValues == null) {
			brownianBridgeValues = new RandomVariableInterface[getNumberOfLibors()][];
		}
		//Lazy Init for each Brownian Bridge
		if(brownianBridgeValues[liborIndex] == null) {
		
			double liborPeriodStart = getLiborPeriod(liborIndex);
			double liborPeriodEnd   = getLiborPeriod(liborIndex + 1);
			TimeDiscretizationInterface completeTimeDiscretization = getTimeDiscretization();
			TimeDiscretizationInterface bridgeDiscretization    = new TimeDiscretization(
					Arrays.copyOfRange(completeTimeDiscretization.getAsDoubleArray(), getTimeIndex(liborPeriodStart), getTimeIndex(liborPeriodEnd) + 1));
			
			BrownianBridgeWithVariance brownianBridge = new BrownianBridgeWithVariance(bridgeDiscretization, interpolationDriver, covarianceModel.getVarianceForInterpolationPeriod(liborIndex), randomVariableFactory);
					//new BrownianBridgeWithVariance(bridgeDiscretization, getProcess().getNumberOfPaths(), getRandomVariableForConstant(0.0), getRandomVariableForConstant(0.0),
					//covarianceModel.getVarianceForInterpolationPeriod(liborIndex));
			
			brownianBridgeValues[liborIndex] = new RandomVariableInterface[bridgeDiscretization.getNumberOfTimes()];
			
			RandomVariableInterface brownianBridgeValue = getRandomVariableForConstant(0.0);
			brownianBridgeValues[liborIndex][0] = brownianBridgeValue;
			int bridgeTimeIndex = bridgeDiscretization.getTimeIndex(time);
			for (int timeIndex = 0; timeIndex < brownianBridgeValues[liborIndex].length - 1; timeIndex++) {
				brownianBridgeValue = brownianBridge.getIncrement(timeIndex, 0).add(brownianBridgeValue);
				brownianBridgeValues[liborIndex][timeIndex + 1] =  brownianBridgeValue;
			}
			return brownianBridgeValues[liborIndex][bridgeTimeIndex];
		}
		int bridgeTimeIndex = getTimeIndex(time) - getTimeIndex(getLiborPeriod(liborIndex));
		return brownianBridgeValues[liborIndex][bridgeTimeIndex];
	}
	
	@Override
	public int getNumberOfComponents() {
		return liborPeriodDiscretization.getNumberOfTimeSteps();
	}

	@Override
	public int getNumberOfLibors()
	{
		// This is just a synonym to number of components
		return getNumberOfComponents();
	}

	/* (non-Javadoc)
	 * @see net.finmath.montecarlo.interestrate.LIBORMarketModelInterface#getLiborPeriod(int)
	 */
	@Override
	public double getLiborPeriod(int timeIndex) {
		if(timeIndex >= liborPeriodDiscretization.getNumberOfTimes() || timeIndex < 0) {
			throw new ArrayIndexOutOfBoundsException("Index for LIBOR period discretization out of bounds: " + timeIndex + ".");
		}
		return liborPeriodDiscretization.getTime(timeIndex);
	}

	/* (non-Javadoc)
	 * @see net.finmath.montecarlo.interestrate.LIBORMarketModelInterface#getLiborPeriodIndex(double)
	 */
	@Override
	public int getLiborPeriodIndex(double time) {
		return liborPeriodDiscretization.getTimeIndex(time);
	}

	/* (non-Javadoc)
	 * @see net.finmath.montecarlo.interestrate.LIBORMarketModelInterface#getLiborPeriodDiscretization()
	 */
	@Override
	public TimeDiscretizationInterface getLiborPeriodDiscretization() {
		return liborPeriodDiscretization;
	}

	/**
	 * @return Returns the measure.
	 */
	public Measure getMeasure() {
		return measure;
	}

	/* (non-Javadoc)
	 * @see net.finmath.montecarlo.interestrate.LIBORMarketModelInterface#getIntegratedLIBORCovariance()
	 */
	@Override
	public double[][][] getIntegratedLIBORCovariance() {
		synchronized (integratedLIBORCovarianceLazyInitLock) {
			if(integratedLIBORCovariance == null) {
				TimeDiscretizationInterface liborPeriodDiscretization = getLiborPeriodDiscretization();
				TimeDiscretizationInterface simulationTimeDiscretization = getCovarianceModel().getTimeDiscretization();

				integratedLIBORCovariance = new double[simulationTimeDiscretization.getNumberOfTimeSteps()][liborPeriodDiscretization.getNumberOfTimeSteps()][liborPeriodDiscretization.getNumberOfTimeSteps()];
				for(int timeIndex = 0; timeIndex < simulationTimeDiscretization.getNumberOfTimeSteps(); timeIndex++) {
					double dt = simulationTimeDiscretization.getTime(timeIndex+1) - simulationTimeDiscretization.getTime(timeIndex);
					RandomVariableInterface[][] factorLoadings = new RandomVariableInterface[liborPeriodDiscretization.getNumberOfTimeSteps()][];
					// Prefetch factor loadings
					for(int componentIndex = 0; componentIndex < liborPeriodDiscretization.getNumberOfTimeSteps(); componentIndex++) {
						factorLoadings[componentIndex] = getCovarianceModel().getFactorLoading(timeIndex, componentIndex, null);
					}
					for(int componentIndex1 = 0; componentIndex1 < liborPeriodDiscretization.getNumberOfTimeSteps(); componentIndex1++) {
						RandomVariableInterface[] factorLoadingOfComponent1 = factorLoadings[componentIndex1];
						// Sum the libor cross terms (use symmetry)
						for(int componentIndex2 = componentIndex1; componentIndex2 < liborPeriodDiscretization.getNumberOfTimeSteps(); componentIndex2++) {
							double integratedLIBORCovarianceValue = 0.0;
							if(getLiborPeriod(componentIndex1) > getTime(timeIndex)) {
								RandomVariableInterface[] factorLoadingOfComponent2 = factorLoadings[componentIndex2];
								for(int factorIndex = 0; factorIndex < getNumberOfFactors(); factorIndex++) {
									integratedLIBORCovarianceValue += factorLoadingOfComponent1[factorIndex].get(0) * factorLoadingOfComponent2[factorIndex].get(0) * dt;
								}
							}
							integratedLIBORCovariance[timeIndex][componentIndex1][componentIndex2] = integratedLIBORCovarianceValue;
						}
					}
				}

				// Integrate over time (i.e. sum up).
				for(int timeIndex = 1; timeIndex < simulationTimeDiscretization.getNumberOfTimeSteps(); timeIndex++) {
					double[][] prevIntegratedLIBORCovariance = integratedLIBORCovariance[timeIndex-1];
					double[][] thisIntegratedLIBORCovariance = integratedLIBORCovariance[timeIndex];
					for(int componentIndex1 = 0; componentIndex1 < liborPeriodDiscretization.getNumberOfTimeSteps(); componentIndex1++) {
						for(int componentIndex2 = componentIndex1; componentIndex2 < liborPeriodDiscretization.getNumberOfTimeSteps(); componentIndex2++) {
							thisIntegratedLIBORCovariance[componentIndex1][componentIndex2] = prevIntegratedLIBORCovariance[componentIndex1][componentIndex2] + thisIntegratedLIBORCovariance[componentIndex1][componentIndex2];
							thisIntegratedLIBORCovariance[componentIndex2][componentIndex1] = thisIntegratedLIBORCovariance[componentIndex1][componentIndex2];
						}
					}
				}
			}
		}

		return integratedLIBORCovariance;
	}
	
	@Override
	public Object clone() {
		try {
			Map<String, Object> properties = new HashMap<String, Object>();
			properties.put("measure",		measure.name());
			properties.put("stateSpace",	stateSpace.name());
			return new LIBORMarketModelWithBridge(getLiborPeriodDiscretization(), getAnalyticModel(),
					getForwardRateCurve(), getDiscountCurve(), randomVariableFactory, covarianceModel,
					new CalibrationItem[0], properties, interpolationDriver);
		} catch (CalculationException e) {
			return null;
		}
	}
	
	@Override
	public LIBORMarketModelWithBridge getCloneWithModifiedCovarianceModel(AbstractLIBORCovarianceModel covarianceModel) {
		LIBORMarketModelWithBridge model = (LIBORMarketModelWithBridge)this.clone();
		try {
			model.covarianceModel = (AbstractLiborCovarianceModelWithInterpolation) covarianceModel;
		} catch (Exception exception) {
			throw exception;
		}
		return model;
	}

	@Override
	public LIBORMarketModelWithBridge getCloneWithModifiedData(Map<String, Object> dataModified) throws CalculationException {
		TimeDiscretizationInterface		liborPeriodDiscretization	= this.liborPeriodDiscretization;
		AnalyticModelInterface			analyticModel				= this.curveModel;
		ForwardCurveInterface			forwardRateCurve			= this.forwardRateCurve;
		DiscountCurveInterface			discountCurve				= this.discountCurve;
		AbstractLiborCovarianceModelWithInterpolation	covarianceModel				= this.covarianceModel;
		AbstractSwaptionMarketData		swaptionMarketData			= null;		// No recalibration, unless new swaption data is specified
		Map<String, Object>				properties					= new HashMap<String, Object>();
		properties.put("measure",		measure.name());
		properties.put("stateSpace",	stateSpace.name());
		BrownianMotionInterface interpolationDriver = this.interpolationDriver;

		if(dataModified.containsKey("liborPeriodDiscretization")) {
			liborPeriodDiscretization = (TimeDiscretizationInterface)dataModified.get("liborPeriodDiscretization");
		}
		if(dataModified.containsKey("forwardRateCurve")) {
			forwardRateCurve = (ForwardCurveInterface)dataModified.get("forwardRateCurve");
		}
		if(dataModified.containsKey("discountCurve")) {
			discountCurve = (DiscountCurveInterface)dataModified.get("discountCurve");
		}
		if(dataModified.containsKey("forwardRateShift")) {
			throw new RuntimeException("Forward rate shift clone currently disabled.");
		}
		if(dataModified.containsKey("covarianceModel")) {
			covarianceModel = (AbstractLiborCovarianceModelWithInterpolation)dataModified.get("covarianceModel");
		}
		if(dataModified.containsKey("swaptionMarketData")) {
			swaptionMarketData = (AbstractSwaptionMarketData)dataModified.get("swaptionMarketData");
		}
		if(dataModified.containsKey("interpolationDriver")) {
			interpolationDriver = (BrownianMotionInterface)dataModified.get("interpolationDriver");
		}
		
		LIBORMarketModelWithBridge newModel = new LIBORMarketModelWithBridge(liborPeriodDiscretization, forwardRateCurve, discountCurve, covarianceModel, swaptionMarketData, properties, interpolationDriver);
		newModel.curveModel = analyticModel;
		return newModel;
	}

	@Override
	public String toString() {
		return "LIBORMarketModelWithBridgeInterpolation [liborPeriodDiscretization="
				+ liborPeriodDiscretization + ", forwardCurveName="
				+ forwardCurveName + ", curveModel=" + curveModel
				+ ", forwardRateCurve=" + forwardRateCurve + ", discountCurve="
				+ discountCurve + ", covarianceModel=" + covarianceModel
				+ ", driftApproximationMethod=" + driftApproximationMethod
				+ ", measure=" + measure + ", stateSpace=" + stateSpace + "]";
	}
	
	@Override
	public AnalyticModelInterface getAnalyticModel() {
		return curveModel;
	}

	@Override
	public DiscountCurveInterface getDiscountCurve() {
		return discountCurve;
	}

	@Override
	public ForwardCurveInterface getForwardRateCurve() {
		return forwardRateCurve;
	}

	/**
	 * Return the swaption market data used for calibration (if any, may be null).
	 * 
	 * @return The swaption market data used for calibration (if any, may be null).
	 */
	public AbstractSwaptionMarketData getSwaptionMarketData() {
		return swaptionMarketData;
	}

	/* (non-Javadoc)
	 * @see net.finmath.montecarlo.interestrate.LIBORMarketModelInterface#getCovarianceModel()
	 */
	@Override
	public AbstractLIBORCovarianceModel getCovarianceModel() {
		return covarianceModel;
	}
	
	public RandomVariableInterface getUnAdjustedNumeraire(double time) throws CalculationException {
		int liborTimeIndex = getLiborPeriodIndex(time);
		if(liborTimeIndex < 0) {
			throw new UnsupportedOperationException("Unadjusted Numeraire only available for Tenortimes. Requested " + time + " not supported!" );
		}
		RandomVariableInterface unAdjustedNumeraire = numeraires.get(liborTimeIndex);
		if(unAdjustedNumeraire == null) {
			getNumeraire(time);
			unAdjustedNumeraire = numeraires.get(liborTimeIndex);
		}
		return unAdjustedNumeraire;
	}
}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
