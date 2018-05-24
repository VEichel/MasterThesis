package montecarlo.interestrates;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import montecarlo.interestrates.modelplugins.AbstractLiborCovarianceModelWithInterpolation;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.AnalyticModelInterface;
import net.finmath.marketdata.model.curves.DiscountCurveInterface;
import net.finmath.marketdata.model.curves.ForwardCurveInterface;
import net.finmath.montecarlo.AbstractRandomVariableFactory;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORMarketModel.CalibrationItem;
import net.finmath.montecarlo.interestrate.modelplugins.AbstractLIBORCovarianceModel;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationInterface;

public class LiborMarketModelWithBridgeInterpolation extends LIBORMarketModel {

	enum InterpolationScheme { LINEAR, LOGLINEAR }
	InterpolationScheme interpolationScheme = InterpolationScheme.LOGLINEAR;
	
	private RandomVariableInterface[][] brownianBridgeValues;

	private int seed;
	
	AbstractLiborCovarianceModelWithInterpolation covarianceModel;
	

	public LiborMarketModelWithBridgeInterpolation(TimeDiscretizationInterface liborPeriodDiscretization,
			AnalyticModelInterface analyticModel, ForwardCurveInterface forwardRateCurve,
			DiscountCurveInterface discountCurve, AbstractLiborCovarianceModelWithInterpolation covarianceModel,
			CalibrationItem[] calibrationItems, Map<String, ?> properties, int seed) throws CalculationException {
		super(liborPeriodDiscretization, analyticModel, forwardRateCurve, discountCurve, covarianceModel, calibrationItems,
				properties);
		this.covarianceModel = covarianceModel;
		this.seed = seed;
	}
	
	public LiborMarketModelWithBridgeInterpolation(
			TimeDiscretizationInterface			liborPeriodDiscretization,
			AnalyticModelInterface				analyticModel,
			ForwardCurveInterface				forwardRateCurve,
			DiscountCurveInterface				discountCurve,
			AbstractLiborCovarianceModelWithInterpolation		covarianceModel,
			CalibrationItem[]					calibrationItems,
			Map<String, ?>						properties
			) throws CalculationException {
		this(liborPeriodDiscretization, analyticModel, forwardRateCurve, discountCurve, covarianceModel, calibrationItems, properties, /*seed*/ 231317);
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

			RandomVariableInterface numeraire = getNumeraire(getLiborPeriod(upperIndex)).div( 
					getInterpolatedLibor(timeIndex , timeIndex).mult(getLiborPeriod(upperIndex) - time).add(1.0) );
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
		return super.getLIBOR(timeIndex, liborIndex);
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
		
		if(previousLiborIndex<0)	previousLiborIndex = (-previousLiborIndex-1)-1;	//i
		double previousLiborTime 			= getLiborPeriod(previousLiborIndex);   //T_i
		double nextLiborTime     			= getLiborPeriod(previousLiborIndex+1); //T_i+1
		double periodLenght     			= nextLiborTime - previousLiborTime;
		double shortPeriodLenght			= nextLiborTime	- processTime;
		double alpha            			= (nextLiborTime - processTime) / periodLenght;
		RandomVariableInterface startLibor	= getLIBOR(evaluationTimeIndex, previousLiborIndex); //L(T_i,T_i+1)
		RandomVariableInterface bridge		= getBrownianBridge(previousLiborIndex, processTime);
		RandomVariableInterface libor;
		if(interpolationScheme == InterpolationScheme.LINEAR) {
			libor = startLibor.mult(periodLenght).add(1.0).mult(alpha).add(1-alpha).add(bridge).sub(1.0).div(shortPeriodLenght);
		}
		
		if(interpolationScheme == InterpolationScheme.LOGLINEAR) {
			libor = startLibor.mult(periodLenght).add(1.0).log().mult(alpha).exp().add(bridge).sub(1.0).div(shortPeriodLenght);
		}
		else { throw new UnsupportedOperationException("InterpolationScheme not supported!"); }
		
		if(getDiscountCurve() != null) {
			double analyticLibor				= getForwardRateCurve().getForward(getAnalyticModel(), previousLiborTime, shortPeriodLenght);
			double analyticLiborShortPeriod		= getForwardRateCurve().getForward(getAnalyticModel(), previousLiborTime, periodLenght);
			double analyticInterpolatedOnePlusLiborDt		= (1 + analyticLiborShortPeriod * periodLenght) / Math.exp(Math.log(1 + analyticLiborShortPeriod * periodLenght) * alpha);
			double analyticOnePlusLiborDt					= (1 + analyticLibor * (processTime - previousLiborTime));
			double adjustment = analyticOnePlusLiborDt / analyticInterpolatedOnePlusLiborDt;
			libor = libor.mult(shortPeriodLenght).add(1.0).div(adjustment).sub(1.0).div(shortPeriodLenght);
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
			
			BrownianBridgeWithVariance brownianBridge = new BrownianBridgeWithVariance(bridgeDiscretization, getProcess().getNumberOfPaths(), seed, getRandomVariableForConstant(0.0), getRandomVariableForConstant(0.0),
					covarianceModel.getVarianceForInterpolationPeriod(liborIndex));
			seed = seed + 1; //for future Bridges.
			brownianBridgeValues[liborIndex] = new RandomVariableInterface[bridgeDiscretization.getNumberOfTimes()];
			
			RandomVariableInterface brownianBridgeValue = getRandomVariableForConstant(0.0);
			brownianBridgeValues[liborIndex][0] = brownianBridgeValue;
			int bridgeTimeIndex = bridgeDiscretization.getTimeIndex(time);
			for (int timeIndex = 0; timeIndex < brownianBridgeValues[liborIndex].length - 1; timeIndex++) {
				brownianBridgeValue = brownianBridge.getIncrement(timeIndex, 0)/*.mult(covarianceModel.getFactorLoadingForInterpolation(globalTimeIndex)[0])*/.add(brownianBridgeValue);
				brownianBridgeValues[liborIndex][timeIndex + 1] = brownianBridgeValue;
			}
			return brownianBridgeValues[liborIndex][bridgeTimeIndex];
		}
		int bridgeTimeIndex = getTimeIndex(time) - getTimeIndex(getLiborPeriod(liborIndex));
		return brownianBridgeValues[liborIndex][bridgeTimeIndex];
	}
	
	@Override
	public LiborMarketModelWithBridgeInterpolation getCloneWithModifiedCovarianceModel(AbstractLIBORCovarianceModel covarianceModel) {
		Map<String, Object> properties = new HashMap<String, Object>();
		properties.put("measure",		measure.name());
		properties.put("stateSpace",	stateSpace.name());
		try {
			return new LiborMarketModelWithBridgeInterpolation(liborPeriodDiscretization, getAnalyticModel(), getForwardRateCurve(), getDiscountCurve(),(AbstractLiborCovarianceModelWithInterpolation) covarianceModel, new CalibrationItem[0],properties);
		} catch (CalculationException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}
	
}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
