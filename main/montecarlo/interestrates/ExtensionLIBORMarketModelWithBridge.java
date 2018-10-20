package montecarlo.interestrates;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;


import montecarlo.interestrates.modelplugins.AbstractLiborCovarianceModelWithInterpolation;
import montecarlo.interestrates.modelplugins.LiborCovarianceModelWithInterpolation;
import montecarlo.interestrates.modelplugins.LiborCovarianceModelWithInterpolation.EvaluationTimeScalingScheme;
import montecarlo.interestrates.modelplugins.LiborCovarianceModelWithInterpolation.InterpolationVarianceScheme;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.AnalyticModelInterface;
import net.finmath.marketdata.model.curves.DiscountCurveInterface;
import net.finmath.marketdata.model.curves.ForwardCurveInterface;
import net.finmath.marketdata.model.volatilities.AbstractSwaptionMarketData;
import net.finmath.montecarlo.BrownianMotionInterface;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORMarketModel.CalibrationItem;
import net.finmath.montecarlo.interestrate.modelplugins.AbstractLIBORCovarianceModel;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationInterface;




public class ExtensionLIBORMarketModelWithBridge extends LIBORMarketModel{

	public enum InterpolationScheme { LINEAR, LOGLINEAR }
	private InterpolationScheme interpolationScheme = InterpolationScheme.LOGLINEAR;

	private RandomVariableInterface[][] brownianBridgeValues;

	private final BrownianMotionInterface interpolationDriver;
	
	
	
	public ExtensionLIBORMarketModelWithBridge() throws CalculationException {
	}
	
	public static ExtensionLIBORMarketModelWithBridge getLIBORMarketModelExtendedWithBridge(LIBORMarketModel liborMarketModel,
			double[] interpolationParameters, double[] evaluationTimeScalingParameters,
			InterpolationVarianceScheme interpolationVarianceScheme, EvaluationTimeScalingScheme evaluationTimeScalingScheme) {
		
		AbstractLiborCovarianceModelWithInterpolation newCovarianceModel = new LiborCovarianceModelWithInterpolation(liborMarketModel.getCovarianceModel(), evaluationTimeScalingParameters,
																					evaluationTimeScalingParameters, interpolationVarianceScheme, evaluationTimeScalingScheme, true);
		ExtensionLIBORMarketModelWithBridge extensionLIBORMarketModelWithBridge = (ExtensionLIBORMarketModelWithBridge) liborMarketModel.getCloneWithModifiedCovarianceModel(newCovarianceModel);
		
		
		
		return extensionLIBORMarketModelWithBridge;
	}
	
//	@Override
//	public RandomVariableInterface getNumeraire(double time) throws CalculationException {
//		int liborTimeIndex = getLiborPeriodIndex(time);
//		int timeIndex	   = Math.max(getTimeIndex(time), - getTimeIndex(time) - 1);
//		
//		if(liborTimeIndex < 0) {
//			// Interpolation of Numeraire: log linear interpolation.
//			int upperIndex = -liborTimeIndex-1;
//			int lowerIndex = upperIndex-1;
//			if(lowerIndex < 0) throw new IllegalArgumentException("Numeraire requested for time " + time + ". Unsupported");
//
//			/*RandomVariableInterface numeraire = getNumeraire(getLiborPeriod(upperIndex)).div( 
//					getInterpolatedLibor(timeIndex , timeIndex).mult(getLiborPeriod(upperIndex) - time).add(1.0) );*/
//			/*
//			 * Adjust for discounting, i.e. funding or collateralization
//			 */
//			double upperTime    				= getLiborPeriod(upperIndex);   //T_i+1
//			
//			RandomVariableInterface numeraire   = getUnAdjustedNumeraire(upperTime).div(getLIBOR(time, time, upperTime).mult(upperTime - time).add(1.0));
//
//			
//			
//			RandomVariableInterface bridge		= getBrownianBridge(lowerIndex, time);
//			RandomVariableInterface evaluationTimeScalingFactor = ((AbstractLiborCovarianceModelWithInterpolation) getCovarianceModel()).getEvaluationTimeScalingFactor(timeIndex);
//			double lowerTime = getLiborPeriod(lowerIndex);
//			double alpha = (upperTime - time) / (upperTime - lowerTime);
//			
//			double analyticLiborShortPeriod				= getForwardRateCurve().getForward(getAnalyticModel(), time, upperTime -time);
//			double analyticLibor					 	= getForwardRateCurve().getForward(getAnalyticModel(), lowerTime, upperTime -lowerTime);
//			double analyticInterpolatedOnePlusLiborDt		= Math.exp(Math.log(1 + analyticLibor * (upperTime - lowerTime)) * alpha);
//			double analyticOnePlusLiborDt					= (1 + analyticLiborShortPeriod * (upperTime -time));
//			double adjustment = analyticOnePlusLiborDt / analyticInterpolatedOnePlusLiborDt;
//			RandomVariableInterface numeraire2 = (getUnAdjustedNumeraire(upperTime).log().mult(1 - alpha).add(getUnAdjustedNumeraire(lowerTime).log().mult(alpha))
//					.sub(bridge.mult(evaluationTimeScalingFactor))).exp();
//		    //numeraire = numeraire2;
//			
//			
//			if(getDiscountCurve() != null) {
//				// This includes a control for zero bonds
//				double deterministicNumeraireAdjustment = numeraire.invert().getAverage() / getDiscountCurve().getDiscountFactor(getAnalyticModel(), time);
//				numeraire = numeraire.mult(deterministicNumeraireAdjustment);
//				//System.out.println("inLMM numeraireAdj interpol " + time + "\t" + deterministicNumeraireAdjustment);
//				
//				/*
//				RandomVariableInterface bridge		= getBrownianBridge(lowerIndex, time);
//				RandomVariableInterface evaluationTimeScalingFactor = covarianceModel.getEvaluationTimeScalingFactor(timeIndex);
//				double lowerTime = getLiborPeriod(lowerIndex);
//				double alpha = (upperTime - time) / (upperTime - lowerTime);
//				
//				double analyticLiborShortPeriod				= getForwardRateCurve().getForward(getAnalyticModel(), time, upperTime -time);
//				double analyticLibor					 	= getForwardRateCurve().getForward(getAnalyticModel(), lowerTime, upperTime -lowerTime);
//				double analyticInterpolatedOnePlusLiborDt		= Math.exp(Math.log(1 + analyticLibor * (upperTime - lowerTime)) * alpha);
//				double analyticOnePlusLiborDt					= (1 + analyticLiborShortPeriod * (upperTime -time));
//				double adjustment = analyticOnePlusLiborDt / analyticInterpolatedOnePlusLiborDt;
//				RandomVariableInterface numeraire2 = (getUnAdjustedNumeraire(upperTime).log().mult(1 - alpha).add(getUnAdjustedNumeraire(lowerTime).log().mult(alpha))
//						.sub(bridge.mult(evaluationTimeScalingFactor))).exp().div(adjustment);
//				System.out.println(numeraire);
//				System.out.println(numeraire2);
//				*/
//			}
//
//			return numeraire;
//		}
//		return super.getNumeraire(time);
//	}
	
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
		RandomVariableInterface evaluationTimeScalingFactor = ((AbstractLiborCovarianceModelWithInterpolation) getCovarianceModel()).getEvaluationTimeScalingFactor(evaluationTimeIndex);
		
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
			
			BrownianBridgeWithVariance brownianBridge = new BrownianBridgeWithVariance(bridgeDiscretization, interpolationDriver, ((AbstractLiborCovarianceModelWithInterpolation) getCovarianceModel()).getVarianceForInterpolationPeriod(liborIndex),
					getRandomVariableFactory());
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
	
	/*
	@Override
	public Object clone() {
		try {
			Map<String, Object> properties = new HashMap<String, Object>();
			properties.put("measure",		getMeasure().name());
			properties.put("stateSpace",	getStateSpace().name());
			return new ExtensionLIBORMarketModelWithBridge(getLiborPeriodDiscretization(), getAnalyticModel(),
					getForwardRateCurve(), getDiscountCurve(), getRandomVariableFactory(), getCovarianceModel(),
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
	*/
	
	public RandomVariableInterface getUnAdjustedNumeraire(double time) throws CalculationException {
		int liborTimeIndex = getLiborPeriodIndex(time);
		if(liborTimeIndex < 0) {
			throw new UnsupportedOperationException("Unadjusted Numeraire only available for Tenortimes. Requested " + time + " not supported!" );
		}
		RandomVariableInterface unAdjustedNumeraire = getNumeraires().get(liborTimeIndex);
		if(unAdjustedNumeraire == null) {
			getNumeraire(time);
			unAdjustedNumeraire = getNumeraires().get(liborTimeIndex);
		}
		return unAdjustedNumeraire;
	}
	
}
