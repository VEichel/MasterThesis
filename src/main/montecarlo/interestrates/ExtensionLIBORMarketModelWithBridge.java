package montecarlo.interestrates;

import montecarlo.interestrates.modelplugins.AbstractLiborCovarianceModelWithInterpolation;
import montecarlo.interestrates.modelplugins.LiborCovarianceModelWithInterpolation;
import montecarlo.interestrates.modelplugins.LiborCovarianceModelWithInterpolation.EvaluationTimeScalingScheme;
import montecarlo.interestrates.modelplugins.LiborCovarianceModelWithInterpolation.InterpolationVarianceScheme;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.AnalyticModelInterface;
import net.finmath.marketdata.model.curves.DiscountCurveInterface;
import net.finmath.marketdata.model.curves.ForwardCurveInterface;
import net.finmath.montecarlo.AbstractRandomVariableFactory;
import net.finmath.montecarlo.BrownianMotionInterface;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.modelplugins.AbstractLIBORCovarianceModel;
import net.finmath.montecarlo.process.AbstractProcess;
import net.finmath.montecarlo.process.ProcessEulerScheme;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationInterface;

import java.util.Arrays;
import java.util.Map;


public class ExtensionLIBORMarketModelWithBridge extends LIBORMarketModel{

	public enum InterpolationScheme { LINEAR, LOGLINEAR }
	private InterpolationScheme interpolationScheme = InterpolationScheme.LOGLINEAR;

	private RandomVariableInterface[][] brownianBridgeValues;

	private final BrownianMotionInterface interpolationDriver;


	public ExtensionLIBORMarketModelWithBridge(TimeDiscretizationInterface liborPeriodDiscretization,
											   AnalyticModelInterface analyticModel,
											   ForwardCurveInterface forwardRateCurve,
											   DiscountCurveInterface discountCurve,
                                               AbstractRandomVariableFactory randomVariableFactory,
                                               AbstractLIBORCovarianceModel covarianceModel,
                                               CalibrationItem[] calibrationItems,
                                               Map<String, ?> properties,
                                               BrownianMotionInterface interpolationDriver) throws CalculationException {
		super(liborPeriodDiscretization, analyticModel, forwardRateCurve,
                discountCurve, randomVariableFactory, covarianceModel, calibrationItems, properties);
		this.interpolationDriver = interpolationDriver;
        if(((ProcessEulerScheme) interpolationDriver).getScheme() == ProcessEulerScheme.Scheme.PREDICTOR_CORRECTOR)
            throw new IllegalArgumentException("Interpolation Process is not allowed to have Predictor Corrector Scheme!");
	}

	public static ExtensionLIBORMarketModelWithBridge getLIBORMarketModelExtendedWithBridge(LIBORMarketModel liborMarketModel,
																							double[] interpolationParameters, double[] evaluationTimeScalingParameters,
																							InterpolationVarianceScheme interpolationVarianceScheme, EvaluationTimeScalingScheme evaluationTimeScalingScheme) {
		
		AbstractLiborCovarianceModelWithInterpolation newCovarianceModel = new LiborCovarianceModelWithInterpolation(liborMarketModel.getCovarianceModel(), evaluationTimeScalingParameters,
																					evaluationTimeScalingParameters, interpolationVarianceScheme, evaluationTimeScalingScheme, true);
		ExtensionLIBORMarketModelWithBridge extensionLIBORMarketModelWithBridge = (ExtensionLIBORMarketModelWithBridge) liborMarketModel.getCloneWithModifiedCovarianceModel(newCovarianceModel);
		
		
		
		return extensionLIBORMarketModelWithBridge;
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
		//RandomVariableInterface evaluationTimeScalingFactor = ((AbstractLiborCovarianceModelWithInterpolation) getCovarianceModel()).getEvaluationTimeScalingFactor(evaluationTimeIndex);
		
		if(previousLiborIndex<0)	previousLiborIndex = (-previousLiborIndex-1)-1;	//i
		double previousLiborTime 			= getLiborPeriod(previousLiborIndex);   //T_i
		double nextLiborTime     			= getLiborPeriod(previousLiborIndex+1); //T_i+1
		double periodLenght     			= nextLiborTime - previousLiborTime;
		double shortPeriodLenght			= nextLiborTime	- processTime;
		double alpha            			= (nextLiborTime - processTime) / periodLenght;
		RandomVariableInterface startLibor	= getLIBOR(evaluationTimeIndex, previousLiborIndex); //L(T_i,T_i+1)
		RandomVariableInterface bridge		= getBrownianBridgeValue(previousLiborIndex, processTime);
		RandomVariableInterface libor;
		/*if(interpolationScheme == InterpolationScheme.LINEAR) {
			libor = startLibor.mult(periodLenght).add(1.0).mult(alpha).add(1-alpha).mult(bridge.mult(evaluationTimeScalingFactor)).sub(1.0).div(shortPeriodLenght);
		}*/
		
		if(interpolationScheme == InterpolationScheme.LOGLINEAR) {
			libor = startLibor.mult(periodLenght).add(1.0).log().mult(alpha).add(bridge).exp().sub(1.0).div(shortPeriodLenght);
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
	public RandomVariableInterface getBrownianBridgeValue(int liborIndex, double time) throws CalculationException {
		
		if(brownianBridgeValues == null) {
			brownianBridgeValues = new RandomVariableInterface[getNumberOfLibors()][];
		}
		//Lazy Init for each Brownian Bridge
		if(brownianBridgeValues[liborIndex] == null) {

			double liborPeriodStart = getLiborPeriod(liborIndex);
			double liborPeriodEnd   = getLiborPeriod(liborIndex + 1);
			TimeDiscretizationInterface tempBridgeDiscretization    = new TimeDiscretization(
					Arrays.copyOfRange(getTimeDiscretization().getAsDoubleArray(), getTimeIndex(liborPeriodStart), getTimeIndex(liborPeriodEnd) + 1));
            brownianBridgeValues[liborIndex] = new RandomVariableInterface[tempBridgeDiscretization.getNumberOfTimeSteps()];

            brownianBridgeValues[liborIndex][0] = getRandomVariableForConstant(0.0);
            for(int bridgeTimeIndex = 1; bridgeTimeIndex<tempBridgeDiscretization.getNumberOfTimeSteps(); bridgeTimeIndex++) {
                double currentTime	= tempBridgeDiscretization.getTime(bridgeTimeIndex);
                double nextTime		= tempBridgeDiscretization.getTime(bridgeTimeIndex+1);
                double alpha		= (nextTime-currentTime)/(liborPeriodEnd-currentTime);

                int    generatorTimeIndex        = interpolationDriver.getTimeDiscretization().getTimeIndex(currentTime);
                RandomVariableInterface variance = ((AbstractLiborCovarianceModelWithInterpolation)getCovarianceModel()).getVarianceForInterpolation(currentTime);

                // Calculate the next point using the "scheme" of the Brownian bridge
                brownianBridgeValues[liborIndex][bridgeTimeIndex] = brownianBridgeValues[liborIndex][bridgeTimeIndex-1].mult(1.0-alpha)
                        .add(interpolationDriver.getBrownianIncrement(generatorTimeIndex, 0)
                        .mult(variance.sqrt()).mult(Math.sqrt(1.0-alpha)));
            }
		}
		int bridgeTimeIndex = getTimeIndex(time) - getTimeIndex(getLiborPeriod(liborIndex));
		return brownianBridgeValues[liborIndex][bridgeTimeIndex];
	}

    public RandomVariableInterface getBrownianBridgeAnalyticCovariance(double time1, double time2) {
	    int liborPeriodEndIndex = getLiborPeriodIndex(time1);
        {
            int templiborPeriodEndIndex2 = getLiborPeriodIndex(time2);
            if (liborPeriodEndIndex > 0 || templiborPeriodEndIndex2 > 0 || liborPeriodEndIndex != templiborPeriodEndIndex2)
            	return getRandomVariableForConstant(0.0);
        }
        	liborPeriodEndIndex   = - liborPeriodEndIndex - 1;
		double liborPeriodEndTime = getLiborPeriod(liborPeriodEndIndex);
		int liborPeriodStartIndex = liborPeriodEndIndex - 1;
		int lowerIndex 			  = Math.min(getTimeIndex(time1), getTimeIndex(time2));
		int startTimeIndex		  = getTimeIndex(getLiborPeriod(liborPeriodStartIndex));
		RandomVariableInterface covariance 		  = getRandomVariableForConstant(
													(liborPeriodEndTime - time1) * (liborPeriodEndTime - time2) );


		for (int j = startTimeIndex; j < lowerIndex; j++) {
			RandomVariableInterface variance = ((AbstractLiborCovarianceModelWithInterpolation)getCovarianceModel()).getVarianceForInterpolation(j);
			double currentTime = getTime(j);
			double nextTime	   = getTime(j+1);
			covariance = covariance.add(
					variance.mult((nextTime - currentTime)/((liborPeriodEndTime - nextTime)*(liborPeriodEndTime - currentTime))) );
		}
		return covariance;
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
}
