package montecarlo.interestrate;

import java.util.Arrays;
import java.util.Map;

import montecarlo.interestrate.modelplugins.LIBORWithBrownianBridgeCovarianceModelInterface;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.AnalyticModelInterface;
import net.finmath.marketdata.model.curves.DiscountCurveInterface;
import net.finmath.marketdata.model.curves.ForwardCurveInterface;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORMarketModelInterface;
import net.finmath.montecarlo.interestrate.LIBORModelInterface;
import net.finmath.montecarlo.interestrate.TermStructureModelInterface;
import net.finmath.montecarlo.interestrate.modelplugins.AbstractLIBORCovarianceModel;
import net.finmath.montecarlo.interestrate.modelplugins.TermStructureCovarianceModelInterface;
import net.finmath.montecarlo.model.AbstractModel;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretizationInterface;

public class LIBORMarketModelWithBrownianBridge extends LIBORMarketModel implements TermStructureModelInterface {

	public LIBORMarketModelWithBrownianBridge(TimeDiscretizationInterface liborPeriodDiscretization,
			ForwardCurveInterface forwardRateCurve, DiscountCurveInterface discountCurve,
			AbstractLIBORCovarianceModel covarianceModel) throws CalculationException {
		super(liborPeriodDiscretization, forwardRateCurve, discountCurve, covarianceModel);
		
	}

	
	@Override
	public RandomVariableInterface[] getInitialState() {
		
		RandomVariableInterface[] initialLiborStates = super.getInitialState();

		RandomVariableInterface[] initialStates = new RandomVariableInterface[initialLiborStates.length + 1];
		for (int factorIndex = 0; factorIndex < initialLiborStates.length; factorIndex++) {
			initialStates[factorIndex] = initialLiborStates[factorIndex];
		}
		//testing
		initialStates[initialStates.length - 1] = getRandomVariableForConstant(0.0);

		return initialStates;
	}
	
	
	@Override
	public RandomVariableInterface[] getDrift(int timeIndex, RandomVariableInterface[] realizationAtTimeIndex, RandomVariableInterface[] realizationPredictor) {
		RandomVariableInterface[] liborDrifts = super.getDrift(timeIndex, realizationAtTimeIndex, realizationPredictor);
		//TODO: implement BB drift
		
		RandomVariableInterface[] drifts = new RandomVariableInterface[liborDrifts.length + 1];
		for (int factorIndex = 0; factorIndex < liborDrifts.length; factorIndex++) {
			drifts[factorIndex] = liborDrifts[factorIndex];
		}
		//testing
		drifts[drifts.length - 1] = getRandomVariableForConstant(0.0);

		return drifts;
	}
	
	@Override
	public RandomVariableInterface[] getFactorLoading(int timeIndex, int componentIndex,
			RandomVariableInterface[] realizationAtTimeIndex) {
		
		if(componentIndex<getNumberOfLibors()) {
			
			RandomVariableInterface[] liborFactorLoadings = super.getFactorLoading(timeIndex, componentIndex, realizationAtTimeIndex);
			RandomVariableInterface[] factorLoadings = new RandomVariableInterface[liborFactorLoadings.length + 1];
			for (int factorIndex = 0; factorIndex < liborFactorLoadings.length; factorIndex++) {
				factorLoadings[factorIndex] = liborFactorLoadings[factorIndex];
			}
			factorLoadings[factorLoadings.length - 1] = getRandomVariableForConstant(0.0); //no correlation to BB
			return factorLoadings;
		}
		return new RandomVariableInterface[] {getRandomVariableForConstant(0.0)};
	}

	@Override
	public RandomVariableInterface getLIBOR(int evaluationTimeIndex, int processTimeIndex) throws CalculationException {
		
		double processTime    = getTime(processTimeIndex);
		int    liborIndex = getLiborPeriodIndex(processTime);
		//processTimeIndex at LiborIndex
		if(liborIndex>=0) return getProcessValue(evaluationTimeIndex, liborIndex);
		
		//start interpolation
		int previousLiborIndex = (-liborIndex-1)-1;						 //i
		double previousLiborTime = getLiborPeriod(previousLiborIndex);   //T_i
		double nextLiborTime     = getLiborPeriod(previousLiborIndex+1); //T_i+1
		double periodLenght      = nextLiborTime - previousLiborIndex;
		return getLIBOR(evaluationTimeIndex, previousLiborIndex).mult((nextLiborTime-processTime)/periodLenght) //L^i(t)*((T_i+1-t)/(T_i+1-T_i))
				.add(getLIBOR(evaluationTimeIndex, previousLiborIndex + 1).mult((processTime-previousLiborTime)/periodLenght)) //L^i+1(t)*((t-T_i)/(T_i+1-T_i))
				.add(getBrownianBridge(processTimeIndex)/*.mult((nextLiborTime-processTime)/periodLenght)*/);
	}

	@Override
	public RandomVariableInterface getLIBOR(double time, double periodStart, double periodEnd)
			throws CalculationException {
		
			int periodStartIndex = getTimeIndex(periodStart);
			int periodEndIndex	= getTimeIndex(periodEnd);
	
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
				double subPeriodLength = getTime(periodIndex+1) - getTime(periodIndex); //including BB period lenght
				RandomVariableInterface liborOverSubPeriod = getLIBOR(timeIndex, periodIndex);

				accrualAccount = accrualAccount == null ? liborOverSubPeriod.mult(subPeriodLength).add(1.0) : accrualAccount.accrue(liborOverSubPeriod, subPeriodLength);
			}

			RandomVariableInterface libor = accrualAccount.sub(1.0).div(periodEnd - periodStart);

			return libor;
	}

	private RandomVariableInterface getBrownianBridge(int processTimeIndex) throws CalculationException {
		int brownianBridgeProcessIndex = getNumberOfLibors() - 1;
		return getProcessValue(processTimeIndex, brownianBridgeProcessIndex);
	}
	
	@Override
	public int getNumberOfComponents() {
		return getNumberOfLibors() + 1;  //BB is last component;
	}
	@Override
	public int getNumberOfLibors() {
		return getLiborPeriodDiscretization().getNumberOfTimeSteps();
	}
	
}

