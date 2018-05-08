package montecarlo.interestrate;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import montecarlo.interestrate.LIBORModelWithBrownianBridge.InterpolationScheme;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.AnalyticModelInterface;
import net.finmath.marketdata.model.curves.DiscountCurveFromForwardCurve;
import net.finmath.marketdata.model.curves.DiscountCurveInterface;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurveInterface;
import net.finmath.montecarlo.BrownianBridge;
import net.finmath.montecarlo.BrownianMotionInterface;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORMarketModelInterface;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulation;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationInterface;
import net.finmath.montecarlo.interestrate.LIBORMarketModel.CalibrationItem;
import net.finmath.montecarlo.interestrate.modelplugins.AbstractLIBORCovarianceModel;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORCorrelationModelExponentialDecay;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORCovarianceModelFromVolatilityAndCorrelation;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORVolatilityModel;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORVolatilityModelFourParameterExponentialForm;
import net.finmath.montecarlo.process.ProcessEulerScheme;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationInterface;

public class LiborMarketModelWithBridgeInterpolation extends LIBORMarketModel {

	enum InterpolationScheme { LINEAR, LOGLINEAR }
	InterpolationScheme interpolationScheme = InterpolationScheme.LOGLINEAR;
	
	private RandomVariableInterface[][] brownianBridgeValues;

	private int seed;

	public LiborMarketModelWithBridgeInterpolation(TimeDiscretizationInterface liborPeriodDiscretization,
			AnalyticModelInterface analyticModel, ForwardCurveInterface forwardRateCurve,
			DiscountCurveInterface discountCurve, AbstractLIBORCovarianceModel covarianceModel,
			CalibrationItem[] calibrationItems, Map<String, ?> properties, int seed) throws CalculationException {
		super(liborPeriodDiscretization, analyticModel, forwardRateCurve, discountCurve, covarianceModel, calibrationItems,
				properties);
		
		this.seed = seed;
	}
	
	public LiborMarketModelWithBridgeInterpolation(
			TimeDiscretizationInterface			liborPeriodDiscretization,
			AnalyticModelInterface				analyticModel,
			ForwardCurveInterface				forwardRateCurve,
			DiscountCurveInterface				discountCurve,
			AbstractLIBORCovarianceModel		covarianceModel,
			CalibrationItem[]					calibrationItems,
			Map<String, ?>						properties
			) throws CalculationException {
		this(liborPeriodDiscretization, analyticModel, forwardRateCurve, discountCurve, covarianceModel, calibrationItems, properties, /*seed*/ 231317);
	}
	
	/**
	 * 
	 * @param evaluationTimeIndex 
	 * @param processTimeIndex    t
	 * @return L(t,m(t)+1;evaluationTimeIndex), where m(t) is the smallest LiborIndex greater than t.
	 * @throws CalculationException
	 */
	public RandomVariableInterface getInterpolatedValue(int evaluationTimeIndex, int processTimeIndex) throws CalculationException {
		double processTime        = getTime(processTimeIndex);
		int    previousLiborIndex = getLiborPeriodIndex(processTime); 
		
		//if(previousLiborIndex>=0) throw new UnsupportedOperationException("This method is only for inner period LIBORs!");
		
		if(previousLiborIndex<0)previousLiborIndex 					= (-previousLiborIndex-1)-1;	//i
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
		
		double analyticStartLibor				= getForwardRateCurve().getForward(getAnalyticModel(), previousLiborTime, nextLiborTime-previousLiborTime);
		double adjustment 						= (Math.exp(Math.log(1 + analyticStartLibor * periodLenght) * alpha) - 1) / shortPeriodLenght;
		libor = libor.mult(periodLenght).add(1.0).div(adjustment).sub(1.0).div(periodLenght);
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
			
			BrownianBridgeWithVariance brownianBridge = new BrownianBridgeWithVariance(bridgeDiscretization, getProcess().getNumberOfPaths(), seed, getRandomVariableForConstant(0.0), getRandomVariableForConstant(0.0));
			seed = seed + 1; //for future Bridges.
			brownianBridgeValues = new RandomVariableInterface[getNumberOfLibors()][bridgeDiscretization.getNumberOfTimes()];
			
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
}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
