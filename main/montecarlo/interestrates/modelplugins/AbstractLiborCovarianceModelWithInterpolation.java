package montecarlo.interestrates.modelplugins;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORMarketModelInterface;
import net.finmath.montecarlo.interestrate.modelplugins.AbstractLIBORCovarianceModelParametric;
import net.finmath.montecarlo.interestrate.products.AbstractLIBORMonteCarloProduct;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretizationInterface;

public abstract class AbstractLiborCovarianceModelWithInterpolation extends AbstractLIBORCovarianceModelParametric {

	public AbstractLiborCovarianceModelWithInterpolation(TimeDiscretizationInterface timeDiscretization, TimeDiscretizationInterface liborPeriodDiscretization, int numberOfFactors) {
		super(timeDiscretization, liborPeriodDiscretization, numberOfFactors);
	}
	
	public abstract RandomVariableInterface getVarianceForInterpolation(double time);
	
	public RandomVariableInterface 			getVarianceForInterpolation(int timeIndex) {
		return getVarianceForInterpolation(getTimeDiscretization().getTime(timeIndex));
	}
	
	public RandomVariableInterface[] getVarianceForInterpolationPeriod(int liborIndex) {
		int		liborPeriodStartTimeIndex = getTimeDiscretization().getTimeIndex(getLiborPeriodDiscretization().getTime(liborIndex));
		if(liborPeriodStartTimeIndex < 0) liborPeriodStartTimeIndex = -liborPeriodStartTimeIndex - 1;
		int 	liborPeriodEndTimeIndex	  = getTimeDiscretization().getTimeIndex(getLiborPeriodDiscretization().getTime(liborIndex + 1));
		if(liborPeriodStartTimeIndex < 0) liborPeriodEndTimeIndex = -liborPeriodEndTimeIndex - 1;
		
		RandomVariableInterface[] varianceForInterpolationPeriod = new RandomVariableInterface[liborPeriodEndTimeIndex - liborPeriodStartTimeIndex];
		for (int timeIndex = 0; timeIndex < varianceForInterpolationPeriod.length; timeIndex++) {
			varianceForInterpolationPeriod[timeIndex] = getVarianceForInterpolation(liborPeriodStartTimeIndex + timeIndex);
		}
		return varianceForInterpolationPeriod;
	}

}
