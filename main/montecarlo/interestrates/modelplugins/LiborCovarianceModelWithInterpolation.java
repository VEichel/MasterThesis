package montecarlo.interestrates.modelplugins;

import java.util.Arrays;

import org.apache.commons.lang3.ArrayUtils;

import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.interestrate.modelplugins.AbstractLIBORCovarianceModelParametric;
import net.finmath.stochastic.RandomVariableInterface;

public class LiborCovarianceModelWithInterpolation extends AbstractLiborCovarianceModelWithInterpolation {

	public enum InterpolationVarianceScheme { PIECEWISECONSTANT, CONSTANT, FINEST }
		
	
	private final AbstractLIBORCovarianceModelParametric nonInterpolationModel;
	private final InterpolationVarianceScheme interpolationVarianceScheme;
	private final double[] interpolationParameters;
	
	/**Standard Constructor
	 * 
	 * @param nonInterpolationModel
	 * @param interpolationParameters
	 * @param interpolationVarianceScheme
	 */
	public LiborCovarianceModelWithInterpolation(
			AbstractLIBORCovarianceModelParametric nonInterpolationModel,
			double[] interpolationParameters,
			InterpolationVarianceScheme interpolationVarianceScheme) {
		super(nonInterpolationModel.getTimeDiscretization(), nonInterpolationModel.getLiborPeriodDiscretization(), nonInterpolationModel.getNumberOfFactors());
		
		this.nonInterpolationModel = nonInterpolationModel;
		this.interpolationVarianceScheme = interpolationVarianceScheme;
		
		//negative variance not allowed:
		for (int index = 0; index < interpolationParameters.length; index++) {
			if(interpolationParameters[index] < 0) interpolationParameters[index] = 0;
		}
		
		
		switch (interpolationVarianceScheme) {
		
		case PIECEWISECONSTANT:
			if(interpolationParameters.length != getLiborPeriodDiscretization().getNumberOfTimeSteps()) {
				throw new ArrayIndexOutOfBoundsException("For piecewise constant interpolation model the parameters lenght must match the number of LIBORs!");
			}
			break;

		case CONSTANT:
			if(interpolationParameters.length != 1) {
				throw new ArrayIndexOutOfBoundsException("For constant interpolation model the parameters lenght must be 1!");
			}
			break;
			
		case FINEST:
			if(interpolationParameters.length != getTimeDiscretization().getNumberOfTimeSteps()) {
				throw new ArrayIndexOutOfBoundsException("For finest interpolation model the parameters lenght must be the number of time steps!");
			}
			break;
			
		default:
			throw new IllegalArgumentException("Interpolation Variance Scheme not known!");	
		}
		
		this.interpolationParameters = interpolationParameters;
	}
	
	/**Constructor using the Interpolation Variance Scheme: constant.
	 * 
	 * @param nonInterpolationModel
	 * @param interpolationParameters
	 */
	public LiborCovarianceModelWithInterpolation(
			AbstractLIBORCovarianceModelParametric nonInterpolationModel,
			double[] interpolationParameters) {
		this(nonInterpolationModel, interpolationParameters, InterpolationVarianceScheme.CONSTANT);
	}

	
	@Override
	public RandomVariableInterface getVarianceForInterpolation(double time) {
		
		switch (interpolationVarianceScheme) {
		
		case CONSTANT:
			return new RandomVariable(interpolationParameters[0]);

		case PIECEWISECONSTANT:
			int liborPeriod = getLiborPeriodDiscretization().getTimeIndex(time);
			if(liborPeriod<0) liborPeriod = (-liborPeriod - 1) - 1;
			
			return new RandomVariable(interpolationParameters[liborPeriod]);
			
		case FINEST:
			int timeIndex = getTimeDiscretization().getTimeIndex(time);
			if(timeIndex<0) timeIndex = -timeIndex - 1;
			
			return new RandomVariable(interpolationParameters[timeIndex]);
			
		default:
			break;
		}
		throw new IllegalArgumentException("Interpolation Variance Scheme not known!");	
	}

	
	@Override
	public RandomVariableInterface[] getVarianceForInterpolationPeriod(int liborIndex) {
		int		liborPeriodStartTimeIndex = getTimeDiscretization().getTimeIndex(getLiborPeriodDiscretization().getTime(liborIndex));
		if(liborPeriodStartTimeIndex < 0) liborPeriodStartTimeIndex = -liborPeriodStartTimeIndex - 1;
		int 	liborPeriodEndTimeIndex	  = getTimeDiscretization().getTimeIndex(getLiborPeriodDiscretization().getTime(liborIndex + 1));
		if(liborPeriodStartTimeIndex < 0) liborPeriodEndTimeIndex = -liborPeriodEndTimeIndex - 1;
		
		
		RandomVariableInterface[] varianceForInterpolationPeriod;
		
		switch (interpolationVarianceScheme) {
		
		case CONSTANT:
			varianceForInterpolationPeriod = new RandomVariableInterface[liborPeriodEndTimeIndex - liborPeriodStartTimeIndex];
			for (int timeIndex = 0; timeIndex < varianceForInterpolationPeriod.length; timeIndex++) {
				varianceForInterpolationPeriod[timeIndex] = new RandomVariable(interpolationParameters[0]);
			}
			break;

		case PIECEWISECONSTANT:
			
			varianceForInterpolationPeriod = new RandomVariableInterface[liborPeriodEndTimeIndex - liborPeriodStartTimeIndex];
			for (int timeIndex = 0; timeIndex < varianceForInterpolationPeriod.length; timeIndex++) {
				varianceForInterpolationPeriod[timeIndex] = new RandomVariable(interpolationParameters[liborIndex]);
			}
			break;
			
		case FINEST:
			
			varianceForInterpolationPeriod = new RandomVariableInterface[liborPeriodEndTimeIndex - liborPeriodStartTimeIndex];
			
			int startTimeIndex = getTimeDiscretization().getTimeIndex(getLiborPeriodDiscretization().getTime(liborIndex));
			//if(startTimeIndex<0) startTimeIndex = -startTimeIndex - 1; can not happen
			
			for (int timeIndex = 0; timeIndex < varianceForInterpolationPeriod.length; timeIndex++) {
				varianceForInterpolationPeriod[timeIndex] = new RandomVariable(interpolationParameters[startTimeIndex + timeIndex]);
			}
			break;
			
		default:
			//should never happen
			throw new IllegalArgumentException("Interpolation Variance Scheme not known!");	
		}
		
		return varianceForInterpolationPeriod;
	}
	
	
	@Override
	public double[] getParameter() {
		
		double[] nonInterpolationModelParameters = nonInterpolationModel.getParameter();
		return ArrayUtils.addAll(nonInterpolationModelParameters, interpolationParameters);
	}

	@Override
	public Object clone() {
		return new LiborCovarianceModelWithInterpolation(nonInterpolationModel, interpolationParameters, interpolationVarianceScheme);
	}

	@Override
	public AbstractLIBORCovarianceModelParametric getCloneWithModifiedParameters(double[] parameters) {
		int numberOfNonInterpolationParameters = nonInterpolationModel.getParameter().length;
		double[] newNonInterpolationParameters = Arrays.copyOfRange(parameters, 0, numberOfNonInterpolationParameters);
		double[] newInterpolationParameters    = Arrays.copyOfRange(parameters, numberOfNonInterpolationParameters, parameters.length);
		AbstractLIBORCovarianceModelParametric newNonInterpolationModel = nonInterpolationModel.getCloneWithModifiedParameters(newNonInterpolationParameters);
		return new LiborCovarianceModelWithInterpolation(newNonInterpolationModel, newInterpolationParameters, interpolationVarianceScheme);
	}

	@Override
	public RandomVariableInterface[] getFactorLoading(int timeIndex, int component,
			RandomVariableInterface[] realizationAtTimeIndex) {
		return nonInterpolationModel.getFactorLoading(timeIndex, component, realizationAtTimeIndex);
	}

	@Override
	public RandomVariableInterface getFactorLoadingPseudoInverse(int timeIndex, int component, int factor,
			RandomVariableInterface[] realizationAtTimeIndex) {
		return nonInterpolationModel.getFactorLoadingPseudoInverse(timeIndex, component, factor, realizationAtTimeIndex);
	}

}
