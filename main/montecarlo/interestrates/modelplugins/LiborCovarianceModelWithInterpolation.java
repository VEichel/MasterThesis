package montecarlo.interestrates.modelplugins;

import java.util.Arrays;

import org.apache.commons.lang3.ArrayUtils;

import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.interestrate.modelplugins.AbstractLIBORCovarianceModel;
import net.finmath.montecarlo.interestrate.modelplugins.AbstractLIBORCovarianceModelParametric;
import net.finmath.stochastic.RandomVariableInterface;
import sun.reflect.generics.tree.ArrayTypeSignature;

public class LiborCovarianceModelWithInterpolation extends AbstractLiborCovarianceModelWithInterpolation {

	public enum InterpolationVarianceScheme { PIECEWISECONSTANT, CONSTANT, FINEST }
	public enum EvaluationTimeScalingScheme	{ PIECEWISECONSTANT, CONSTANT, FINEST }
		
	
	private final AbstractLIBORCovarianceModel	nonInterpolationModel;
	private final InterpolationVarianceScheme interpolationVarianceScheme;
	private final EvaluationTimeScalingScheme       evaluationTimeScalingScheme;
	private final double[] interpolationParameters;
	private final double[] evaluationTimeScalingParameters;
	
	private final boolean nonInterpolationModelIsCalibrated;
	
	/**Standard Constructor
	 * 
	 * @param nonInterpolationModel
	 * @param interpolationParameters
	 * @param interpolationVarianceScheme
	 */
	public LiborCovarianceModelWithInterpolation(
			AbstractLIBORCovarianceModel nonInterpolationModel,
			double[] interpolationParameters,
			double[] evaluationTimeScalingParameters,
			InterpolationVarianceScheme interpolationVarianceScheme,
			EvaluationTimeScalingScheme	evaluationTimeScalingScheme,
			boolean 					nonInterpolationModelIsCalibrated) {
		super(nonInterpolationModel.getTimeDiscretization(), nonInterpolationModel.getLiborPeriodDiscretization(), nonInterpolationModel.getNumberOfFactors());
		
		this.nonInterpolationModel = nonInterpolationModel;
		this.interpolationVarianceScheme = interpolationVarianceScheme;
		this.evaluationTimeScalingScheme = evaluationTimeScalingScheme;
		this.nonInterpolationModelIsCalibrated = nonInterpolationModelIsCalibrated;
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
		
		switch (evaluationTimeScalingScheme) {
				
		case PIECEWISECONSTANT:
			if(evaluationTimeScalingParameters.length != getLiborPeriodDiscretization().getNumberOfTimeSteps()) {
				throw new ArrayIndexOutOfBoundsException("For piecewise constant evaluation time scaling model the parameters lenght must match the number of LIBORs!");
			}
			break;

		case CONSTANT:
			if(evaluationTimeScalingParameters.length != 1) {
				throw new ArrayIndexOutOfBoundsException("For constant evaluation time scaling model the parameters lenght must be 1!");
			}
			break;
			
		case FINEST:
			if(evaluationTimeScalingParameters.length != getTimeDiscretization().getNumberOfTimeSteps()) {
				throw new ArrayIndexOutOfBoundsException("For finest evaluation time scaling model the parameters lenght must be the number of time steps!");
			}
			break;
			
		default:
			throw new IllegalArgumentException("Evaluation time scaling model not known!");	
		}
		
		this.evaluationTimeScalingParameters = evaluationTimeScalingParameters;
	}
	
	/**Constructor using the Interpolation Variance Scheme: constant.
	 * 
	 * @param nonInterpolationModel
	 * @param interpolationParameters
	 */
	public LiborCovarianceModelWithInterpolation(
			AbstractLIBORCovarianceModelParametric nonInterpolationModel,
			double[] interpolationParameters) {
		this(nonInterpolationModel, interpolationParameters, new double[] {1.0}, InterpolationVarianceScheme.CONSTANT, EvaluationTimeScalingScheme.CONSTANT, false);
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
	public RandomVariableInterface getEvaluationTimeScalingFactor(double time) {
		
		switch (interpolationVarianceScheme) {
				
		case CONSTANT:
			return new RandomVariable(evaluationTimeScalingParameters[0]);

		case PIECEWISECONSTANT:
			int liborPeriod = getLiborPeriodDiscretization().getTimeIndex(time);
			if(liborPeriod<0) liborPeriod = (-liborPeriod - 1) - 1;
			
			return new RandomVariable(evaluationTimeScalingParameters[liborPeriod]);
			
		case FINEST:
			int timeIndex = getTimeDiscretization().getTimeIndex(time);
			if(timeIndex<0) timeIndex = -timeIndex - 1;
			
			return new RandomVariable(evaluationTimeScalingParameters[timeIndex]);
			
		default:
			break;
		}
		
		throw new IllegalArgumentException("Interpolation Variance Scheme not known!");	
	}
	
	@Override
	public double[] getParameter() {
		
		double[] temp;
		
		if(!nonInterpolationModelIsCalibrated) {
			double[] nonInterpolationModelParameters = ((AbstractLIBORCovarianceModelParametric) nonInterpolationModel).getParameter();
			temp = ArrayUtils.addAll(nonInterpolationModelParameters, interpolationParameters);
			return ArrayUtils.addAll(temp, evaluationTimeScalingParameters);
		} else {
		//If the covariance Model for the coarse Libor Tenors is already calibrated, we do not return those parameter.
			temp = interpolationParameters;
		}
		return ArrayUtils.addAll(temp, evaluationTimeScalingParameters);
	}

	@Override
	public Object clone() {
		return new LiborCovarianceModelWithInterpolation(nonInterpolationModel, interpolationParameters, evaluationTimeScalingParameters, interpolationVarianceScheme, evaluationTimeScalingScheme, nonInterpolationModelIsCalibrated);
	}

	@Override
	public AbstractLIBORCovarianceModelParametric getCloneWithModifiedParameters(double[] parameters) {
		
		int numberOfInterpolationParameters	= interpolationParameters.length;
		
		double[] newInterpolationParameters;
		double[] newEvaluationTimeScalingParameters;
		AbstractLIBORCovarianceModelParametric newNonInterpolationModel;
		
		if(!nonInterpolationModelIsCalibrated) {
			int numberOfNonInterpolationParameters		= ((AbstractLIBORCovarianceModelParametric) nonInterpolationModel).getParameter().length;
			double[] newNonInterpolationParameters		= Arrays.copyOfRange(parameters, 0, numberOfNonInterpolationParameters);
			newInterpolationParameters    				= Arrays.copyOfRange(parameters, numberOfNonInterpolationParameters, numberOfNonInterpolationParameters + numberOfInterpolationParameters);
			newEvaluationTimeScalingParameters 			= Arrays.copyOfRange(parameters, numberOfNonInterpolationParameters + numberOfInterpolationParameters, parameters.length);
			newNonInterpolationModel 					= ((AbstractLIBORCovarianceModelParametric) nonInterpolationModel).getCloneWithModifiedParameters(newNonInterpolationParameters);
		} else {
			newInterpolationParameters    	   = Arrays.copyOfRange(parameters, 0, numberOfInterpolationParameters);
			newEvaluationTimeScalingParameters = Arrays.copyOfRange(parameters, numberOfInterpolationParameters, parameters.length);	
			newNonInterpolationModel		   = (AbstractLIBORCovarianceModelParametric) nonInterpolationModel;
		}
		
		return new LiborCovarianceModelWithInterpolation(newNonInterpolationModel, newInterpolationParameters, newEvaluationTimeScalingParameters, interpolationVarianceScheme, evaluationTimeScalingScheme, nonInterpolationModelIsCalibrated);
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
