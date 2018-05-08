package montecarlo.interestrate;

import java.util.Arrays;

import montecarlo.interestrate.modelplugins.AbstractLiborWithInterpolationCovarianceModel;
import montecarlo.interestrate.modelplugins.LiborWithInterpolationCovarianceModelFromArray;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.BrownianBridge;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORModelInterface;
import net.finmath.montecarlo.process.AbstractProcess;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationInterface;

public class LIBORModelWithBrownianBridge extends AbstractLIBORModelWithInterpolation {

	enum InterpolationScheme { LINEAR, LOGLINEAR }
	
	InterpolationScheme interpolationScheme = InterpolationScheme.LOGLINEAR;
	private BrownianBridge[]        brownianBridges;
	private AbstractLiborWithInterpolationCovarianceModel covarianceModel;
	
	
	public LIBORModelWithBrownianBridge(LIBORModelInterface liborModel, AbstractProcess liborProcess, int seed, double[] interpolationFactorArray) {
		super(liborModel, liborProcess);
		this.covarianceModel = new LiborWithInterpolationCovarianceModelFromArray(((LIBORMarketModel) liborModel).getCovarianceModel(), interpolationFactorArray);
		this.brownianBridges = new BrownianBridge[liborModel.getNumberOfLibors() - 1]; //last one not interpolated
		TimeDiscretizationInterface completeTimeDiscretization = liborModel.getTimeDiscretization();
		TimeDiscretizationInterface liborPeriodDiscretization  = liborModel.getLiborPeriodDiscretization();
		TimeDiscretizationInterface timeDiscretization;
		for (int liborIndex = 0; liborIndex < brownianBridges.length; liborIndex++) {
			double intialLiborTime = liborPeriodDiscretization.getTime(liborIndex);
			double lastLiborTime   = liborPeriodDiscretization.getTime(liborIndex + 1);
			int initialTimeIndex   = completeTimeDiscretization.getTimeIndex(intialLiborTime);
			int lastTimeIndex      = completeTimeDiscretization.getTimeIndex(lastLiborTime);
			timeDiscretization     = new TimeDiscretization(Arrays.copyOfRange(completeTimeDiscretization.getAsDoubleArray(), initialTimeIndex, lastTimeIndex + 1));
			brownianBridges[liborIndex] = new BrownianBridge(timeDiscretization, liborModel.getProcess().getNumberOfPaths(), seed, liborModel.getRandomVariableForConstant(0.0), liborModel.getRandomVariableForConstant(0.0));
			seed = seed + 1;
		}
	}

	@Override
	public RandomVariableInterface getInterpolatedValue(int evaluationTimeIndex, int processTimeIndex) throws CalculationException {
		double processTime    = getTime(processTimeIndex);
		int    previousLiborIndex = getLiborPeriodIndex(processTime); 
		previousLiborIndex = (-previousLiborIndex-1)-1;	//i
		if(previousLiborIndex == getLiborPeriodDiscretization().getNumberOfTimeSteps()) {
			throw new IndexOutOfBoundsException("There is no interpolation inside the last LIBOR period!");
		}
		double previousLiborTime = getLiborPeriod(previousLiborIndex);   //T_i
		double nextLiborTime     = getLiborPeriod(previousLiborIndex+1); //T_i+1
		double periodLenght      = nextLiborTime - previousLiborTime;
		
		if(interpolationScheme == InterpolationScheme.LINEAR) {
			
			return (getLIBORWithInterpolation(evaluationTimeIndex, getTimeIndex(previousLiborTime)).mult((nextLiborTime-processTime)/periodLenght)) //L^i(t)*((T_i+1-t)/(T_i+1-T_i))
					.add(getLIBORWithInterpolation(evaluationTimeIndex, getTimeIndex(nextLiborTime)).mult((processTime-previousLiborTime)/periodLenght)) //L^i+1(t)*((t-T_i)/(T_i+1-T_i))
					.add(getBrownianBridge(previousLiborIndex, processTime))/*.mult((nextLiborTime-processTime)/periodLenght)*/;
		}
		
		if(interpolationScheme == InterpolationScheme.LOGLINEAR) {
			//not meaningfull
			return  (getLIBORWithInterpolation(evaluationTimeIndex, getTimeIndex(previousLiborTime)).log().mult((nextLiborTime-processTime)/periodLenght)) //L^i(t)*((T_i+1-t)/(T_i+1-T_i))
					.add(getLIBORWithInterpolation(evaluationTimeIndex, getTimeIndex(nextLiborTime)).log().mult((processTime-previousLiborTime)/periodLenght)).exp() //L^i+1(t)*((t-T_i)/(T_i+1-T_i))
					.add(getBrownianBridge(previousLiborIndex, processTime))/*.mult((nextLiborTime-processTime)/periodLenght)*/;
		}
		throw new UnsupportedOperationException();
	}

	private RandomVariableInterface getBrownianBridge(int liborIndex, double time) {
		RandomVariableInterface brownianBridgeValue = getRandomVariableForConstant(0.0);
		TimeDiscretizationInterface brownianTimeDiscretizazion = brownianBridges[liborIndex].getTimeDiscretization();
		int brownianIndex  = brownianTimeDiscretizazion.getTimeIndex(time);
		int globalTimeIndex = brownianTimeDiscretizazion.getTimeIndex(getLiborPeriod(liborIndex)) + brownianIndex;
		for (int timeIndex = 0; timeIndex < brownianIndex; timeIndex++) {
			brownianBridgeValue = brownianBridges[liborIndex].getIncrement(timeIndex, 0).mult(0.0)/*.mult(covarianceModel.getFactorLoadingForInterpolation(globalTimeIndex)[0])*/.add(brownianBridgeValue);
		}
		return brownianBridgeValue;
	}	
	
}
