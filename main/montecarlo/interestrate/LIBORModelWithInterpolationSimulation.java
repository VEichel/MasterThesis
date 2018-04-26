package montecarlo.interestrate;

import java.util.Arrays;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.BrownianBridge;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.BrownianMotionInterface;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORMarketModelInterface;
import net.finmath.montecarlo.interestrate.LIBORModelInterface;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulation;
import net.finmath.montecarlo.process.AbstractProcess;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationInterface;

public class LIBORModelWithInterpolationSimulation extends LIBORModelMonteCarloSimulation{

	private AbstractLIBORModelWithInterpolation interpolatedModel;

	public LIBORModelWithInterpolationSimulation(AbstractLIBORModelWithInterpolation model, AbstractProcess process) {
		super(model.getLiborModel(), process);
		this.interpolatedModel = model;
	}

	public RandomVariableInterface getLIBORWithInterpolation(int evaluationTimeIndex, int processTimeIndex) throws CalculationException {
		return interpolatedModel.getLIBORWithInterpolation(evaluationTimeIndex, processTimeIndex);
	}
	
	@Override
	public RandomVariableInterface getNumeraire(double time) throws CalculationException {
		return interpolatedModel.getNumeraire(time);
	}
	
	@Override
	public RandomVariableInterface getLIBOR(double time, double periodStart, double periodEnd) throws CalculationException {
		return interpolatedModel.getLIBOR(time, periodStart, periodEnd);
	}

}
