package montecarlo.interestrate;

import java.util.Arrays;

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

public class LIBORMarketModelWithBrownianBridgeSimulation extends LIBORModelMonteCarloSimulation{

	private BrownianBridge[]        brownianBridges;

	public LIBORMarketModelWithBrownianBridgeSimulation(LIBORModelInterface model, AbstractProcess process, int seed) {
		super(model, process);
		this.brownianBridges = new BrownianBridge[model.getNumberOfLibors()];
		TimeDiscretizationInterface completeTimeDiscretization = model.getTimeDiscretization();
		TimeDiscretizationInterface liborPeriodDiscretization  = model.getLiborPeriodDiscretization();
		TimeDiscretizationInterface timeDiscretization;
		for (int liborIndex = 0; liborIndex < brownianBridges.length; liborIndex++) {
			double intialLiborTime = liborPeriodDiscretization.getTime(liborIndex);
			double lastLiborTime   = liborPeriodDiscretization.getTime(liborIndex + 1);
			int initialTimeIndex   = completeTimeDiscretization.getTimeIndex(intialLiborTime);
			int lastTimeIndex      = completeTimeDiscretization.getTimeIndex(lastLiborTime);
			timeDiscretization     = new TimeDiscretization(Arrays.copyOfRange(completeTimeDiscretization.getAsDoubleArray(), initialTimeIndex, lastTimeIndex + 1));
			brownianBridges[liborIndex] = new BrownianBridge(timeDiscretization, process.getNumberOfPaths(), seed, model.getRandomVariableForConstant(0.0), model.getRandomVariableForConstant(0.0));
			seed = seed + 1;
		}
		
	}

	public RandomVariableInterface getBrownianBridge(int timeIndex) {
		interpolationBrownianMotion.
	}
}
