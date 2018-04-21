package montecarlo.interestrate.modelplugins;

import net.finmath.montecarlo.interestrate.modelplugins.TermStructureCovarianceModelInterface;
import net.finmath.stochastic.RandomVariableInterface;

public interface LIBORWithBrownianBridgeCovarianceModelInterface extends TermStructureCovarianceModelInterface {

	public RandomVariableInterface getBrownianBridgeVariance(double time);
}
