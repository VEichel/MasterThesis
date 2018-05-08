package montecarlo.interestrate.modelplugins;

import java.util.Arrays;

import net.finmath.montecarlo.AbstractRandomVariableFactory;
import net.finmath.montecarlo.interestrate.TermStructureModelInterface;
import net.finmath.montecarlo.interestrate.modelplugins.AbstractLIBORCovarianceModel;
import net.finmath.montecarlo.interestrate.modelplugins.AbstractLIBORCovarianceModelParametric;
import net.finmath.montecarlo.interestrate.modelplugins.TermStructureCovarianceModelParametric;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretizationInterface;


/**non correlated
 * 
 * @author vince
 *
 */
public class LiborWithInterpolationCovarianceModelFromArray extends AbstractLiborWithInterpolationCovarianceModel{

	double[] interpolationFactorArray;
	private final AbstractRandomVariableFactory	randomVariableFactory;
	
	public LiborWithInterpolationCovarianceModelFromArray(AbstractLIBORCovarianceModel liborCovarianceModel, double[] interpolationFactorArray) {
		super(liborCovarianceModel);
		this.interpolationFactorArray = interpolationFactorArray;
		this.randomVariableFactory = null;
	}

	@Override
	public RandomVariableInterface[] getFactorLoadingForInterpolation(int interpolationIndex) {
		return new RandomVariableInterface[] { randomVariableFactory.createRandomVariable(interpolationFactorArray[interpolationIndex]) };
	}

	@Override
	public double[] getParameter() {
		double[] liborParameter = ((AbstractLIBORCovarianceModelParametric) liborCovarianceModel).getParameter();
		double[] parameters = new double[interpolationFactorArray.length + liborParameter.length];
		for (int index = 0; index < interpolationFactorArray.length; index++) {
			parameters[index] = interpolationFactorArray[index];
		}
		for (int index = interpolationFactorArray.length; index < parameters.length; index++) {
			parameters[index] = interpolationFactorArray[index];
		}
		return parameters;
	}

	@Override
	public AbstractLiborWithInterpolationCovarianceModel getCloneWithModifiedParameter(double[] parameters) {
		return new LiborWithInterpolationCovarianceModelFromArray(liborCovarianceModel, parameters);
	}
}
