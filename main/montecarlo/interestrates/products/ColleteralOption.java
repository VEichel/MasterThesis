package montecarlo.interestrates.products;

import java.util.Arrays;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationInterface;
import net.finmath.montecarlo.interestrate.products.AbstractLIBORMonteCarloProduct;
import net.finmath.stochastic.RandomVariableInterface;


public class ColleteralOption extends AbstractLIBORMonteCarloProduct {
	
	

	private final double  fixingDate;
	private final double  paymentDate;
	private final double  strike;
	private final boolean useAnalyticFormula;
	;

	public ColleteralOption(
			double fixingDate,
			double paymentDate,
			double strike,
			boolean useAnalyticFormula
			) {
		super();
		this.fixingDate			= fixingDate;
		this.paymentDate		= paymentDate;
		this.strike				= strike;
		this.useAnalyticFormula = useAnalyticFormula;
	}



	@Override
	public RandomVariableInterface getValue(double evaluationTime, LIBORModelMonteCarloSimulationInterface model) throws CalculationException {
	
		
		RandomVariableInterface numeraireAtFixing  	  				= model.getNumeraire(fixingDate);
		RandomVariableInterface numeraireAtPayment	 				= model.getNumeraire(paymentDate);
		RandomVariableInterface numeraireAtEvaluation				= model.getNumeraire(evaluationTime);
		
		RandomVariableInterface monteCarloProbabilitiesAtEvaluation = model.getMonteCarloWeights(evaluationTime);
		RandomVariableInterface monteCarloProbabilitiesAtPayment    = model.getMonteCarloWeights(paymentDate);
		
		double periodLenght 					    			  	= paymentDate - fixingDate;
		
		RandomVariableInterface values = numeraireAtFixing.div(numeraireAtPayment).sub(1).div(periodLenght).sub(strike).floor(0.0);
								values = values.div(numeraireAtPayment).mult(monteCarloProbabilitiesAtPayment);
								
								values = values.mult(numeraireAtEvaluation).div(monteCarloProbabilitiesAtEvaluation);
								
		return values;
	}

	@Override
	public String toString() {
		return super.toString()
				+ "\n" + "fixingDate: " + fixingDate
				+ "\n" + "paymentDate: " + paymentDate
				+ "\n" + "strike: " + strike
				+ "\n" + "useAnalyticFormula: " + useAnalyticFormula;
	}
}
