package montecarlo.interestrates.products;

import java.util.Arrays;

import org.apache.commons.math3.analysis.function.Abs;

import montecarlo.interestrates.LiborMarketModelWithBridgeInterpolation;
import montecarlo.interestrates.modelplugins.AbstractLiborCovarianceModelWithInterpolation;
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
	
		if(useAnalyticFormula) {
			//currently only working for paymentDate and fixingDate in same period.
			LiborMarketModelWithBridgeInterpolation interpolationModel = (LiborMarketModelWithBridgeInterpolation) model.getModel();
			
			int fixingIndex 	= interpolationModel.getTimeIndex(fixingDate);
			int paymentIndex 	= interpolationModel.getTimeIndex(paymentDate);
			int evaluationIndex = interpolationModel.getTimeIndex(evaluationTime);
			//Calculate BB Variance:
			AbstractLiborCovarianceModelWithInterpolation covarianceModel = (AbstractLiborCovarianceModelWithInterpolation) interpolationModel.getCovarianceModel();
			double bridgeVarianceAtFixing  = 0.0;
			double bridgeVarianceAtPayment = 0.0;
			double brdgeCovariance         = 0.0;
			
			for(int i = 1; i < paymentIndex + 1; i++) {
				double varianceI    = covarianceModel.getVarianceForInterpolation(i).getAverage();
				int    timeI			
				int    timeBeforeI
				if(i < fixingIndex + 1) {
					bridgeVarianceAtFixing += varianceI*()
				}
				for (int j = 0; j < paymentIndex + 1; j++) {
					
				}
			}
			
		}
		RandomVariableInterface numeraireAtFixing  	  				= model.getNumeraire(fixingDate);
		RandomVariableInterface numeraireAtPayment	 				= model.getNumeraire(paymentDate);
		RandomVariableInterface numeraireAtEvaluation				= model.getNumeraire(evaluationTime);
		
		RandomVariableInterface monteCarloProbabilitiesAtEvaluation = model.getMonteCarloWeights(evaluationTime);
		RandomVariableInterface monteCarloProbabilitiesAtPayment    = model.getMonteCarloWeights(paymentDate);
		
		double periodLenght 					    			  	= paymentDate - fixingDate;
		
		RandomVariableInterface values = numeraireAtPayment.div(numeraireAtFixing).sub(1).div(periodLenght).sub(strike).floor(0.0);
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
