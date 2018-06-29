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
			
			int fixingIndex 	  			 = interpolationModel.getTimeIndex(fixingDate);
			int paymentIndex 	 			 = interpolationModel.getTimeIndex(paymentDate);
			int evaluationIndex   		 	 = interpolationModel.getTimeIndex(evaluationTime);
			
			int liborStartIndex  			 = interpolationModel.getLiborPeriodIndex(fixingDate);
			if(liborStartIndex<0) liborStartIndex = (-liborStartIndex - 1) - 1;
		    int liborEndIndex	 		 	 = liborStartIndex + 1;
		    double liborStartTime			 = interpolationModel.getLiborPeriod(liborStartIndex);
		    double liborEndTime	 			 = interpolationModel.getLiborPeriod(liborEndIndex);
		    
		    int periodStartToPaymentIndex    = paymentIndex - interpolationModel.getTimeIndex(liborStartTime);
		    int periodStartToFixingIndex     = fixingIndex  - interpolationModel.getTimeIndex(liborStartTime);
		    
			//Calculate BB Variance:
			AbstractLiborCovarianceModelWithInterpolation covarianceModel = (AbstractLiborCovarianceModelWithInterpolation) interpolationModel.getCovarianceModel();
			double bridgeVarianceAtFixing  = 0.0;
			double bridgeVarianceAtPayment = 0.0;
			double bridgeCovariance        = 0.0;
			for(int i = 1; i < periodStartToPaymentIndex + 1 ; i++) {
				double varianceI    = covarianceModel.getVarianceForInterpolation(i-1).getAverage();
				double timeI		= interpolationModel.getTime(i + liborStartIndex);	
				double timeBeforeI  = interpolationModel.getTime(i + liborStartIndex - 1);
				if(i < periodStartToFixingIndex + 1) {
					bridgeVarianceAtFixing  += varianceI * ((timeI - timeBeforeI) - (timeI - timeBeforeI)*(timeI - timeBeforeI)/(liborEndTime - liborStartTime));
				} else {
					bridgeVarianceAtPayment += varianceI * ((timeI - timeBeforeI) - (timeI - timeBeforeI)*(timeI - timeBeforeI)/(liborEndTime - liborStartTime));
				}
				for (int j = 1; j < periodStartToPaymentIndex + 1; j++) {
					double varianceJ    = covarianceModel.getVarianceForInterpolation(j-1).getAverage();
					double timeJ		= interpolationModel.getTime(j + liborStartIndex);	
					double timeBeforeJ  = interpolationModel.getTime(j + liborStartIndex - 1);	
					if(i < periodStartToFixingIndex + 1 && j < periodStartToFixingIndex + 1 && i != j) {
						bridgeVarianceAtFixing  -= Math.sqrt(varianceI) * Math.sqrt(varianceJ) * (timeI - timeBeforeI) * (timeJ - timeBeforeJ) / (liborEndTime - liborStartTime);
					} else if(j < periodStartToFixingIndex + 1 && i != j){
						bridgeCovariance        -= Math.sqrt(varianceI) * Math.sqrt(varianceJ) * (timeI - timeBeforeI) * (timeJ - timeBeforeJ) / (liborEndTime - liborStartTime);
						bridgeVarianceAtPayment -= Math.sqrt(varianceI) * Math.sqrt(varianceJ) * (timeI - timeBeforeI) * (timeJ - timeBeforeJ) / (liborEndTime - liborStartTime);
					} else if(i != j) {
						bridgeVarianceAtPayment -= Math.sqrt(varianceI) * Math.sqrt(varianceJ) * (timeI - timeBeforeI) * (timeJ - timeBeforeJ) / (liborEndTime - liborStartTime);
					}
				}
			}
			double timeScalingFactorSqaured      = Math.pow(covarianceModel.getEvaluationTimeScalingFactor(evaluationTime).getAverage(), 2);
			bridgeVarianceAtPayment = (bridgeVarianceAtPayment + bridgeVarianceAtFixing) * timeScalingFactorSqaured;
			bridgeCovariance		= (bridgeCovariance        + bridgeVarianceAtFixing) * timeScalingFactorSqaured;
			bridgeVarianceAtFixing  = bridgeVarianceAtFixing * timeScalingFactorSqaured;
			//Finished BB Variance calculations
			
			
			
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
