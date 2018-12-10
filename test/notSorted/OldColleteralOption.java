package notSorted;

import montecarlo.interestrates.LIBORMarketModelWithBridge;
import montecarlo.interestrates.modelplugins.AbstractLiborCovarianceModelWithInterpolation;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationInterface;
import net.finmath.montecarlo.interestrate.products.AbstractLIBORMonteCarloProduct;
import net.finmath.stochastic.RandomVariableInterface;
import org.apache.commons.math3.distribution.NormalDistribution;


public class OldColleteralOption extends AbstractLIBORMonteCarloProduct {
	
	

	private final double  fixingDate;
	private final double  paymentDate;
	private final double  strike;
	private final boolean useAnalyticFormula;
	;

	public OldColleteralOption(
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
	
		RandomVariableInterface values;
		
		if(useAnalyticFormula) {
			//currently only working for paymentDate and fixingDate in same period.
			LIBORMarketModelWithBridge interpolationModel = (LIBORMarketModelWithBridge) model.getModel();
			
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
			
			//Calculate needed Variables (complete covariance & residuelRandomVariable)
			/*
			double discountCurveAtFixing    = interpolationModel.getDiscountCurve().getDiscountFactor(fixingDate);
			double discountCurveAtPayment   = interpolationModel.getDiscountCurve().getDiscountFactor(paymentDate);
			double discountCurveAtPeriodEnd = interpolationModel.getDiscountCurve().getDiscountFactor(liborEndTime);
			*/
			RandomVariableInterface	adjustedNumeraireAtEnd     = interpolationModel.getNumeraire(liborEndTime);
			RandomVariableInterface logLiborOverPeriod         = interpolationModel.getLIBOR(interpolationModel.getTimeIndex(interpolationModel.getLiborPeriod(liborStartIndex)), liborStartIndex)
																	.mult(liborEndTime - liborStartTime).add(1.0).log();
			RandomVariableInterface fixingExpLiborOverPeriod   = logLiborOverPeriod.mult((liborEndTime-fixingDate)/(liborEndTime-liborStartTime)).exp();
			RandomVariableInterface paymentExpLiborOverPeriod  = logLiborOverPeriod.mult((liborEndTime-paymentDate)/(liborEndTime-liborStartTime)).exp();

			//for now no adjustment is made
			RandomVariableInterface residuelRandomVariable     = fixingExpLiborOverPeriod.sub( paymentExpLiborOverPeriod.mult( 1 + strike * (paymentDate - fixingDate)) )
																	.div(adjustedNumeraireAtEnd.mult(paymentDate - fixingDate));
			
			RandomVariableInterface coefficientBridgeAtFixing  = adjustedNumeraireAtEnd.mult(paymentDate - fixingDate).invert();
			RandomVariableInterface coefficientBridgeAtPayment = adjustedNumeraireAtEnd.mult(paymentDate - fixingDate).div(1 + strike * (paymentDate - fixingDate)).invert();
			
			RandomVariableInterface completeVariance = coefficientBridgeAtFixing.pow(2.0).mult(bridgeVarianceAtFixing)
															.add(coefficientBridgeAtPayment.pow(2.0).mult(bridgeVarianceAtPayment))
															.sub(coefficientBridgeAtFixing.mult(coefficientBridgeAtPayment).mult(2*bridgeCovariance));
			//end Calculate needed Variables
			
			//Calculate Average
			double[] doubleValues = new double[residuelRandomVariable.size()];
			for (int pathIndex = 0; pathIndex < doubleValues.length; pathIndex++) {
				double residuelVariable = residuelRandomVariable.get(pathIndex);
				double varianceVariable = completeVariance.get(pathIndex);
				if(varianceVariable>0) {
				NormalDistribution normalDistribution = new NormalDistribution(0.0, Math.sqrt(varianceVariable));
				double cumulativeProbability = normalDistribution.cumulativeProbability(- residuelVariable);
				doubleValues[pathIndex] = residuelVariable * (1 - 2 * cumulativeProbability) 
											+ varianceVariable * normalDistribution.density( - residuelVariable);
				} else {
					doubleValues[pathIndex] = residuelVariable;
				}
			}
			values = new RandomVariable(evaluationTime, doubleValues);
			/*
			//temp test:
			double test = (interpolationModel.getBrownianBridge(liborStartIndex, fixingDate).mult(coefficientBridgeAtFixing.get(1)))
							.sub(interpolationModel.getBrownianBridge(liborStartIndex, paymentDate).mult(coefficientBridgeAtPayment.get(1))).getVariance();
			//RandomVariableInterface test = residuelRandomVariable.add(coefficientBridgeAtFixing.mult(interpolationModel.getBrownianBridge(liborStartIndex, fixingDate)))
			//								.sub(coefficientBridgeAtPayment.mult(interpolationModel.getBrownianBridge(liborStartIndex, paymentDate)));
			//RandomVariableInterface test = adjustedNumeraireAtEnd.mult(paymentDate-fixingDate).invert().mult(model.getLIBOR(fixingDate, fixingDate, liborEndTime).mult(liborEndTime-fixingDate).add(1)
			//		.sub(model.getLIBOR(paymentDate, paymentDate, liborEndTime).mult(liborEndTime-paymentDate).add(1).mult(1+strike*(paymentDate-fixingDate)) ) );
			System.out.println("test: " + test + "   " + completeVariance.get(1));
			System.out.println();
			//end temp test
			*/				
		} else {
			RandomVariableInterface numeraireAtFixing  	  				= model.getNumeraire(fixingDate);
			RandomVariableInterface numeraireAtPayment	 				= model.getNumeraire(paymentDate);
		
			RandomVariableInterface monteCarloProbabilitiesAtPayment    = model.getMonteCarloWeights(paymentDate);
			
			double periodLenght 					    			  	= paymentDate - fixingDate;
			
			values = numeraireAtPayment.div(numeraireAtFixing).sub(1).div(periodLenght).sub(strike).floor(0.0);
			values = values.div(numeraireAtPayment).mult(monteCarloProbabilitiesAtPayment);
		}
		RandomVariableInterface numeraireAtEvaluation				= model.getNumeraire(evaluationTime);
		RandomVariableInterface monteCarloProbabilitiesAtEvaluation = model.getMonteCarloWeights(evaluationTime);
		if(useAnalyticFormula) monteCarloProbabilitiesAtEvaluation = model.getRandomVariableForConstant(1.0);
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
