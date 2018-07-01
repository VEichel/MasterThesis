package montecarlo.interestrates.products;

import java.util.Arrays;

import org.apache.commons.math3.analysis.function.Abs;
import org.apache.commons.math3.distribution.NormalDistribution;

import montecarlo.interestrates.LiborMarketModelWithBridgeInterpolation;
import montecarlo.interestrates.modelplugins.AbstractLiborCovarianceModelWithInterpolation;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORMarketModelInterface;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationInterface;
import net.finmath.montecarlo.interestrate.products.AbstractLIBORMonteCarloProduct;
import net.finmath.stochastic.RandomVariableInterface;


public class ColleteralOption extends AbstractLIBORMonteCarloProduct {
	
	

	private final double  fixingDate;
	private final double  paymentDate;
	private final double  strike;
	private       boolean useAnalyticFormula;
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
	
		RandomVariableInterface values;
		
		if(useAnalyticFormula) {
			//currently only working for paymentDate and fixingDate in same period.
			LiborMarketModelWithBridgeInterpolation interpolationModel = (LiborMarketModelWithBridgeInterpolation) model.getModel();
			
			int fixingIndex 	  			 = interpolationModel.getTimeIndex(fixingDate);
			int paymentIndex 	 			 = interpolationModel.getTimeIndex(paymentDate);
			int evaluationIndex   		 	 = interpolationModel.getTimeIndex(evaluationTime);
			
			int liborStartIndexFixing  			 = interpolationModel.getLiborPeriodIndex(fixingDate);
			int liborStartIndexPayment 			 = interpolationModel.getLiborPeriodIndex(paymentDate);
			
			//needed times for Fixing
			boolean fixingDateOnLiborStart = true;
			if(liborStartIndexFixing<0) {
				liborStartIndexFixing = (-liborStartIndexFixing - 1) - 1;
				fixingDateOnLiborStart = false;
			}
		    int liborEndIndexFixing	 		 = liborStartIndexFixing + 1;		    
		    double liborStartTimeFixing		 = interpolationModel.getLiborPeriod(liborStartIndexFixing);
		    double liborEndTimeFixing	 	 = interpolationModel.getLiborPeriod(liborEndIndexFixing);
		    
		    //needed times for payment
		    boolean paymentDateOnLiborStart = true;
			if(liborStartIndexPayment<0) {
				liborStartIndexPayment = (-liborStartIndexPayment - 1) - 1;
				paymentDateOnLiborStart = false;
			}
			boolean samePeriod = true;
			if(liborStartIndexFixing != liborStartIndexPayment) {
				samePeriod = false;
			}
		    int liborEndIndexPayment 		 = liborStartIndexPayment + 1;		    
		    double liborStartTimePayment	 = interpolationModel.getLiborPeriod(liborStartIndexPayment);
		    double liborEndTimePayment	 	 = interpolationModel.getLiborPeriod(liborEndIndexPayment);
			

		    
		    if(fixingDateOnLiborStart && paymentDateOnLiborStart) {
		    	this.useAnalyticFormula = false;
		    	return this.getValue(evaluationTime, model);
		    }
		    
		    int periodStartToFixingIndex    = fixingIndex - interpolationModel.getTimeIndex(liborStartTimeFixing);
		    int periodStartToPaymentIndex   = paymentIndex  - interpolationModel.getTimeIndex(liborStartTimePayment);

		   
		    
			//Calculate BB Variance:
		    AbstractLiborCovarianceModelWithInterpolation covarianceModel = (AbstractLiborCovarianceModelWithInterpolation) interpolationModel.getCovarianceModel();
		    double bridgeVarianceAtFixing  = 0.0;
			double bridgeVarianceAtPayment = 0.0;
			double bridgeCovariance        = 0.0;
			
		    if(samePeriod) {
				for(int i = 1; i < periodStartToPaymentIndex + 1 ; i++) {
					double varianceI    = covarianceModel.getVarianceForInterpolation(i + liborStartIndexFixing - 1).getAverage();
					double timeI		= interpolationModel.getTime(i + liborStartIndexFixing);	
					double timeBeforeI  = interpolationModel.getTime(i + liborStartIndexFixing - 1);
					if(i < periodStartToFixingIndex + 1) {
						bridgeVarianceAtFixing  += varianceI * ((timeI - timeBeforeI) - (timeI - timeBeforeI)*(timeI - timeBeforeI)/(liborEndTimeFixing - liborStartTimeFixing));
					} else {
						bridgeVarianceAtPayment += varianceI * ((timeI - timeBeforeI) - (timeI - timeBeforeI)*(timeI - timeBeforeI)/(liborEndTimeFixing - liborStartTimeFixing));
					}
					for (int j = 1; j < periodStartToPaymentIndex + 1; j++) {
						double varianceJ    = covarianceModel.getVarianceForInterpolation(j-1).getAverage();
						double timeJ		= interpolationModel.getTime(j + liborStartIndexFixing);	
						double timeBeforeJ  = interpolationModel.getTime(j + liborStartIndexFixing - 1);	
						if(i < periodStartToFixingIndex + 1 && j < periodStartToFixingIndex + 1 && i != j) {
							bridgeVarianceAtFixing  -= Math.sqrt(varianceI) * Math.sqrt(varianceJ) * (timeI - timeBeforeI) * (timeJ - timeBeforeJ) / (liborEndTimeFixing - liborStartTimeFixing);
						} else if(j < periodStartToFixingIndex + 1 && i != j){
							bridgeCovariance        -= Math.sqrt(varianceI) * Math.sqrt(varianceJ) * (timeI - timeBeforeI) * (timeJ - timeBeforeJ) / (liborEndTimeFixing - liborStartTimeFixing);
							bridgeVarianceAtPayment -= Math.sqrt(varianceI) * Math.sqrt(varianceJ) * (timeI - timeBeforeI) * (timeJ - timeBeforeJ) / (liborEndTimeFixing - liborStartTimeFixing);
						} else if(i != j) {
							bridgeVarianceAtPayment -= Math.sqrt(varianceI) * Math.sqrt(varianceJ) * (timeI - timeBeforeI) * (timeJ - timeBeforeJ) / (liborEndTimeFixing - liborStartTimeFixing);
						}
					}
				}
				double timeScalingFactorSqaured      = Math.pow(covarianceModel.getEvaluationTimeScalingFactor(evaluationTime).getAverage(), 2);
				bridgeVarianceAtPayment = (bridgeVarianceAtPayment + bridgeVarianceAtFixing) * timeScalingFactorSqaured;
				bridgeCovariance		= (bridgeCovariance        + bridgeVarianceAtFixing) * timeScalingFactorSqaured;
				bridgeVarianceAtFixing  = bridgeVarianceAtFixing * timeScalingFactorSqaured;
				
		    } else {
				
				for(int i = 1; i < periodStartToFixingIndex + 1 ; i++) {
					double varianceIFixing     = covarianceModel.getVarianceForInterpolation(i + liborStartIndexFixing - 1).getAverage();
					double timeIFixing		   = interpolationModel.getTime(i + liborStartIndexFixing);	
					double timeBeforeIFixing   = interpolationModel.getTime(i + liborStartIndexFixing - 1);
					
					bridgeVarianceAtFixing  += varianceIFixing * ((timeIFixing - timeBeforeIFixing) - (timeIFixing - timeBeforeIFixing)*(timeIFixing - timeBeforeIFixing)/(liborEndTimeFixing - liborStartTimeFixing));
					
					for (int j = 1; j < i; j++) {
						double varianceJFixing    = covarianceModel.getVarianceForInterpolation(j + liborStartIndexFixing - 1).getAverage();
						double timeJFixing		  = interpolationModel.getTime(j + liborStartIndexFixing);	
						double timeBeforeJFixing  = interpolationModel.getTime(j + liborStartIndexFixing - 1);	
						
					    bridgeVarianceAtFixing  -=  2 * Math.sqrt(varianceIFixing) * Math.sqrt(varianceJFixing) * (timeIFixing - timeBeforeIFixing) * (timeJFixing - timeBeforeJFixing) / (liborEndTimeFixing - liborStartTimeFixing);
					}
				}
				
				for(int i = 1; i < periodStartToPaymentIndex + 1 ; i++) {
					double varianceIPayment    = covarianceModel.getVarianceForInterpolation(i + liborStartIndexPayment - 1).getAverage();
					double timeIPayment		   = interpolationModel.getTime(i + liborStartIndexPayment);	
					double timeBeforeIPayment  = interpolationModel.getTime(i + liborStartIndexPayment - 1);
					
					bridgeVarianceAtPayment += varianceIPayment * ((timeIPayment - timeBeforeIPayment) - (timeIPayment - timeBeforeIPayment)*(timeIPayment - timeBeforeIPayment)/(liborEndTimePayment - liborStartTimePayment));
					
					for (int j = 1; j < i; j++) {
						double varianceJPayment    = covarianceModel.getVarianceForInterpolation(j + liborStartIndexPayment - 1).getAverage();
						double timeJPayment		  = interpolationModel.getTime(j + liborStartIndexPayment);	
						double timeBeforeJPayment  = interpolationModel.getTime(j + liborStartIndexPayment - 1);	
						
					    bridgeVarianceAtPayment  -=  2 * Math.sqrt(varianceIPayment) * Math.sqrt(varianceJPayment) * (timeIPayment - timeBeforeIPayment) * (timeJPayment - timeBeforeJPayment) / (liborEndTimePayment - liborStartTimePayment);
					}
				}
		    	
		    }
			//Finished BB Variance calculations
			
		  //temp test 3:
			double test = interpolationModel.getBrownianBridge(liborStartIndexFixing, fixingDate).getVariance();
			//RandomVariableInterface test = residuelRandomVariable.add(coefficientBridgeAtFixing.mult(interpolationModel.getBrownianBridge(liborStartIndex, fixingDate)))
			//								.sub(coefficientBridgeAtPayment.mult(interpolationModel.getBrownianBridge(liborStartIndex, paymentDate)));
			//RandomVariableInterface test = adjustedNumeraireAtEnd.mult(paymentDate-fixingDate).invert().mult(model.getLIBOR(fixingDate, fixingDate, liborEndTime).mult(liborEndTime-fixingDate).add(1)
			//		.sub(model.getLIBOR(paymentDate, paymentDate, liborEndTime).mult(liborEndTime-paymentDate).add(1).mult(1+strike*(paymentDate-fixingDate)) ) );
			System.out.println("test: " + (Math.abs(bridgeVarianceAtFixing - test) / bridgeVarianceAtFixing)*100 + "%" );
			System.out.println();
			//end temp test 3
		    
		    
			//Calculate needed Variables (complete covariance & residuelRandomVariable)
			RandomVariableInterface	adjustedNumeraireAtEndFixing  = interpolationModel.getNumeraire(liborEndTimeFixing);
			RandomVariableInterface	adjustedNumeraireAtEndPayment =  samePeriod ? adjustedNumeraireAtEndFixing : 
																		interpolationModel.getNumeraire(liborEndTimePayment);
			
			
			double numeraireAdjustmentForFixing  = 1.0;
			double numeraireAdjustmentForPayment = 1.0;
			
			/*
			//adjustment not correct jet. N(t_1) N(t_2) is used and the analytic average is taken.
			if(interpolationModel.getDiscountCurve() != null) {
				double discountCurveAtFixing    = interpolationModel.getDiscountCurve().getDiscountFactor(fixingDate);
				double discountCurveAtPayment   = interpolationModel.getDiscountCurve().getDiscountFactor(paymentDate);
				double discountCurveAtPeriodEnd = interpolationModel.getDiscountCurve().getDiscountFactor(liborEndTime);
				
				numeraireAdjustmentForFixing  = discountCurveAtPeriodEnd *(interpolationModel.getNumeraire(fixingDate).invert().getAverage()) / (discountCurveAtFixing * adjustedNumeraireAtEnd.invert().getAverage());
				numeraireAdjustmentForPayment = discountCurveAtPeriodEnd *(interpolationModel.getNumeraire(paymentDate).invert().getAverage()) / (discountCurveAtPayment * adjustedNumeraireAtEnd.invert().getAverage());
				System.out.println(numeraireAdjustmentForFixing);
			}*/
			
			RandomVariableInterface logLiborOverPeriodFixing   = interpolationModel.getLIBOR(interpolationModel.getTimeIndex(liborStartTimeFixing), liborStartIndexFixing)
																	.mult(liborEndTimeFixing - liborStartTimeFixing).add(1.0).log();
			RandomVariableInterface logLiborOverPeriodPayment  = samePeriod ? logLiborOverPeriodFixing : interpolationModel.getLIBOR(interpolationModel.getTimeIndex(liborStartTimePayment), liborStartIndexPayment)
																											.mult(liborEndTimePayment - liborStartTimePayment).add(1.0).log();
															
			RandomVariableInterface fixingExpLiborOverPeriod   = logLiborOverPeriodFixing.mult((liborEndTimeFixing-fixingDate)/(liborEndTimeFixing-liborStartTimeFixing)).exp();
			RandomVariableInterface paymentExpLiborOverPeriod  = logLiborOverPeriodPayment.mult((liborEndTimePayment-paymentDate)/(liborEndTimePayment-liborStartTimePayment)).exp();

			RandomVariableInterface coefficientBridgeAtFixing  = adjustedNumeraireAtEndFixing.mult((paymentDate - fixingDate) * numeraireAdjustmentForFixing).invert();
			RandomVariableInterface coefficientBridgeAtPayment = adjustedNumeraireAtEndPayment.mult((paymentDate - fixingDate) * numeraireAdjustmentForPayment).div(1 + strike * (paymentDate - fixingDate)).invert();
			
			
			//analytic ForwardCurve adjustment
			if(interpolationModel.getForwardRateCurve() != null) {
				double forwardCurveFixingEnd    = 1 + interpolationModel.getForwardRateCurve().getForward(interpolationModel.getAnalyticModel(), fixingDate, liborEndTimeFixing - fixingDate) * (liborEndTimeFixing - fixingDate);
				double forwardCurvePaymentEnd   = 1 + interpolationModel.getForwardRateCurve().getForward(interpolationModel.getAnalyticModel(), paymentDate, liborEndTimePayment - paymentDate) * (liborEndTimePayment - paymentDate);
				double forwardCurveLongFixing   = 1 + interpolationModel.getForwardRateCurve().getForward(interpolationModel.getAnalyticModel(), liborStartTimeFixing, liborEndTimeFixing - liborStartTimeFixing) * (liborEndTimeFixing - liborStartTimeFixing);
				double forwardCurveLongPayment  = samePeriod ? forwardCurveLongFixing : 
													1 + interpolationModel.getForwardRateCurve().getForward(interpolationModel.getAnalyticModel(), liborStartTimePayment, liborEndTimePayment - liborStartTimePayment) * (liborEndTimePayment - liborStartTimePayment); 
				
				double forwardCurveFixingEndLog  = Math.exp(Math.log(forwardCurveLongFixing)*(liborEndTimeFixing - fixingDate)/(liborEndTimeFixing - liborStartTimeFixing));
				double forwardCurvePaymentEndLog = Math.exp(Math.log(forwardCurveLongPayment)*(liborEndTimePayment - paymentDate)/(liborEndTimePayment - liborStartTimePayment)); 
				
				double forwardFixingAdjustment  = forwardCurveFixingEnd  / forwardCurveFixingEndLog;
				double forwardPaymentAdjustment = forwardCurvePaymentEnd / forwardCurvePaymentEndLog;
				
				coefficientBridgeAtFixing  = coefficientBridgeAtFixing.mult(forwardFixingAdjustment);
				coefficientBridgeAtPayment = coefficientBridgeAtPayment.mult(forwardPaymentAdjustment);
				fixingExpLiborOverPeriod   = fixingExpLiborOverPeriod.mult(forwardFixingAdjustment);
				paymentExpLiborOverPeriod  = paymentExpLiborOverPeriod.mult(forwardPaymentAdjustment);
			}
			
			RandomVariableInterface completeVariance 	   = coefficientBridgeAtFixing.pow(2.0).mult(bridgeVarianceAtFixing)
																.add(coefficientBridgeAtPayment.pow(2.0).mult(bridgeVarianceAtPayment));
			if(samePeriod) completeVariance = completeVariance.sub(coefficientBridgeAtFixing.mult(coefficientBridgeAtPayment).mult(2*bridgeCovariance));
			
			RandomVariableInterface residuelRandomVariable = fixingExpLiborOverPeriod.div( adjustedNumeraireAtEndFixing.mult((paymentDate - fixingDate)*numeraireAdjustmentForFixing) )
					.sub( paymentExpLiborOverPeriod.mult( (1 + strike * (paymentDate - fixingDate)) )
					.div(adjustedNumeraireAtEndPayment.mult((paymentDate - fixingDate)*numeraireAdjustmentForPayment)) );

			//end Calculate needed Variables
			
			
			//Calculate Average
			long start = System.currentTimeMillis();
			
			double[] doubleValues = new double[residuelRandomVariable.size()];
			for (int pathIndex = 0; pathIndex < doubleValues.length; pathIndex++) {
				double residuelVariable = residuelRandomVariable.get(pathIndex);
				double varianceVariable = completeVariance.get(pathIndex);
				NormalDistribution normalDistribution = new NormalDistribution(0.0, Math.sqrt(varianceVariable));
				double cumulativeProbability = normalDistribution.cumulativeProbability(- residuelVariable);
				doubleValues[pathIndex] = residuelVariable * (1 - 2 * cumulativeProbability) 
											+ varianceVariable * normalDistribution.density( - residuelVariable);
			}
			values = new RandomVariable(evaluationTime, doubleValues);
			System.out.println(System.currentTimeMillis() - start);
			
			/*
			//temp test2:
			NormalDistribution nd = new NormalDistribution(0.0, Math.sqrt(completeVariance.getAverage()));
			double cumulativeProbability = nd.cumulativeProbability(- residuelRandomVariable.getAverage());
			values = interpolationModel.getRandomVariableForConstant(residuelRandomVariable.getAverage() * (1 - 2 * cumulativeProbability) 
										+ completeVariance.getAverage() * nd.density( - residuelRandomVariable.getAverage()));
			//end temp test2
			*/
			
			/*
			//temp test:
			double test = (interpolationModel.getBrownianBridge(liborStartIndexFixing, fixingDate).mult(coefficientBridgeAtFixing.get(1)))
							.sub(interpolationModel.getBrownianBridge(liborStartIndexPayment, paymentDate).mult(coefficientBridgeAtPayment.get(1))).getVariance();
			//RandomVariableInterface test = residuelRandomVariable.add(coefficientBridgeAtFixing.mult(interpolationModel.getBrownianBridge(liborStartIndex, fixingDate)))
			//								.sub(coefficientBridgeAtPayment.mult(interpolationModel.getBrownianBridge(liborStartIndex, paymentDate)));
			//RandomVariableInterface test = adjustedNumeraireAtEnd.mult(paymentDate-fixingDate).invert().mult(model.getLIBOR(fixingDate, fixingDate, liborEndTime).mult(liborEndTime-fixingDate).add(1)
			//		.sub(model.getLIBOR(paymentDate, paymentDate, liborEndTime).mult(liborEndTime-paymentDate).add(1).mult(1+strike*(paymentDate-fixingDate)) ) );
			System.out.println("test: " + (Math.abs(completeVariance.get(1) - test) / completeVariance.get(1))*100 + "%" );
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
