package montecarlo.interestrates.products;

//TODO: for fixingOnLiborStart use change liborStartToEnd and dont use any libor. Numeraire is enough!

import java.util.Arrays;

import org.apache.commons.math3.analysis.function.Abs;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.genetics.GeneticAlgorithm;
import org.apache.commons.math3.special.Erf;

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
				double timeScalingFactorSqauredFixing      = Math.pow(covarianceModel.getEvaluationTimeScalingFactor(evaluationTime).getAverage(), 2);
				bridgeVarianceAtPayment = (bridgeVarianceAtPayment + bridgeVarianceAtFixing) * timeScalingFactorSqauredFixing;
				bridgeCovariance		= (bridgeCovariance        + bridgeVarianceAtFixing) * timeScalingFactorSqauredFixing;
				bridgeVarianceAtFixing  = bridgeVarianceAtFixing * timeScalingFactorSqauredFixing;
				
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
						
					    bridgeVarianceAtFixing   -=  2 * Math.sqrt(varianceIFixing) * Math.sqrt(varianceJFixing) * (timeIFixing - timeBeforeIFixing) * (timeJFixing - timeBeforeJFixing) / (liborEndTimeFixing - liborStartTimeFixing);
					}
				}
				
				for(int i = 1; i < periodStartToPaymentIndex + 1 ; i++) {
					double varianceIPayment    = covarianceModel.getVarianceForInterpolation(i + liborStartIndexPayment - 1).getAverage();
					double timeIPayment		   = interpolationModel.getTime(i + liborStartIndexPayment);	
					double timeBeforeIPayment  = interpolationModel.getTime(i + liborStartIndexPayment - 1);
					
					bridgeVarianceAtPayment += varianceIPayment * ((timeIPayment - timeBeforeIPayment) - (timeIPayment - timeBeforeIPayment)*(timeIPayment - timeBeforeIPayment)/(liborEndTimePayment - liborStartTimePayment));
					
					for (int j = 1; j < i; j++) {
						double varianceJPayment    = covarianceModel.getVarianceForInterpolation(j + liborStartIndexPayment - 1).getAverage();
						double timeJPayment		   = interpolationModel.getTime(j + liborStartIndexPayment);	
						double timeBeforeJPayment  = interpolationModel.getTime(j + liborStartIndexPayment - 1);	
						
					    bridgeVarianceAtPayment   -=  2 * Math.sqrt(varianceIPayment) * Math.sqrt(varianceJPayment) * (timeIPayment - timeBeforeIPayment) * (timeJPayment - timeBeforeJPayment) / (liborEndTimePayment - liborStartTimePayment);
					}
				}
				double timeScalingFactorSqauredFixing   = Math.pow(covarianceModel.getEvaluationTimeScalingFactor(evaluationTime).getAverage(), 2);
				bridgeVarianceAtFixing  *= timeScalingFactorSqauredFixing;
				bridgeVarianceAtPayment *= timeScalingFactorSqauredFixing;
		    }
			//Finished BB Variance calculations
			
		    //temporär:
		    bridgeVarianceAtFixing = interpolationModel.getBrownianBridge(liborStartIndexFixing, fixingDate).getVariance();
		    bridgeVarianceAtPayment = interpolationModel.getBrownianBridge(liborStartIndexPayment, paymentDate).getVariance();
		    bridgeCovariance = interpolationModel.getBrownianBridge(liborStartIndexFixing, fixingDate).mult(interpolationModel.getBrownianBridge(liborStartIndexPayment, paymentDate)).getAverage();
		    /*
		   //temp test 3:
			double test = covarianceModel.getEvaluationTimeScalingFactor(evaluationTime).mult(interpolationModel.getBrownianBridge(liborStartIndexPayment, paymentDate)).getVariance();
			//RandomVariableInterface test = residuelRandomVariable.add(coefficientBridgeAtFixing.mult(interpolationModel.getBrownianBridge(liborStartIndex, fixingDate)))
			//								.sub(coefficientBridgeAtPayment.mult(interpolationModel.getBrownianBridge(liborStartIndex, paymentDate)));
			//RandomVariableInterface test = adjustedNumeraireAtEnd.mult(paymentDate-fixingDate).invert().mult(model.getLIBOR(fixingDate, fixingDate, liborEndTime).mult(liborEndTime-fixingDate).add(1)
			//		.sub(model.getLIBOR(paymentDate, paymentDate, liborEndTime).mult(liborEndTime-paymentDate).add(1).mult(1+strike*(paymentDate-fixingDate)) ) );
			System.out.println("test 3: " + "\tcalculated: " + bridgeVarianceAtPayment +"\t not calc: " + test +"\t Error: " + (Math.abs(bridgeVarianceAtPayment - test) / bridgeVarianceAtPayment)*100 + "%" );
			System.out.println();
			
			double test4 = covarianceModel.getEvaluationTimeScalingFactor(evaluationTime).mult(interpolationModel.getBrownianBridge(liborStartIndexPayment, paymentDate)
					.mult(interpolationModel.getBrownianBridge(liborStartIndexPayment, fixingDate))).getAverage();
			//RandomVariableInterface test = residuelRandomVariable.add(coefficientBridgeAtFixing.mult(interpolationModel.getBrownianBridge(liborStartIndex, fixingDate)))
			//								.sub(coefficientBridgeAtPayment.mult(interpolationModel.getBrownianBridge(liborStartIndex, paymentDate)));
			//RandomVariableInterface test = adjustedNumeraireAtEnd.mult(paymentDate-fixingDate).invert().mult(model.getLIBOR(fixingDate, fixingDate, liborEndTime).mult(liborEndTime-fixingDate).add(1)
			//		.sub(model.getLIBOR(paymentDate, paymentDate, liborEndTime).mult(liborEndTime-paymentDate).add(1).mult(1+strike*(paymentDate-fixingDate)) ) );
			System.out.println("test 4: " + "\tcalculated: " + bridgeCovariance +"\t not calc: " + test4 +"\t Error: " + (Math.abs(bridgeCovariance - test4) / bridgeCovariance)*100 + "%" );
			System.out.println();
			//end temp test 3
		    */
		    
			//Calculate needed Variables (complete covariance & residuelRandomVariable)
			RandomVariableInterface	numeraireAtEndFixing       = interpolationModel.getUnAdjustedNumeraire(liborEndTimeFixing);
			RandomVariableInterface	numeraireAtEndPayment      =  samePeriod ? numeraireAtEndFixing : interpolationModel.getUnAdjustedNumeraire(liborEndTimePayment);
			
			RandomVariableInterface liborOverPeriodFixing      = interpolationModel.getLIBOR(interpolationModel.getTimeIndex(liborStartTimeFixing), liborStartIndexFixing)
																	.mult(liborEndTimeFixing - liborStartTimeFixing).add(1.0);
			RandomVariableInterface liborOverPeriodPayment 	   = samePeriod ? liborOverPeriodFixing : interpolationModel.getLIBOR(interpolationModel.getTimeIndex(liborStartTimePayment), liborStartIndexPayment)
																											.mult(liborEndTimePayment - liborStartTimePayment).add(1.0);
			RandomVariableInterface fixingExpLiborOverPeriod   = fixingDateOnLiborStart ? liborOverPeriodFixing :
																	liborOverPeriodFixing.log().mult((liborEndTimeFixing-fixingDate)/(liborEndTimeFixing-liborStartTimeFixing)).exp();
			RandomVariableInterface paymentExpLiborOverPeriod  = paymentDateOnLiborStart ? liborOverPeriodPayment :
																	liborOverPeriodPayment.log().mult((liborEndTimePayment-paymentDate)/(liborEndTimePayment-liborStartTimePayment)).exp();
			//multi curve adjustment
			if(interpolationModel.getDiscountCurve() != null && false) {
				double numeraireAtFixingAdjustment           = fixingDateOnLiborStart  ? numeraireAtEndFixing.invert().getAverage()  / interpolationModel.getDiscountCurve().getDiscountFactor(interpolationModel.getAnalyticModel(), fixingDate)  :
																			fixingExpLiborOverPeriod.div(numeraireAtEndFixing).getAverage() / interpolationModel.getDiscountCurve().getDiscountFactor(interpolationModel.getAnalyticModel(), fixingDate);
				double numeraireAtPaymentAdjustment          = paymentDateOnLiborStart ? numeraireAtEndPayment.invert().getAverage() / interpolationModel.getDiscountCurve().getDiscountFactor(interpolationModel.getAnalyticModel(), paymentDate) :
																			paymentExpLiborOverPeriod.div(numeraireAtEndPayment).getAverage() / interpolationModel.getDiscountCurve().getDiscountFactor(interpolationModel.getAnalyticModel(), paymentDate);
				numeraireAtEndFixing  = numeraireAtEndFixing.mult(numeraireAtFixingAdjustment);
				numeraireAtEndPayment = numeraireAtEndPayment.mult(numeraireAtPaymentAdjustment);
			}
			//end multi curve adjustment
			
			RandomVariableInterface coefficientBridgeAtFixing  = numeraireAtEndFixing.mult(paymentDate - fixingDate).invert();
			RandomVariableInterface coefficientBridgeAtPayment = numeraireAtEndPayment.mult(paymentDate - fixingDate).div(1 + strike * (paymentDate - fixingDate)).invert();
			
			//analytic ForwardCurve adjustment
			if(interpolationModel.getForwardRateCurve() != null && false) {
				double forwardCurveFixingEnd    = 1 + interpolationModel.getForwardRateCurve().getForward(interpolationModel.getAnalyticModel(), fixingDate, liborEndTimeFixing - fixingDate) * (liborEndTimeFixing - fixingDate);
				double forwardCurvePaymentEnd   = 1 + interpolationModel.getForwardRateCurve().getForward(interpolationModel.getAnalyticModel(), paymentDate, liborEndTimePayment - paymentDate) * (liborEndTimePayment - paymentDate);
				double forwardCurveLongFixing   = 1 + interpolationModel.getForwardRateCurve().getForward(interpolationModel.getAnalyticModel(), liborStartTimeFixing, liborEndTimeFixing - liborStartTimeFixing) * (liborEndTimeFixing - liborStartTimeFixing);
				double forwardCurveLongPayment  = samePeriod ? forwardCurveLongFixing : 
													1 + interpolationModel.getForwardRateCurve().getForward(interpolationModel.getAnalyticModel(), liborStartTimePayment, liborEndTimePayment - liborStartTimePayment) * (liborEndTimePayment - liborStartTimePayment); 
				
				double forwardCurveFixingEndLog  = Math.exp(Math.log(forwardCurveLongFixing)*(liborEndTimeFixing - fixingDate)/(liborEndTimeFixing - liborStartTimeFixing));
				double forwardCurvePaymentEndLog = Math.exp(Math.log(forwardCurveLongPayment)*(liborEndTimePayment - paymentDate)/(liborEndTimePayment - liborStartTimePayment)); 
				
				double forwardFixingAdjustment   = forwardCurveFixingEnd  / forwardCurveFixingEndLog;
				double forwardPaymentAdjustment  = forwardCurvePaymentEnd / forwardCurvePaymentEndLog;
				
				if(!fixingDateOnLiborStart) {
					fixingExpLiborOverPeriod    = fixingExpLiborOverPeriod.mult(forwardFixingAdjustment);
					coefficientBridgeAtFixing   = coefficientBridgeAtFixing.mult(forwardFixingAdjustment);
				}
				if(!paymentDateOnLiborStart) {
					coefficientBridgeAtPayment  = coefficientBridgeAtPayment.mult(forwardPaymentAdjustment);
					paymentExpLiborOverPeriod   = paymentExpLiborOverPeriod.mult(forwardPaymentAdjustment);	
				}	
			}
			//end analytic ForwardCurve adjustment
			
			RandomVariableInterface completeVariance;
			if(fixingDateOnLiborStart) {
				completeVariance = coefficientBridgeAtPayment.pow(2.0).mult(bridgeVarianceAtPayment);
			} else if(paymentDateOnLiborStart) {
				completeVariance = coefficientBridgeAtFixing.pow(2.0).mult(bridgeVarianceAtFixing);		
			} else {
				completeVariance = coefficientBridgeAtFixing.pow(2.0).mult(bridgeVarianceAtFixing)
						.add(coefficientBridgeAtPayment.pow(2.0).mult(bridgeVarianceAtPayment));
				if(samePeriod) completeVariance = completeVariance.sub(coefficientBridgeAtFixing.mult(coefficientBridgeAtPayment).mult(2*bridgeCovariance));
			}
			
			RandomVariableInterface residuelRandomVariable = fixingExpLiborOverPeriod.div( numeraireAtEndFixing.mult(paymentDate - fixingDate) )
					.sub( paymentExpLiborOverPeriod.mult( (1 + strike * (paymentDate - fixingDate)) )
					.div(numeraireAtEndPayment.mult(paymentDate - fixingDate)) );

			//end Calculate needed Variables
			
			
			//Calculate Average
			/*temp*/ long start = System.currentTimeMillis();
			
			double[] doubleValues = new double[residuelRandomVariable.size()];
			for (int pathIndex = 0; pathIndex < doubleValues.length; pathIndex++) {
				double residuelVariable = residuelRandomVariable.get(pathIndex);
				double varianceVariable = completeVariance.get(pathIndex);
				//TODO: Tritt fall ein?
				if(varianceVariable>0) {
					double erf = Erf.erf(residuelVariable / Math.sqrt(2 * varianceVariable)); 
					double lt  = 0.5 * residuelVariable * (erf - 1) 
								+ Math.sqrt(varianceVariable / (2 * Math.PI)) * Math.exp(- residuelVariable * residuelVariable / (2 * varianceVariable)); 
					doubleValues[pathIndex] = lt + residuelVariable;
				} else {
					doubleValues[pathIndex] = ((- residuelVariable) < 0) ?  residuelVariable : 0.0;
				}
			}
			values = new RandomVariable(evaluationTime, doubleValues);
			/*temp*/ System.out.println("Time for cumulative Distribution Part in Analytic Formula: " + (System.currentTimeMillis() - start) + "ms");
			
			
			//temp test 7 (bb)
			RandomVariableInterface test7 = (interpolationModel.getBrownianBridge(liborStartIndexFixing, fixingDate).mult(coefficientBridgeAtFixing))
					   .sub(interpolationModel.getBrownianBridge(liborStartIndexPayment, paymentDate).mult(coefficientBridgeAtPayment)).add(residuelRandomVariable);
			System.out.println("test7: " + test7.floor(0.0).getAverage());
			
			System.out.println( "BB: " + interpolationModel.getBrownianBridge(liborStartIndexFixing, fixingDate).mult(coefficientBridgeAtFixing)
					   .sub(interpolationModel.getBrownianBridge(liborStartIndexPayment, paymentDate).mult(coefficientBridgeAtPayment)).getVariance());
			System.out.println("comp: " + completeVariance);
			
			double[] test8Array = new double[residuelRandomVariable.size()];
			double test8 = 0.0;
			for (int i = 0; i < test8Array.length; i++) {
				test8Array[i] = interpolationModel.getBrownianBridge(liborStartIndexFixing, fixingDate).mult(coefficientBridgeAtFixing.mult(numeraireAtEndFixing))
						   .sub(interpolationModel.getBrownianBridge(liborStartIndexPayment, paymentDate).mult(coefficientBridgeAtPayment.mult(numeraireAtEndFixing)))
								   .add(numeraireAtEndFixing.mult(residuelRandomVariable).get(i)).floor(0.0).getAverage();
				test8 += test8Array[i] / numeraireAtEndFixing.get(i);
			}
			System.out.println("test8: " + test8 / test8Array.length);
			System.out.println(interpolationModel.getBrownianBridge(liborStartIndexFixing, fixingDate).mult(coefficientBridgeAtFixing)
					   .sub(interpolationModel.getBrownianBridge(liborStartIndexPayment, paymentDate).mult(coefficientBridgeAtPayment))
							   .add(residuelRandomVariable).floor(0.0).getAverage());
			/*
			//temp test 6 (residuel)
			RandomVariableInterface liborT1 = interpolationModel.getLIBOR(fixingDate, fixingDate, liborEndTimeFixing).mult(liborEndTimePayment - fixingDate).add(1.0).sub(interpolationModel.getBrownianBridge(liborStartIndexFixing, fixingDate));
			RandomVariableInterface liborT2 = interpolationModel.getLIBOR(paymentDate, paymentDate, liborEndTimeFixing).mult(liborEndTimePayment - paymentDate).add(1.0).sub(interpolationModel.getBrownianBridge(liborStartIndexFixing, paymentDate));
			RandomVariableInterface residuel2 = liborT1.sub( liborT2.mult(1+strike*(paymentDate - fixingDate))  ).div(interpolationModel.getNumeraire(liborEndTimeFixing).mult(paymentDate - fixingDate));
			System.out.println("Test6: " + residuel2.getAverage() +"\t"+ residuelRandomVariable.getAverage());
			*/
			
			
			//temp test 5:
			/*RandomVariableInterface coefficientBridgeAtFixing2 = interpolationModel.getNumeraire(liborEndTimeFixing).mult(paymentDate - fixingDate).invert(); 
			RandomVariableInterface coefficientBridgeAtPayment2 = interpolationModel.getNumeraire(liborEndTimeFixing).mult(paymentDate - fixingDate).div(1 + strike * (paymentDate - fixingDate)).invert(); 
			double test5 = (interpolationModel.getBrownianBridge(liborStartIndexFixing, fixingDate).mult(coefficientBridgeAtFixing2.get(5)))
							.sub(interpolationModel.getBrownianBridge(liborStartIndexPayment, paymentDate).mult(coefficientBridgeAtPayment2.get(5))).getVariance();
			//RandomVariableInterface test = residuelRandomVariable.add(coefficientBridgeAtFixing.mult(interpolationModel.getBrownianBridge(liborStartIndex, fixingDate)))
			//								.sub(coefficientBridgeAtPayment.mult(interpolationModel.getBrownianBridge(liborStartIndex, paymentDate)));
			//RandomVariableInterface test = adjustedNumeraireAtEnd.mult(paymentDate-fixingDate).invert().mult(model.getLIBOR(fixingDate, fixingDate, liborEndTime).mult(liborEndTime-fixingDate).add(1)
			//		.sub(model.getLIBOR(paymentDate, paymentDate, liborEndTime).mult(liborEndTime-paymentDate).add(1).mult(1+strike*(paymentDate-fixingDate)) ) );
			System.out.println("test5: " + (Math.abs(completeVariance.get(5) - test5) / completeVariance.get(5))*100 + "%" );
			System.out.println();*/
			//end temp test
			
			
		} else {
			RandomVariableInterface numeraireAtFixing  	  				= model.getNumeraire(fixingDate);
			RandomVariableInterface numeraireAtPayment	 				= model.getNumeraire(paymentDate);
			RandomVariableInterface monteCarloProbabilitiesAtPayment    = model.getMonteCarloWeights(paymentDate);
			double periodLenght 					    			  	= paymentDate - fixingDate;
			System.out.println("nAF:" + numeraireAtPayment);
			values = numeraireAtPayment.div(numeraireAtFixing).sub(1).div(periodLenght).sub(strike).floor(0.0);
			values = values.div(numeraireAtPayment).mult(monteCarloProbabilitiesAtPayment);
		}
		RandomVariableInterface numeraireAtEvaluation				= model.getNumeraire(evaluationTime);
		RandomVariableInterface monteCarloProbabilitiesAtEvaluation = model.getMonteCarloWeights(evaluationTime);
		if(useAnalyticFormula) monteCarloProbabilitiesAtEvaluation  = model.getRandomVariableForConstant(1.0);
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
