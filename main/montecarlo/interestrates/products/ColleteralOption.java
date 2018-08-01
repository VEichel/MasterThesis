package montecarlo.interestrates.products;

//TODO: for fixingOnLiborStart use change liborStartToEnd and dont use any libor. Numeraire is enough!

import java.util.Arrays;

import org.apache.commons.math3.analysis.function.Abs;
import org.apache.commons.math3.genetics.GeneticAlgorithm;
import org.apache.commons.math3.special.Erf;

import montecarlo.interestrates.LIBORMarketModelWithBridge;
import montecarlo.interestrates.modelplugins.AbstractLiborCovarianceModelWithInterpolation;
import net.finmath.exception.CalculationException;
import net.finmath.functions.AnalyticFormulas;
import net.finmath.functions.NormalDistribution;
import net.finmath.marketdata.model.AnalyticModelInterface;
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
			LIBORMarketModelWithBridge interpolationModel = (LIBORMarketModelWithBridge) model.getModel();
			
			int fixingIndex 	  			 = interpolationModel.getTimeIndex(fixingDate);
			int paymentIndex 	 			 = interpolationModel.getTimeIndex(paymentDate);
			int evaluationIndex   		 	 = interpolationModel.getTimeIndex(evaluationTime);
			
			int liborStartIndexFixing  			 = interpolationModel.getLiborPeriodIndex(fixingDate);
			int liborStartIndexPayment 			 = interpolationModel.getLiborPeriodIndex(paymentDate);
			
			double periodLenght = paymentDate - fixingDate;
			
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
			

		    double longPeriodLenght = liborEndTimeFixing - liborStartTimeFixing;
		    
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
			
			/*RandomVariableInterface coefficientBridgeAtFixing  = numeraireAtEndFixing.mult(paymentDate - fixingDate).invert();
			RandomVariableInterface coefficientBridgeAtPayment = numeraireAtEndPayment.mult(paymentDate - fixingDate).div(1 + strike * (paymentDate - fixingDate)).invert();
			*/
			double coefficientBridgeAtFixing  =  1 / (paymentDate - fixingDate);
			double coefficientBridgeAtPayment =  (1 + strike * (paymentDate - fixingDate)) / (paymentDate - fixingDate);
			
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
					coefficientBridgeAtFixing   = coefficientBridgeAtFixing * forwardFixingAdjustment;
				}
				if(!paymentDateOnLiborStart) {
					coefficientBridgeAtPayment  = coefficientBridgeAtPayment * forwardPaymentAdjustment;
					paymentExpLiborOverPeriod   = paymentExpLiborOverPeriod.mult(forwardPaymentAdjustment);	
				}	
			}
			//end analytic ForwardCurve adjustment
			
			double completeVariance;
			if(fixingDateOnLiborStart) {
				completeVariance = coefficientBridgeAtPayment * coefficientBridgeAtPayment * bridgeVarianceAtPayment;
			} else if(paymentDateOnLiborStart) {
				completeVariance = coefficientBridgeAtFixing * coefficientBridgeAtFixing * bridgeVarianceAtFixing;		
			} else {
				completeVariance = coefficientBridgeAtFixing * coefficientBridgeAtFixing * bridgeVarianceAtFixing
						 + coefficientBridgeAtPayment * coefficientBridgeAtPayment * bridgeVarianceAtPayment ;
				if(samePeriod) completeVariance = completeVariance - coefficientBridgeAtFixing * coefficientBridgeAtPayment * 2 * bridgeCovariance ;
			}
			
			RandomVariableInterface residuelRandomVariable = fixingExpLiborOverPeriod
					.sub( paymentExpLiborOverPeriod.mult( (1 + strike * (paymentDate - fixingDate))) )
					.div(paymentDate - fixingDate);
			
			//end Calculate needed Variables
			
			
			//Calculate Average
			/*temp*/ long start = System.currentTimeMillis();
			
			/*
			double[] doubleValues = new double[residuelRandomVariable.size()];
			for (int pathIndex = 0; pathIndex < doubleValues.length; pathIndex++) {
				double residuelVariable = residuelRandomVariable.get(pathIndex);
				double varianceVariable = completeVariance;
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
			values = new RandomVariable(evaluationTime, doubleValues).div(numeraireAtEndFixing);
			/*temp*/ System.out.println("Time for cumulative Distribution Part in Analytic Formula: " + (System.currentTimeMillis() - start) + "ms");
			
			/*
			//temp test 7 (bb)
			RandomVariableInterface test7 = (interpolationModel.getBrownianBridge(liborStartIndexFixing, fixingDate).mult(coefficientBridgeAtFixing))
					   .sub(interpolationModel.getBrownianBridge(liborStartIndexPayment, paymentDate).mult(coefficientBridgeAtPayment)).add(residuelRandomVariable);
			System.out.println("test7: " + test7.floor(0.0).getAverage());
			
			System.out.println( "BB: " + interpolationModel.getBrownianBridge(liborStartIndexFixing, fixingDate).mult(coefficientBridgeAtFixing)
					   .sub(interpolationModel.getBrownianBridge(liborStartIndexPayment, paymentDate).mult(coefficientBridgeAtPayment)).getVariance());
			System.out.println("comp: " + completeVariance);
			*/
			/*
			RandomVariableInterface numeraireAtStart = model.getNumeraire(liborStartTimeFixing);
			RandomVariableInterface numeraireAtEnd 	 = model.getNumeraire(liborEndTimeFixing);
			double[] test8Array = new double[residuelRandomVariable.size()];
			double test8 = 0.0;
			double alpha1 = (liborEndTimeFixing - fixingDate) / (liborEndTimeFixing - liborStartTimeFixing);
			double alpha2 = (liborEndTimeFixing - paymentDate) / (liborEndTimeFixing - liborStartTimeFixing);
			RandomVariableInterface coefficient1Rv = numeraireAtStart.log().mult(- alpha1).sub( numeraireAtEnd.log().mult(1 - alpha1)).exp().div(paymentDate - fixingDate) ;
			RandomVariableInterface coefficient2Rv = numeraireAtStart.log().mult(- alpha2).sub( numeraireAtEnd.log().mult(1 - alpha2)).exp()
									.div(paymentDate - fixingDate).mult(1 + strike * (paymentDate - fixingDate)) ;
			System.out.println("coeff1: " + coefficient1Rv.invert());
			System.out.println("coeff2: " + coefficient2Rv.invert());
			
			for (int i = 0; i < test8Array.length; i++) {
				double coefficient1 = coefficient1Rv.get(i);
				double coefficient2 = coefficient2Rv.get(i);
				
				test8Array[i] = interpolationModel.getBrownianBridge(liborStartIndexFixing, fixingDate).mult(-1.0).exp().mult(coefficient1)
						   .sub(interpolationModel.getBrownianBridge(liborStartIndexPayment, paymentDate).mult(-1.0).exp().mult(coefficient2))
								   .floor(0.0).getAverage();
				test8 += test8Array[i];
			}
			System.out.println("test8: " + test8 / test8Array.length);
			System.out.println( interpolationModel.getBrownianBridge(liborStartIndexFixing, fixingDate).mult(-1.0).exp().mult(coefficient1Rv)
						   .sub(interpolationModel.getBrownianBridge(liborStartIndexPayment, paymentDate).mult(-1.0).exp().mult(coefficient2Rv))
								   .floor(0.0).getAverage() );
			*/
			/*System.out.println(values.getAverage());
			System.out.println(interpolationModel.getBrownianBridge(liborStartIndexFixing, fixingDate).mult(coefficientBridgeAtFixing)
					   .sub(interpolationModel.getBrownianBridge(liborStartIndexPayment, paymentDate).mult(coefficientBridgeAtPayment))
							   .add(residuelRandomVariable).floor(0.0).div(numeraireAtEndFixing).getAverage());
							   */
			/*values = interpolationModel.getBrownianBridge(liborStartIndexFixing, fixingDate).mult(coefficientBridgeAtFixing)
					   .sub(interpolationModel.getBrownianBridge(liborStartIndexPayment, paymentDate).mult(coefficientBridgeAtPayment))
					   .add(residuelRandomVariable).floor(0.0).div(numeraireAtEndFixing);*/
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
			
			
			/**
			 * NEU
			 */
			double alphaFixing = (liborEndTimeFixing - fixingDate) / (liborEndTimeFixing - liborStartTimeFixing);
			double alphaPayment = (liborEndTimeFixing - paymentDate) / (liborEndTimeFixing - liborStartTimeFixing);
			/**
			 * Ajustment for Numeraire Change:
			 */
			double discountFactorFixing = interpolationModel.getDiscountCurve().getDiscountFactor(interpolationModel.getAnalyticModel(), fixingDate);
			double discountFactorPayment = interpolationModel.getDiscountCurve().getDiscountFactor(interpolationModel.getAnalyticModel(), paymentDate);
			double discountFactorLiborStart = interpolationModel.getDiscountCurve().getDiscountFactor(interpolationModel.getAnalyticModel(), liborStartTimeFixing);
			double discountFactorLiborEnd = interpolationModel.getDiscountCurve().getDiscountFactor(interpolationModel.getAnalyticModel(), liborEndTimeFixing);
			double integratedVarianceOfLibor = interpolationModel.getIntegratedLIBORCovariance()[model.getTimeIndex(liborStartTimeFixing)][liborStartIndexFixing][liborStartIndexFixing];
			System.out.println("integratedVarianceOfLibor: " + integratedVarianceOfLibor);
			/*double numeraireChangeAdjustment = (discountFactorFixing / (discountFactorPayment * periodLenght) - 1) * ( - 1 +
					Math.exp(-(alphaFixing - alphaPayment) * alphaPayment * integratedVarianceOfLibor * Math.pow(1 - discountFactorLiborEnd/discountFactorLiborStart, 2) ) );*/
			/**
			 * END Adjustment
			 */
			
			
			//double completeAlpha = (fixingDate - paymentDate) / (liborEndTimeFixing - liborStartTimeFixing);

			//RandomVariableInterface completeNumeraire = model.getNumeraire(liborStartTimeFixing).log().mult(completeAlpha).add(model.getNumeraire(liborEndTimeFixing).log().mult(- completeAlpha)).exp();
			RandomVariableInterface interpolationNumeraireFixing  = model.getNumeraire(liborStartTimeFixing).log().mult(alphaFixing).add(model.getNumeraire(liborEndTimeFixing).log().mult(1 - alphaFixing)).exp();
			RandomVariableInterface interpolationNumerairePayment = model.getNumeraire(liborStartTimeFixing).log().mult(alphaPayment).add(model.getNumeraire(liborEndTimeFixing).log().mult(1 - alphaPayment)).exp();
			RandomVariableInterface interpolationNumeraire		  = interpolationNumerairePayment.div(interpolationNumeraireFixing);
	
			//test
			RandomVariableInterface brownianBridge1 = interpolationModel.getBrownianBridge(liborStartIndexFixing, fixingDate);
			RandomVariableInterface brownianBridge2 = interpolationModel.getBrownianBridge(liborStartIndexFixing, paymentDate);
			//test END
			
			if(interpolationModel.getDiscountCurve() != null) {
				double numeraireAdjustment   = (interpolationNumerairePayment.invert().getAverage() / interpolationNumeraireFixing.invert().getAverage()) * Math.exp(0.5 * (bridgeVarianceAtPayment - bridgeVarianceAtFixing)) 
												* discountFactorFixing / discountFactorPayment;
				interpolationNumeraire = interpolationNumeraire.mult(numeraireAdjustment);
			}
			
			
			
			
			double numeraire1Moment = interpolationNumeraire.getAverage();
			double numeraire2Moment = interpolationNumeraire.pow(2).getAverage();
			
			double mu           = 2*Math.log(numeraire1Moment) - 0.5*Math.log(numeraire2Moment);
			double sigmaSquared = -2*Math.log(numeraire1Moment) + completeVariance / 2 + Math.log(numeraire2Moment);
			
			
			
			//test START
			RandomVariableInterface brownianBridge = interpolationNumeraire.mult(brownianBridge2.sub(brownianBridge1).exp()).sub(1.0).div(periodLenght); 
			System.out.println(brownianBridge.getAverage() +"\t" + (Math.exp(mu+sigmaSquared/2)) );
			System.out.println(brownianBridge.pow(2).getAverage() +"\t" + (Math.exp(2*mu+2*sigmaSquared)));
			//test END
			
			double modifiedStrike = 1.0 + strike * periodLenght;
			double value = 0.5 * Math.exp(mu + sigmaSquared / 2) * (1 - Erf.erf((Math.log(modifiedStrike) - mu - sigmaSquared) / Math.sqrt(2 * sigmaSquared))) + 0.5 * modifiedStrike * ( -1 + Erf.erf((Math.log(modifiedStrike) - mu) / Math.sqrt(2 * sigmaSquared)) );
			
			value = discountFactorPayment / (paymentDate - fixingDate) * value;
			
			
			/**
			 * new 17.07
			 */
			double varianceOfBB = 0.0;
			double varianceOfBBFixing = 0.0;
			double varianceOfBBPayment = 0.0;
			double varianceOfBBSumPart = 0.0;
			for (int i = interpolationModel.getTimeIndex(liborStartTimeFixing); i < fixingIndex; i++) {
				double varianceI     = covarianceModel.getVarianceForInterpolation(i).getAverage();
				double ti     = interpolationModel.getTime(i);
				double ti2    = interpolationModel.getTime(i+1);
				varianceOfBBSumPart += varianceI * (ti2 -ti) / ((liborEndTimeFixing-ti2)*(liborEndTimeFixing-ti));
			}
			varianceOfBBFixing = varianceOfBBSumPart * Math.pow((liborEndTimeFixing - fixingDate) , 2);
			varianceOfBB       = varianceOfBBFixing + varianceOfBBSumPart *(- 2.0*(liborEndTimeFixing - fixingDate) * (liborEndTimeFixing - paymentDate));	
			for (int i = fixingIndex; i < paymentIndex; i++) {
				double varianceI    = covarianceModel.getVarianceForInterpolation(i).getAverage();
				double ti     = interpolationModel.getTime(i);
				double ti2    = interpolationModel.getTime(i+1);
				varianceOfBBSumPart += varianceI * (ti2 -ti) / ((liborEndTimeFixing-ti2)*(liborEndTimeFixing-ti));
			}
			varianceOfBBPayment =  varianceOfBBSumPart * Math.pow((liborEndTimeFixing - paymentDate) , 2);
			varianceOfBB += varianceOfBBPayment;
			System.out.println("variances: " + varianceOfBBFixing +"\t" + varianceOfBBPayment + "\t" + varianceOfBB);
			
			System.out.println("31.07:" + brownianBridge1.sub(brownianBridge2).getVariance() +"\t" + varianceOfBB);
			
			System.out.println("bs: " + AnalyticFormulas.blackScholesGeneralizedOptionValue( discountFactorFixing / discountFactorPayment  ,
									periodLenght * Math.sqrt(integratedVarianceOfLibor/liborStartTimeFixing) * (1.0 / longPeriodLenght)* (1.0 - discountFactorLiborEnd/discountFactorLiborStart)
									, liborStartTimeFixing
									, modifiedStrike
									, discountFactorPayment / periodLenght));
			values = new RandomVariable(0.0);
			
			double adjustment = 1.0;
			//For the adjustments:
			//if(interpolationModel.getForwardRateCurve() != null) {
				double forwardCurveFixingEnd    = 1 + interpolationModel.getForwardRateCurve().getForward(interpolationModel.getAnalyticModel(), fixingDate, liborEndTimeFixing - fixingDate) * (liborEndTimeFixing - fixingDate);
				double forwardCurvePaymentEnd   = 1 + interpolationModel.getForwardRateCurve().getForward(interpolationModel.getAnalyticModel(), paymentDate, liborEndTimePayment - paymentDate) * (liborEndTimePayment - paymentDate);
				double forwardCurveLongFixing   = 1 + interpolationModel.getForwardRateCurve().getForward(interpolationModel.getAnalyticModel(), liborStartTimeFixing, liborEndTimeFixing - liborStartTimeFixing) * (liborEndTimeFixing - liborStartTimeFixing);
				double forwardCurveLongPayment  = samePeriod ? forwardCurveLongFixing : 
													1 + interpolationModel.getForwardRateCurve().getForward(interpolationModel.getAnalyticModel(), liborStartTimePayment, liborEndTimePayment - liborStartTimePayment) * (liborEndTimePayment - liborStartTimePayment); 
				
				double forwardCurveFixingEndLog  = Math.exp(Math.log(forwardCurveLongFixing)*(liborEndTimeFixing - fixingDate)/(liborEndTimeFixing - liborStartTimeFixing));
				double forwardCurvePaymentEndLog = Math.exp(Math.log(forwardCurveLongPayment)*(liborEndTimePayment - paymentDate)/(liborEndTimePayment - liborStartTimePayment)); 
				
				double forwardFixingAdjustment   = forwardCurveFixingEnd  / forwardCurveFixingEndLog;
				System.out.println(forwardFixingAdjustment);
				double forwardPaymentAdjustment  = forwardCurvePaymentEnd / forwardCurvePaymentEndLog;
				System.out.println(forwardPaymentAdjustment);
				
				adjustment =  forwardFixingAdjustment / forwardPaymentAdjustment ;
			//}
			System.out.println(adjustment);
			if(interpolationModel.getDiscountCurve() != null) {
				double numeraireAtPeriodEndAdjustment   = interpolationModel.getUnAdjustedNumeraire(liborEndTimeFixing).invert().getAverage()   / discountFactorLiborEnd;
				double numeraireAtPeriodStartAdjustment = interpolationModel.getUnAdjustedNumeraire(liborStartTimeFixing).invert().getAverage() / discountFactorLiborStart;		
				System.out.println("numeraireAtPeriodEndAdjustment " + numeraireAtPeriodEndAdjustment);
				System.out.println("numeraireAtPeriodStartAdjustment " + numeraireAtPeriodStartAdjustment);
				double numeraireAtFixingAdjustment   = interpolationModel.getUnAdjustedNumeraire(liborEndTimeFixing).log().mult(alphaFixing - 1)
														.add(interpolationModel.getUnAdjustedNumeraire(liborStartTimeFixing).log().mult(-alphaFixing)).exp().getAverage() 
														* Math.exp(varianceOfBBFixing / 2.0) * forwardFixingAdjustment/ discountFactorFixing;
				System.out.println("numeraireAtFixingAdjustment " + numeraireAtFixingAdjustment);
				double numeraireAtPaymentAdjustment = interpolationModel.getUnAdjustedNumeraire(liborEndTimeFixing).pow(alphaPayment - 1)
														.mult(interpolationModel.getUnAdjustedNumeraire(liborStartTimeFixing).pow(-alphaPayment)).getAverage() 
														* Math.exp(varianceOfBBPayment / 2.0) * forwardPaymentAdjustment/ discountFactorPayment;
				System.out.println("numeraireAtPaymentAdjustment " + numeraireAtPaymentAdjustment);
				//adjustment *= Math.pow(numeraireAtPeriodEndAdjustment, alphaFixing - alphaPayment) * Math.pow(numeraireAtPeriodStartAdjustment, alphaPayment - alphaFixing)
				//		   * numeraireAtPaymentAdjustment / numeraireAtFixingAdjustment;
				
				adjustment *= numeraireAtPaymentAdjustment / numeraireAtFixingAdjustment;
			}
			
			
			System.out.println("01.08 adjusmtent: " + adjustment);

			//test bs mit exp(BB) am ende drauf: 18.07
			System.out.println("31.07 2: " + analyticFormula( /*adjustment* */ discountFactorFixing / discountFactorPayment  ,
									periodLenght * Math.sqrt(integratedVarianceOfLibor/fixingDate) * (1.0 / longPeriodLenght)* (1.0 - discountFactorLiborEnd/discountFactorLiborStart)
									, liborStartTimeFixing
									, modifiedStrike
									, discountFactorPayment / periodLenght, Math.sqrt(varianceOfBB)));
			System.out.println();
			RandomVariableInterface bb = brownianBridge2.sub(brownianBridge1).exp();
			System.out.println(discountFactorFixing / discountFactorPayment);
			double test0108 =0.0;
			 for (int i = 0; i < bb.size(); i++) {
				 test0108 += AnalyticFormulas.blackScholesGeneralizedOptionValue(bb.get(i) * discountFactorFixing / discountFactorPayment  ,
							periodLenght * Math.sqrt(integratedVarianceOfLibor/liborStartTimeFixing) * (1.0 / longPeriodLenght)* (1.0 - discountFactorLiborEnd/discountFactorLiborStart)
							, liborStartTimeFixing
							, modifiedStrike
							, discountFactorPayment / periodLenght);
			}
			 System.out.println("test0108: " + (test0108 / bb.size()) );
			 System.out.println();
			//end test
		} else {
			
			
			RandomVariableInterface numeraireAtFixing  	  				= model.getNumeraire(fixingDate);
			RandomVariableInterface numeraireAtPayment	 				= model.getNumeraire(paymentDate);
			RandomVariableInterface monteCarloProbabilitiesAtPayment    = model.getMonteCarloWeights(paymentDate);
			double periodLenght 					    			  	= paymentDate - fixingDate;
			values = numeraireAtPayment.div(numeraireAtFixing).sub(1.0).div(periodLenght).sub(strike).floor(0.0);
			values = values.div(numeraireAtPayment).mult(monteCarloProbabilitiesAtPayment);
		}
		RandomVariableInterface numeraireAtEvaluation				= model.getNumeraire(evaluationTime);
		RandomVariableInterface monteCarloProbabilitiesAtEvaluation = model.getMonteCarloWeights(evaluationTime);
		values = values.mult(numeraireAtEvaluation).div(monteCarloProbabilitiesAtEvaluation);
		//double discountFactorPayment = model.getModel().getDiscountCurve().getDiscountFactor(model.getModel().getAnalyticModel(), paymentDate);
		
		//values = values.mult(discountFactorPayment);
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
	
	
	public static double analyticFormula(double forward, double volatility, double optionMaturity, double optionStrike, double payoffUnit, double bridgeVolatility) {
		
		//double dPlus = (Math.log(forward / optionStrike) + 0.5 * volatility * volatility * optionMaturity) / (Math.sqrt(optionMaturity));
		//double dMinus = dPlus - volatility * Math.sqrt(optionMaturity);
		
		double phiDenominator = Math.sqrt(2*(bridgeVolatility*bridgeVolatility + volatility*volatility*optionMaturity));
		double phiSInput = (Math.log(forward / optionStrike) + 0.5 * volatility * volatility * optionMaturity + bridgeVolatility * bridgeVolatility) / phiDenominator;
		double phiKInput = (Math.log(forward / optionStrike) - 0.5 * volatility * volatility * optionMaturity) / phiDenominator;
		
		double value = forward * Math.exp(bridgeVolatility*bridgeVolatility*0.5) * NormalDistribution.cumulativeDistribution(phiSInput) 
						- optionStrike  * NormalDistribution.cumulativeDistribution(phiKInput);
		return value * payoffUnit;
	}
}
