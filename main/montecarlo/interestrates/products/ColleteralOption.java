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
			
			boolean fixingDateOnLiborStart = true;
			if(liborStartIndexFixing<0) {
				liborStartIndexFixing = (-liborStartIndexFixing - 1) - 1;
				fixingDateOnLiborStart = false;
			}
		    int liborEndIndexFixing	 		 = liborStartIndexFixing + 1;		    
		    double liborStartTimeFixing		 = interpolationModel.getLiborPeriod(liborStartIndexFixing);
		    double liborEndTimeFixing	 	 = interpolationModel.getLiborPeriod(liborEndIndexFixing);
		    
		    boolean paymentDateOnLiborStart = true;
			if(liborStartIndexPayment<0) {
				liborStartIndexPayment = (-liborStartIndexPayment - 1) - 1;
				paymentDateOnLiborStart = false;
			}
		    int liborEndIndexPayment 		 = liborStartIndexPayment + 1;		    
		    double liborStartTimePayment	 = interpolationModel.getLiborPeriod(liborStartIndexPayment);
		    double liborEndTimePayment	 	 = interpolationModel.getLiborPeriod(liborEndIndexPayment);
		    
		    System.out.println(paymentDateOnLiborStart);
		    System.out.println(fixingDateOnLiborStart);
		    
		    if(fixingDateOnLiborStart && paymentDateOnLiborStart) {
		    	this.useAnalyticFormula = false;
		    	return this.getValue(evaluationTime, model);
		    }
		    
		    int periodStartToPaymentIndex    = paymentIndex - interpolationModel.getTimeIndex(liborStartTimeFixing);
		    int periodStartToFixingIndex     = fixingIndex  - interpolationModel.getTimeIndex(liborStartTimeFixing);
		    
			//Calculate BB Variance:
			AbstractLiborCovarianceModelWithInterpolation covarianceModel = (AbstractLiborCovarianceModelWithInterpolation) interpolationModel.getCovarianceModel();
			double bridgeVarianceAtFixing  = 0.0;
			double bridgeVarianceAtPayment = 0.0;
			double bridgeCovariance        = 0.0;
			for(int i = 1; i < periodStartToPaymentIndex + 1 ; i++) {
				double varianceI    = covarianceModel.getVarianceForInterpolation(i-1).getAverage();
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
			//Finished BB Variance calculations
			
			//Calculate needed Variables (complete covariance & residuelRandomVariable)
			RandomVariableInterface	adjustedNumeraireAtEnd     = interpolationModel.getNumeraire(liborEndTimeFixing);
			
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
			
			
			RandomVariableInterface logLiborOverPeriod         = interpolationModel.getLIBOR(interpolationModel.getTimeIndex(interpolationModel.getLiborPeriod(liborStartIndexFixing)), liborStartIndexFixing)
																	.mult(liborEndTimeFixing - liborStartTimeFixing).add(1.0).log();
			RandomVariableInterface fixingExpLiborOverPeriod   = logLiborOverPeriod.mult((liborEndTimeFixing-fixingDate)/(liborEndTimeFixing-liborStartTimeFixing)).exp();
			RandomVariableInterface paymentExpLiborOverPeriod  = logLiborOverPeriod.mult((liborEndTimeFixing-paymentDate)/(liborEndTimeFixing-liborStartTimeFixing)).exp();

			RandomVariableInterface coefficientBridgeAtFixing  = adjustedNumeraireAtEnd.mult((paymentDate - fixingDate) * numeraireAdjustmentForFixing).invert();
			RandomVariableInterface coefficientBridgeAtPayment = adjustedNumeraireAtEnd.mult((paymentDate - fixingDate) * numeraireAdjustmentForPayment).div(1 + strike * (paymentDate - fixingDate)).invert();
			
			
			
			//analytic ForwardCurve adjustment
			if(interpolationModel.getForwardRateCurve() != null) {
				double forwardCurveFixingEnd    = 1 + interpolationModel.getForwardRateCurve().getForward(interpolationModel.getAnalyticModel(), fixingDate, liborEndTimeFixing - fixingDate) * (liborEndTimeFixing - fixingDate);
				double forwardCurvePaymentEnd   = 1 + interpolationModel.getForwardRateCurve().getForward(interpolationModel.getAnalyticModel(), paymentDate, liborEndTimeFixing - paymentDate) * (liborEndTimeFixing - paymentDate);
				double forwardCurveLongPeriod   = 1 + interpolationModel.getForwardRateCurve().getForward(interpolationModel.getAnalyticModel(), liborStartTimeFixing, liborEndTimeFixing - liborStartTimeFixing) * (liborEndTimeFixing - liborStartTimeFixing);
			
				double forwardCurveFixingEndLog  = Math.exp(Math.log(forwardCurveLongPeriod)*(liborEndTimeFixing - fixingDate)/(liborEndTimeFixing - liborStartTimeFixing));
				double forwardCurvePaymentEndLog = Math.exp(Math.log(forwardCurveLongPeriod)*(liborEndTimeFixing - paymentDate)/(liborEndTimeFixing - liborStartTimeFixing)); 
				
				double forwardFixingAdjustment  = forwardCurveFixingEnd  / forwardCurveFixingEndLog;
				double forwardPaymentAdjustment = forwardCurvePaymentEnd / forwardCurvePaymentEndLog;
				
				coefficientBridgeAtFixing  = coefficientBridgeAtFixing.mult(forwardFixingAdjustment);
				coefficientBridgeAtPayment = coefficientBridgeAtPayment.mult(forwardPaymentAdjustment);
				fixingExpLiborOverPeriod   = fixingExpLiborOverPeriod.mult(forwardFixingAdjustment);
				paymentExpLiborOverPeriod  = paymentExpLiborOverPeriod.mult(forwardPaymentAdjustment);
			}
			
			RandomVariableInterface completeVariance = coefficientBridgeAtFixing.pow(2.0).mult(bridgeVarianceAtFixing)
					.add(coefficientBridgeAtPayment.pow(2.0).mult(bridgeVarianceAtPayment))
					.sub(coefficientBridgeAtFixing.mult(coefficientBridgeAtPayment).mult(2*bridgeCovariance));
			
			RandomVariableInterface residuelRandomVariable     = fixingExpLiborOverPeriod.div(numeraireAdjustmentForFixing).sub( paymentExpLiborOverPeriod.mult( (1 + strike * (paymentDate - fixingDate)) / numeraireAdjustmentForPayment))
					.div(adjustedNumeraireAtEnd.mult(paymentDate - fixingDate));

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
			double test = (interpolationModel.getBrownianBridge(liborStartIndex, fixingDate).mult(coefficientBridgeAtFixing.get(1)))
							.sub(interpolationModel.getBrownianBridge(liborStartIndex, paymentDate).mult(coefficientBridgeAtPayment.get(1))).getVariance();
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
