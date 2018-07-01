package notSorted;
import java.util.HashMap;
import java.util.Map;

import montecarlo.interestrates.LiborMarketModelWithBridgeInterpolation;
import montecarlo.interestrates.modelplugins.AbstractLiborCovarianceModelWithInterpolation;
import montecarlo.interestrates.modelplugins.LiborCovarianceModelWithInterpolation;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.curves.DiscountCurveFromForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.BrownianMotionInterface;
import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORMarketModelInterface;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulation;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORCorrelationModelExponentialDecay;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORCovarianceModelFromVolatilityAndCorrelation;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORVolatilityModel;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORVolatilityModelFourParameterExponentialForm;
import net.finmath.montecarlo.interestrate.products.SimpleSwap;
import net.finmath.montecarlo.process.ProcessEulerScheme;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretization;

public class LMMWIthBB {

	final static double liborPeriodLength	= 10;
	final static double liborRateTimeHorzion	= 40.0;
	final static double lastTime	= 40.0;
	final static double dt		= 0.5;
	
	public static void main(String[] args) throws CalculationException {
		
		LIBORModelMonteCarloSimulation LMMBB = createLIBORMarketModelWithBB(1000, 3, 0.1);
		LIBORModelMonteCarloSimulation LMM = createLIBORMarketModel(1000, 3, 0.1);
		
		System.out.println(LMM.getLIBOR(10.0, 10.0, 20.0).get(3));
		System.out.println();
		
		//Test 1:

		for (int timeIndex = 20; timeIndex < 40; timeIndex++) {
			System.out.println( ((LiborMarketModelWithBridgeInterpolation) LMMBB.getModel()).getInterpolatedLibor(0, timeIndex).get(15) + "\t" +
					LMM.getLIBOR(0.0, LMM.getTime(timeIndex), 20.0).get(15));
		}
		
		
		System.out.println();
		System.out.println("Test 2");
		//Test 2:
		for (int timeIndex = 0; timeIndex < LMMBB.getTimeDiscretization().getNumberOfTimes(); timeIndex++) {
			double time = LMMBB.getTime(timeIndex);
			System.out.println( LMMBB.getNumeraire(time).get(15) + "\t" + LMM.getNumeraire(time).get(15));
		}
		System.out.println();
		/*
		System.out.println("Test 3");
		for (int timeIndex = 0; timeIndex < LMMBB.getTimeDiscretization().getNumberOfTimes(); timeIndex++) {
			double time = LMMBB.getTime(timeIndex);
			System.out.println( LMMBB.getNumeraire(time).get(19) + "\t" + LMM.getNumeraire(time).get(19));
		}
		System.out.println();
		/*
		for (int timeIndex = 0; timeIndex < LMM.getTimeDiscretization().getNumberOfTimeSteps(); timeIndex++) {
			
			double[] fixingDates = {LMM.getTime(timeIndex)};
			double[] paymentDates = {LMM.getTime(timeIndex + 1)};
			double[] swaprates    = {0.01};
			
			SimpleSwap swap = new SimpleSwap(fixingDates, paymentDates, swaprates);
			System.out.println(swap.getValue(LMM) + "\t" + swap.getValue(LMMBB));
			
		}*/
		
	}
	
	public static LIBORModelMonteCarloSimulation createLIBORMarketModelWithBB(
			int numberOfPaths, int numberOfFactors, double correlationDecayParam) throws CalculationException {

		/*
		 * Create the libor tenor structure and the initial values
		 */

		TimeDiscretization liborPeriodDiscretization = new TimeDiscretization(0.0, (int) (liborRateTimeHorzion / liborPeriodLength), liborPeriodLength);

		// Create the forward curve (initial value of the LIBOR market liborModel)
		ForwardCurve forwardCurve = ForwardCurve.createForwardCurveFromForwards(
				"forwardCurve"								/* name of the curve */,
				new double[] {0.5 , 1.0 , 2.0 , 5.0 , 40.0}	/* fixings of the forward */,
				new double[] {0.05, 0.05, 0.05, 0.05, 0.05}	/* forwards */,
				liborPeriodLength							/* tenor / period length */
				);

		/*
		 * Create a simulation time discretization
		 */

		TimeDiscretization timeDiscretization = new TimeDiscretization(0.0, (int) (lastTime / dt), dt);

		/*
		 * Create a volatility structure v[i][j] = sigma_j(t_i)
		 */
		double a = 0.5, b = 0.0, c = 0.25, d = 0.1;
		LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelFourParameterExponentialForm(timeDiscretization, liborPeriodDiscretization, a, b, c, d, false);		

		/*
		 * Create a correlation liborModel rho_{i,j} = exp(-a * abs(T_i-T_j))
		 */
		LIBORCorrelationModelExponentialDecay correlationModel = new LIBORCorrelationModelExponentialDecay(
				timeDiscretization, liborPeriodDiscretization, numberOfFactors,
				correlationDecayParam);


		/*
		 * Combine volatility liborModel and correlation liborModel to a covariance liborModel
		 */
		LIBORCovarianceModelFromVolatilityAndCorrelation liborCovarianceModel =
				new LIBORCovarianceModelFromVolatilityAndCorrelation(timeDiscretization,
						liborPeriodDiscretization, volatilityModel, correlationModel);

		// BlendedLocalVolatlityModel (future extension)
		//		AbstractLIBORCovarianceModel covarianceModel2 = new BlendedLocalVolatlityModel(covarianceModel, 0.00, false);

		// Set liborModel properties
		Map<String, String> properties = new HashMap<String, String>();

		// Choose the simulation measure
		properties.put("measure", LIBORMarketModel.Measure.SPOT.name());

		// Choose log normal liborModel
		properties.put("stateSpace", LIBORMarketModel.StateSpace.LOGNORMAL.name());

		// Empty array of calibration items - hence, liborModel will use given covariance
		LIBORMarketModel.CalibrationItem[] calibrationItems = new LIBORMarketModel.CalibrationItem[0];

		/* For piecewise constant model
		double[] interpolationParameters = new double[liborPeriodDiscretization.getNumberOfTimeSteps()];
		for (int liborIndex = 0; liborIndex < interpolationParameters.length; liborIndex++) {
			interpolationParameters[liborIndex] = 0.02 + liborIndex * 0.02;
		}
		*/
		double[] interpolationParameters = new double[1];
		interpolationParameters[0] = 0.00;
		 
		
		AbstractLiborCovarianceModelWithInterpolation covarianceModel = new LiborCovarianceModelWithInterpolation(liborCovarianceModel, interpolationParameters);
		
		
		/*
		 * Create corresponding LIBOR Market Model
		 */
		BrownianMotionInterface interpolationDriver = new BrownianMotion(timeDiscretization, 1, numberOfPaths, 412421);
		
		LIBORMarketModelInterface liborMarketModel = new LiborMarketModelWithBridgeInterpolation(liborPeriodDiscretization, null, forwardCurve, new DiscountCurveFromForwardCurve(forwardCurve), new RandomVariableFactory(), covarianceModel, calibrationItems, properties, interpolationDriver);

		
		BrownianMotionInterface brownianMotion = new net.finmath.montecarlo.BrownianMotion(timeDiscretization, numberOfFactors, numberOfPaths, 3141);

		ProcessEulerScheme process = new ProcessEulerScheme(brownianMotion, ProcessEulerScheme.Scheme.EULER);

		return new LIBORModelMonteCarloSimulation(liborMarketModel, process);
	}
	
	public static LIBORModelMonteCarloSimulation createLIBORMarketModel(
			int numberOfPaths, int numberOfFactors, double correlationDecayParam) throws CalculationException {

		/*
		 * Create the libor tenor structure and the initial values
		 */

		TimeDiscretization liborPeriodDiscretization = new TimeDiscretization(0.0, (int) (liborRateTimeHorzion / liborPeriodLength), liborPeriodLength);

		// Create the forward curve (initial value of the LIBOR market liborModel)
		ForwardCurve forwardCurve = ForwardCurve.createForwardCurveFromForwards(
				"forwardCurve"								/* name of the curve */,
				new double[] {0.5 , 1.0 , 2.0 , 5.0 , 40.0}	/* fixings of the forward */,
				new double[] {0.05, 0.05, 0.05, 0.05, 0.05}	/* forwards */,
				liborPeriodLength							/* tenor / period length */
				);

		/*
		 * Create a simulation time discretization
		 */


		TimeDiscretization timeDiscretization = new TimeDiscretization(0.0, (int) (lastTime / dt), dt);

		/*
		 * Create a volatility structure v[i][j] = sigma_j(t_i)
		 */
		double a = 0.5, b = 0.0, c = 0.25, d = 0.1;
		LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelFourParameterExponentialForm(timeDiscretization, liborPeriodDiscretization, a, b, c, d, false);		

		/*
		 * Create a correlation liborModel rho_{i,j} = exp(-a * abs(T_i-T_j))
		 */
		LIBORCorrelationModelExponentialDecay correlationModel = new LIBORCorrelationModelExponentialDecay(
				timeDiscretization, liborPeriodDiscretization, numberOfFactors,
				correlationDecayParam);


		/*
		 * Combine volatility liborModel and correlation liborModel to a covariance liborModel
		 */
		LIBORCovarianceModelFromVolatilityAndCorrelation covarianceModel =
				new LIBORCovarianceModelFromVolatilityAndCorrelation(timeDiscretization,
						liborPeriodDiscretization, volatilityModel, correlationModel);

		// BlendedLocalVolatlityModel (future extension)
		//		AbstractLIBORCovarianceModel covarianceModel2 = new BlendedLocalVolatlityModel(covarianceModel, 0.00, false);

		// Set liborModel properties
		Map<String, String> properties = new HashMap<String, String>();

		// Choose the simulation measure
		properties.put("measure", LIBORMarketModel.Measure.SPOT.name());

		// Choose log normal liborModel
		properties.put("stateSpace", LIBORMarketModel.StateSpace.LOGNORMAL.name());

		// Empty array of calibration items - hence, liborModel will use given covariance
		LIBORMarketModel.CalibrationItem[] calibrationItems = new LIBORMarketModel.CalibrationItem[0];
		
		/*
		 * Create corresponding LIBOR Market Model
		 */
		LIBORMarketModelInterface liborMarketModel = new LIBORMarketModel(liborPeriodDiscretization, null, forwardCurve, new DiscountCurveFromForwardCurve(forwardCurve), covarianceModel, calibrationItems, properties);

				
		
		BrownianMotionInterface brownianMotion = new net.finmath.montecarlo.BrownianMotion(timeDiscretization, numberOfFactors, numberOfPaths, 3141);

		ProcessEulerScheme process = new ProcessEulerScheme(brownianMotion, ProcessEulerScheme.Scheme.EULER);

		return new LIBORModelMonteCarloSimulation(liborMarketModel, process);
	}
}

