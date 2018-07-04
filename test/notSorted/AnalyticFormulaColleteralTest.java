package notSorted;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

import montecarlo.interestrates.LiborMarketModelWithBridgeInterpolation;
import montecarlo.interestrates.modelplugins.AbstractLiborCovarianceModelWithInterpolation;
import montecarlo.interestrates.modelplugins.LiborCovarianceModelWithInterpolation;
import montecarlo.interestrates.modelplugins.LiborCovarianceModelWithInterpolation.EvaluationTimeScalingScheme;
import montecarlo.interestrates.modelplugins.LiborCovarianceModelWithInterpolation.InterpolationVarianceScheme;
import montecarlo.interestrates.products.ColleteralOption;
import net.finmath.analytic.model.curves.DiscountCurve;
import net.finmath.analytic.model.curves.DiscountCurveInterface;
import net.finmath.analytic.model.curves.ForwardCurveInterface;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.AnalyticModelInterface;
import net.finmath.marketdata.model.curves.DiscountCurveFromForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.BrownianMotionInterface;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORMarketModelCalibrationTest;
import net.finmath.montecarlo.interestrate.LIBORMarketModelInterface;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulation;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORCorrelationModelExponentialDecay;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORCovarianceModelFromVolatilityAndCorrelation;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORVolatilityModel;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORVolatilityModelFourParameterExponentialForm;
import net.finmath.montecarlo.interestrate.products.components.Numeraire;
import net.finmath.montecarlo.process.ProcessEulerScheme;
import net.finmath.optimizer.SolverException;
import net.finmath.time.TimeDiscretization;

public class AnalyticFormulaColleteralTest {

	static int interpolationSeed = 13425453;
	
	final static double liborPeriodLength	= 10;
	final static double liborRateTimeHorzion	= 40.0;
	final static double lastTime	= 40.0;
	final static double dt		= 0.1;
	
	public static void main(String[] args) throws CalculationException, InterruptedException, SolverException {
		
		//Thread.sleep(20000);
		
		System.out.println("Test1:");
		System.out.println("------");
		test1();
		System.out.println("-------");
		System.out.println();
		System.exit(0);
	}
	
	public static void test1() throws CalculationException, SolverException {
		
		LIBORModelMonteCarloSimulation LMMBB = createLIBORMarketModelWithBB(10000, 3, 0.1);
		LIBORModelMonteCarloSimulation LMM   = createLIBORMarketModel      (10000, 3, 0.1);
		
		//For precalculating the non interpolated Libors
		LMMBB.getNumeraire(LMMBB.getTime(LMMBB.getTimeDiscretization().getNumberOfTimeSteps()));
		
		double evaluationTime = 0.0;
		double fixingDate     = 13.0;
		double paymentDate    = 14.0;
		double strike         = 0.2;
		
		//Non Analytic:
		long beforeNonAnalyticMillis = System.currentTimeMillis();
		ColleteralOption colleteralOptionNonAnalytic = new ColleteralOption(fixingDate, paymentDate, strike, false);
		double colleteralNonAnalyticPrice = colleteralOptionNonAnalytic.getValue(evaluationTime, LMMBB).getAverage();
		long afterNonAnalyticMillis	 = System.currentTimeMillis();
		
		//Analytic
		long beforeAnalyticMillis = System.currentTimeMillis();
		ColleteralOption colleteralOptionAnalytic	 = new ColleteralOption(fixingDate, paymentDate, strike, true);
		double colleteralAnalyticPrice = colleteralOptionAnalytic.getValue(evaluationTime, LMMBB).getAverage();
		long afterAnalyticMillis = System.currentTimeMillis();
		
		
		//OLD Analytic
		long beforeOldAnalyticMillis = System.currentTimeMillis();
		OldColleteralOption colleteralOldOptionAnalytic	 = new OldColleteralOption(fixingDate, paymentDate, strike, true);
		double colleteralOldAnalyticPrice = colleteralOldOptionAnalytic.getValue(evaluationTime, LMMBB).getAverage();
		long afterOldAnalyticMillis = System.currentTimeMillis();
		System.out.println("OldAnalytic: "  + colleteralOldAnalyticPrice);
		
		System.out.println("LMMBB: " + LMMBB.getNumeraire(13.0));
		System.out.println(LMMBB.getLIBOR(13.0, 13.0, 14.0).size());
		
		//printOuts:
		System.out.println("NonAnalyticPrice:\t" + colleteralNonAnalyticPrice);
		System.out.println("   AnalyticPrice:\t" + colleteralAnalyticPrice);
		System.out.println("           Error:\t"  + (Math.abs(colleteralNonAnalyticPrice - colleteralAnalyticPrice)/ colleteralNonAnalyticPrice) * 100 + "%");
		System.out.println("NonAnalytic Computation Time: " + (afterNonAnalyticMillis - beforeNonAnalyticMillis) + "ms");
		System.out.println("   Analytic Computation Time: " + (afterAnalyticMillis - beforeAnalyticMillis) + "ms");

		/*System.out.println(colleteralOptionNonAnalytic.getValue(evaluationTime, LMMBB) + "\t" + colleteralOptionAnalytic.getValue(evaluationTime, LMMBB));*/
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	public static LIBORModelMonteCarloSimulation createLIBORMarketModelWithBB(
			int numberOfPaths, int numberOfFactors, double correlationDecayParam) throws CalculationException, SolverException {

		/*
		 * Create the libor tenor structure and the initial values
		 */

		TimeDiscretization liborPeriodDiscretization = new TimeDiscretization(0.0, (int) (liborRateTimeHorzion / liborPeriodLength), liborPeriodLength);

		// Create the forward curve (initial value of the LIBOR market liborModel)
		net.finmath.marketdata.model.curves.ForwardCurveInterface forwardCurve;
		net.finmath.marketdata.model.curves.DiscountCurveInterface discountCurve;
		AnalyticModelInterface analyticModel = null;
		if(true) {
			 forwardCurve = ForwardCurve.createForwardCurveFromForwards(
					"forwardCurve"								/* name of the curve */,
					new double[] {0.5 , 1.0 , 2.0 , 5.0 , 40.0}	/* fixings of the forward */,
					new double[] {0.05, 0.05, 0.05, 0.05, 0.05}	/* forwards */,
					liborPeriodLength							/* tenor / period length */
					);
			discountCurve = new DiscountCurveFromForwardCurve(forwardCurve);
		} else {
			analyticModel = getCurveModel();
			String forwardCurveName  = "ForwardCurveFromDiscountCurve(discountCurve-EUR,6M)";
			String discountCurveName = "discountCurve-EUR";
			
			forwardCurve =  analyticModel.getForwardCurve(forwardCurveName);
			discountCurve =  analyticModel.getDiscountCurve(discountCurveName);
			discountCurve = new DiscountCurveFromForwardCurve(forwardCurve);
		}
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
		java.util.Random random = new Random(12);
		double[] interpolationParameters = new double[timeDiscretization.getNumberOfTimeSteps()];
		for (int i = 0; i < interpolationParameters.length; i++) {
			
			interpolationParameters[i] = random.nextDouble()/10.0 + 0.1;
		}
		
		double[] evaluationTimeScalingParameters = new double[timeDiscretization.getNumberOfTimeSteps()];
		for (int i = 0; i < evaluationTimeScalingParameters.length; i++) {
			evaluationTimeScalingParameters[i] = 1.0; //+ timeDiscretization.getTimeStep(i) * 0.05;
		}
		 
		
		AbstractLiborCovarianceModelWithInterpolation covarianceModel = new LiborCovarianceModelWithInterpolation(liborCovarianceModel, interpolationParameters, evaluationTimeScalingParameters, InterpolationVarianceScheme.FINEST, EvaluationTimeScalingScheme.FINEST, true);
		
		
		/*
		 * Create corresponding LIBOR Market Model
		 */
		BrownianMotionInterface interpolationDriver = new BrownianMotion(timeDiscretization, 1, numberOfPaths, interpolationSeed);
		
		LIBORMarketModelInterface liborMarketModel = new LiborMarketModelWithBridgeInterpolation(liborPeriodDiscretization,
				analyticModel, forwardCurve, discountCurve, new RandomVariableFactory(), covarianceModel, calibrationItems, properties, interpolationDriver);

		
		BrownianMotionInterface brownianMotion = new net.finmath.montecarlo.BrownianMotion(timeDiscretization, numberOfFactors, numberOfPaths, 3141);

		ProcessEulerScheme process = new ProcessEulerScheme(brownianMotion, ProcessEulerScheme.Scheme.EULER);

		return new LIBORModelMonteCarloSimulation(liborMarketModel, process);
	}	
	
	
	
	
	
	//with covarianceModel: interpolationModel
	public static LIBORModelMonteCarloSimulation createLIBORMarketModel(
			int numberOfPaths, int numberOfFactors, double correlationDecayParam) throws CalculationException {

		/*
		 * Create the libor tenor structure and the initial values
		 */

		TimeDiscretization liborPeriodDiscretization = new TimeDiscretization(0.0, (int) (liborRateTimeHorzion / liborPeriodLength), liborPeriodLength);

		
		/*DiscountCurve discountCurve = DiscountCurve.createDiscountCurveFromAnnualizedZeroRates
				(name, referenceDate, times, givenAnnualizedZeroRates, interpolationMethod, extrapolationMethod, interpolationEntity);
		*/		
				
				
		// Create the forward curve (initial value of the LIBOR market liborModel)
		ForwardCurve forwardCurve = ForwardCurve.createForwardCurveFromForwards(
				"forwardCurve"								/* name of the curve */,
				new double[] {0.5 , 1.0 , 2.0 , 5.0 , 40.0}	/* fixings of the forward */,
				new double[] {0.05, 0.05, 0.05, 0.05, 0.05}	/* forwards */,
				liborPeriodLength							/* tenor / period length */
				);
		DiscountCurveFromForwardCurve discountCurve = new DiscountCurveFromForwardCurve(forwardCurve);
		
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

		double[] interpolationParameters = new double[timeDiscretization.getNumberOfTimeSteps()];
		Arrays.fill(interpolationParameters, 0.02);
		double[] evaluationTimeScalingParameters = new double[timeDiscretization.getNumberOfTimeSteps()];
		Arrays.fill(evaluationTimeScalingParameters, 1.0);
		 
		
		AbstractLiborCovarianceModelWithInterpolation covarianceModel = new LiborCovarianceModelWithInterpolation(liborCovarianceModel, interpolationParameters, evaluationTimeScalingParameters, InterpolationVarianceScheme.FINEST, EvaluationTimeScalingScheme.FINEST, true);
		
		
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
		LIBORMarketModelInterface liborMarketModel = new LIBORMarketModel(liborPeriodDiscretization,
				null, forwardCurve, discountCurve, covarianceModel, calibrationItems, properties);

				
		
		BrownianMotionInterface brownianMotion = new net.finmath.montecarlo.BrownianMotion(timeDiscretization, numberOfFactors, numberOfPaths, 3141);

		ProcessEulerScheme process = new ProcessEulerScheme(brownianMotion, ProcessEulerScheme.Scheme.EULER);

		return new LIBORModelMonteCarloSimulation(liborMarketModel, process);
	}
	
	public static AnalyticModelInterface getCurveModel() throws SolverException {
		PrintStream originalStream = System.out;

		PrintStream dummyStream = new PrintStream(new OutputStream(){
		    public void write(int b) {
		        // NO-OP
		    }
		});
		System.setOut(dummyStream);
		
		LIBORMarketModelCalibrationTest lmmCalTest = new LIBORMarketModelCalibrationTest();
		AnalyticModelInterface analyticModel = lmmCalTest.getCalibratedCurve();
		
		System.setOut(originalStream);
		
		return analyticModel;		
	}
	
}
