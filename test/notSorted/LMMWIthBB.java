package notSorted;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

import montecarlo.interestrates.LIBORMarketModelWithBridge;
import montecarlo.interestrates.modelplugins.AbstractLiborCovarianceModelWithInterpolation;
import montecarlo.interestrates.modelplugins.LiborCovarianceModelWithInterpolation;
import montecarlo.interestrates.modelplugins.LiborCovarianceModelWithInterpolation.EvaluationTimeScalingScheme;
import montecarlo.interestrates.modelplugins.LiborCovarianceModelWithInterpolation.InterpolationVarianceScheme;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.AnalyticModelInterface;
import net.finmath.marketdata.model.curves.DiscountCurveFromForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.BrownianMotionInterface;
import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORMarketModelCalibrationTest;
import net.finmath.montecarlo.interestrate.LIBORMarketModelInterface;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulation;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORCorrelationModelExponentialDecay;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORCovarianceModelFromVolatilityAndCorrelation;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORVolatilityModel;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORVolatilityModelFourParameterExponentialForm;
import net.finmath.montecarlo.interestrate.products.SimpleSwap;
import net.finmath.montecarlo.process.ProcessEulerScheme;
import net.finmath.optimizer.SolverException;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretization;

public class LMMWIthBB {

	final static int    path = 12;
	final static double liborPeriodLength	= 10;
	final static double liborRateTimeHorzion	= 40.0;
	final static double lastTime	= 40.0;
	final static double dt		= 0.5;
	
	static AnalyticModelInterface staticAnalyticModel;
	public static void main(String[] args) throws CalculationException, SolverException {
		
		staticAnalyticModel = getCurveModel();
		LIBORModelMonteCarloSimulation LMMBB = createLIBORMarketModelWithBB(100000, 3, 0.1);
		LIBORModelMonteCarloSimulation LMM = createLIBORMarketModel(100000, 3, 0.1);
		
		//Test 1:
		System.out.println("Test 1 (interpolated Libor)");
		for (int timeIndex = 20; timeIndex < 40; timeIndex++) {
			System.out.println( LMM.getTime(timeIndex) + "\t" + ((LIBORMarketModelWithBridge) LMMBB.getModel()).getInterpolatedLibor(timeIndex, timeIndex).get(path) + "\t" +
					LMM.getLIBOR(0.0, LMM.getTime(timeIndex), 20.0).get(path));
			System.out.println("drunter: " + LMMBB.getNumeraire(20.0).div(LMMBB.getNumeraire(LMM.getTime(timeIndex))).sub(1.0).div(20.0 - LMM.getTime(timeIndex)).get(path));
		}
		
		
		System.out.println();
		//Test 2:
		System.out.println("Test 2 (Numeraire)");
		for (int timeIndex = 0; timeIndex < LMMBB.getTimeDiscretization().getNumberOfTimes(); timeIndex++) {
			double time = LMMBB.getTime(timeIndex);
			System.out.println( timeIndex + "\t" + LMMBB.getNumeraire(time).get(path) + "\t" + LMM.getNumeraire(time).get(path));
		}
		System.out.println();
		/*
		//Test 3
		System.out.println("Test 3");
		for (int timeIndex = 0; timeIndex < LMMBB.getTimeDiscretization().getNumberOfTimes(); timeIndex++) {
			double time = LMMBB.getTime(timeIndex);
			System.out.println( LMMBB.getNumeraire(time).get(19) + "\t" + LMM.getNumeraire(time).get(19));
		}
		*/
		/*
		System.out.println();
		//Test 4
		System.out.println("Test4:");
		
		for (int timeIndex = 0; timeIndex < LMM.getTimeDiscretization().getNumberOfTimeSteps(); timeIndex++) {
			
			double[] fixingDates = {LMM.getTime(timeIndex)};
			double[] paymentDates = {LMM.getTime(timeIndex + 1)};
			double[] swaprates    = {0.01};
			
			SimpleSwap swap = new SimpleSwap(fixingDates, paymentDates, swaprates);
			System.out.println(swap.getValue(LMM) + "\t" + swap.getValue(LMMBB));
			
		}*/
		
		/*
		System.out.println();
		//Test5
		System.out.println("Test5");
		System.out.println(LMMBB.getNumeraire(20.0).getAverage());
		System.out.println( ((LIBORMarketModelWithBridge) (LMMBB.getModel())).getUnAdjustedNumeraire(20.0).getAverage());
		*/
		
		double varia = 0.0;
		for (int i = 0; i < 20; i++) {
			varia = ((LIBORMarketModelWithBridge)LMMBB.getModel()).getBrownianBridge(0, 0.5*i ).get(path);
			System.out.println(varia);
		}
		
		System.exit(0);
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
			analyticModel = staticAnalyticModel;
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
			
			interpolationParameters[i] = 0.0;//0.001 * (-Math.log(random.nextDouble()));
		}
		
		double[] evaluationTimeScalingParameters = new double[timeDiscretization.getNumberOfTimeSteps()];
		for (int i = 0; i < evaluationTimeScalingParameters.length; i++) {
			evaluationTimeScalingParameters[i] = 1.0; //+ timeDiscretization.getTimeStep(i) * 0.05;
		}
		
		AbstractLiborCovarianceModelWithInterpolation covarianceModel = new LiborCovarianceModelWithInterpolation(liborCovarianceModel,
									interpolationParameters, evaluationTimeScalingParameters, InterpolationVarianceScheme.FINEST, EvaluationTimeScalingScheme.FINEST, false);
		
		
		/*
		 * Create corresponding LIBOR Market Model
		 */
		BrownianMotionInterface interpolationDriver = new BrownianMotion(timeDiscretization, 1, numberOfPaths, 197123);
		
		LIBORMarketModelInterface liborMarketModel = new LIBORMarketModelWithBridge(liborPeriodDiscretization, 
				analyticModel, forwardCurve, discountCurve, new RandomVariableFactory(), covarianceModel, calibrationItems, properties, interpolationDriver);

		
		BrownianMotionInterface brownianMotion = new net.finmath.montecarlo.BrownianMotion(timeDiscretization, numberOfFactors, numberOfPaths, 3141);

		ProcessEulerScheme process = new ProcessEulerScheme(brownianMotion, ProcessEulerScheme.Scheme.EULER);

		return new LIBORModelMonteCarloSimulation(liborMarketModel, process);
	}
	
	public static LIBORModelMonteCarloSimulation createLIBORMarketModel(
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
			analyticModel = staticAnalyticModel;
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
		LIBORMarketModelInterface liborMarketModel = new LIBORMarketModel(liborPeriodDiscretization, analyticModel, forwardCurve,
				discountCurve, covarianceModel, calibrationItems, properties);

				
		
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

