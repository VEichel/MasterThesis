package snippet;
/*
 * (c) Copyright Christian P. Fries, Germany. Contact: email@christian-fries.de.
 *
 * Created on 10.02.2004
 */


import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;

import org.junit.Test;

import montecarlo.interestrate.AbstractLIBORModelWithInterpolation;
import montecarlo.interestrate.LIBORMarketModelWithBrownianBridge;
import montecarlo.interestrate.LIBORModelWithBrownianBridge;
import montecarlo.interestrate.LIBORModelWithInterpolationSimulation;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.curves.DiscountCurveFromForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.montecarlo.BrownianMotionInterface;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORMarketModelInterface;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulation;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationInterface;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORCorrelationModelExponentialDecay;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORCovarianceModelFromVolatilityAndCorrelation;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORVolatilityModel;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORVolatilityModelFourParameterExponentialForm;
import net.finmath.montecarlo.process.ProcessEulerScheme;
import net.finmath.time.TimeDiscretization;

/**
 * This class tests the LIBOR market liborModel and products.
 * 
 * @author Christian Fries
 */
public class LiborMarketModelParameterTest {

	private final int numberOfPaths		= 1000;
	private final int numberOfFactors	= 1;

	private LIBORModelMonteCarloSimulationInterface liborMarketModel;
	//private LIBORModelMonteCarloSimulationInterface liborMarketModelWithBrownianBridge; 
	private LIBORModelWithInterpolationSimulation liborModelWithInterpolation;

	private static DecimalFormat formatterNumeraire	= new DecimalFormat("###0.000", new DecimalFormatSymbols(Locale.ENGLISH));

	public LiborMarketModelParameterTest() throws CalculationException {

		// Create a libor market liborModel
		Locale.setDefault(Locale.ENGLISH);
		liborMarketModel = createLIBORMarketModel(numberOfPaths, numberOfFactors, 0.1 /* Correlation */);
		liborModelWithInterpolation = createLIBORModelWithBrownianBridge(numberOfPaths, numberOfFactors, 0.1);
		
		//liborMarketModelWithBrownianBridge = createLIBORMarketModelWithBrownianBridge(numberOfPaths, numberOfFactors, 0.1);
	}

	public static LIBORModelWithInterpolationSimulation createLIBORModelWithBrownianBridge(
			int numberOfPaths, int numberOfFactors, double correlationDecayParam) throws CalculationException {

		/*
		 * Create the libor tenor structure and the initial values
		 */
		double liborPeriodLength	= 1;
		double liborRateTimeHorzion	= 40.0;
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
		double lastTime	= 40.0;
		double dt		= 0.5;

		TimeDiscretization timeDiscretization = new TimeDiscretization(0.0, (int) (lastTime / dt), dt);

		/*
		 * Create a volatility structure v[i][j] = sigma_j(t_i)
		 */
		double a = 0.05, b = 0.0, c = 0.25, d = 0.1;
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
		LIBORMarketModelInterface liborMarketModel = new LIBORMarketModel(liborPeriodDiscretization, forwardCurve, new DiscountCurveFromForwardCurve(forwardCurve), covarianceModel, calibrationItems, properties);

		BrownianMotionInterface brownianMotion = new net.finmath.montecarlo.BrownianMotion(timeDiscretization, numberOfFactors, numberOfPaths, 3141 /* seed */);

		ProcessEulerScheme process = new ProcessEulerScheme(brownianMotion, ProcessEulerScheme.Scheme.PREDICTOR_CORRECTOR);
		
		AbstractLIBORModelWithInterpolation model = new LIBORModelWithBrownianBridge(liborMarketModel, process, 435);
		return new LIBORModelWithInterpolationSimulation(model, process);
	}
	
	public static LIBORModelMonteCarloSimulationInterface createLIBORMarketModel(
			int numberOfPaths, int numberOfFactors, double correlationDecayParam) throws CalculationException {

		/*
		 * Create the libor tenor structure and the initial values
		 */
		double liborPeriodLength	= 5;
		double liborRateTimeHorzion	= 40.0;
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
		double lastTime	= 40.0;
		double dt		= 0.5;

		TimeDiscretization timeDiscretization = new TimeDiscretization(0.0, (int) (lastTime / dt), dt);

		/*
		 * Create a volatility structure v[i][j] = sigma_j(t_i)
		 */
		double a = 0.05, b = 0.0, c = 0.25, d = 0.1;
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
		LIBORMarketModelInterface liborMarketModel = new LIBORMarketModel(liborPeriodDiscretization, forwardCurve, new DiscountCurveFromForwardCurve(forwardCurve), covarianceModel, calibrationItems, properties);

		BrownianMotionInterface brownianMotion = new net.finmath.montecarlo.BrownianMotion(timeDiscretization, numberOfFactors, numberOfPaths, 3141 /* seed */);

		ProcessEulerScheme process = new ProcessEulerScheme(brownianMotion, ProcessEulerScheme.Scheme.PREDICTOR_CORRECTOR);

		return new LIBORModelMonteCarloSimulation(liborMarketModel, process);
	}

	public static LIBORModelMonteCarloSimulationInterface createLIBORMarketModelWithBrownianBridge(
			int numberOfPaths, int numberOfFactors, double correlationDecayParam) throws CalculationException {

		/*
		 * Create the libor tenor structure and the initial values
		 */
		double liborPeriodLength	= 5.0;
		double liborRateTimeHorzion	= 40.0;
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
		double lastTime	= 40.0;
		double dt		= 0.1;

		TimeDiscretization timeDiscretization = new TimeDiscretization(0.0, (int) (lastTime / dt), dt);

		/*
		 * Create a volatility structure v[i][j] = sigma_j(t_i)
		 */
		double a = 0.05, b = 0.0, c = 0.25, d = 0.1;
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
		properties.put("stateSpace", LIBORMarketModel.StateSpace.NORMAL.name());

		// Empty array of calibration items - hence, liborModel will use given covariance
		LIBORMarketModel.CalibrationItem[] calibrationItems = new LIBORMarketModel.CalibrationItem[0];

		/*
		 * Create corresponding LIBOR Market Model
		 */
		LIBORMarketModelInterface liborMarketModelWithBrownianBridge = new LIBORMarketModelWithBrownianBridge(liborPeriodDiscretization, forwardCurve, new DiscountCurveFromForwardCurve(forwardCurve), covarianceModel);

		BrownianMotionInterface brownianMotion = new net.finmath.montecarlo.BrownianMotion(timeDiscretization, numberOfFactors + 1, numberOfPaths, 3141 /* seed */);

		ProcessEulerScheme process = new ProcessEulerScheme(brownianMotion, ProcessEulerScheme.Scheme.PREDICTOR_CORRECTOR);

		return new LIBORModelMonteCarloSimulation(liborMarketModelWithBrownianBridge, process);
	}
	
	@Test
	public void testLMM() throws CalculationException {
		long timeStart = System.currentTimeMillis();
		System.out.println("At time:\tNumeraire:\tLast Libor:\tLast Libor discounted:");
		for (int timeIndex = 1; timeIndex < liborMarketModel.getTimeDiscretization().getNumberOfTimeSteps() - 10; timeIndex++) {
			double time = liborMarketModel.getTime(timeIndex);
			
			double numeraire = liborMarketModel.getNumeraire(time).getAverage();
			double libor     = liborMarketModel.getLIBOR(timeIndex, liborMarketModel.getNumberOfLibors()-1).getAverage();
			double liborDiscounted = liborMarketModel.getLIBOR(timeIndex, liborMarketModel.getNumberOfLibors()-1).div(liborMarketModel.getNumeraire(time)).getAverage();
			System.out.print(time + "\t\t" + formatterNumeraire.format(numeraire) + "\t");
			System.out.print(String.format("%9.3f", libor) + "\t");
			System.out.println(String.format("%9.3f",liborDiscounted));
		}
		/*
		System.out.println();
		for (int timeIndex = 0; timeIndex < liborMarketModel.getTimeDiscretization().getNumberOfTimeSteps() - 10; timeIndex++) 
		{
			System.out.println(timeIndex + "\t" + String.format("%9.5f",liborMarketModel.getLIBOR(0, liborMarketModel.getTime(timeIndex), liborMarketModel.getTime(timeIndex+1)).get(5)));
		}
		System.out.println();
		System.out.println((System.currentTimeMillis() - timeStart));
		*/	
		/*
		System.out.println(liborMarketModelWithBrownianBridge.getNumberOfLibors());
		System.out.println(((LIBORModelMonteCarloSimulation) liborMarketModelWithBrownianBridge).getNumberOfComponents());
		System.out.println("At time:\tNumeraire:\tLast Libor:\tLast Libor discounted:");
		for (int timeIndex = 1; timeIndex < liborMarketModelWithBrownianBridge.getTimeDiscretization().getNumberOfTimeSteps(); timeIndex++) {
			double time = liborMarketModelWithBrownianBridge.getTime(timeIndex);
			
			double numeraire = liborMarketModelWithBrownianBridge.getNumeraire(time).getAverage();
			double libor     = liborMarketModelWithBrownianBridge.getLIBOR(timeIndex, liborMarketModelWithBrownianBridge.getNumberOfLibors()-1).getAverage();
			double liborDiscounted = liborMarketModelWithBrownianBridge.getLIBOR(timeIndex, liborMarketModelWithBrownianBridge.getNumberOfLibors()-1).div(liborMarketModelWithBrownianBridge.getNumeraire(time)).getAverage();
			System.out.print(time + "\t\t" + formatterNumeraire.format(numeraire) + "\t");
			System.out.print(String.format("%9.3f", libor) + "\t");
			System.out.println(String.format("%9.3f",liborDiscounted));
		}
		
		for (int timeIndex = 0; timeIndex < liborMarketModelWithBrownianBridge.getTimeDiscretization().getNumberOfTimeSteps(); timeIndex++) {
			double time = liborMarketModelWithBrownianBridge.getTime(timeIndex);
			
			double libor     = liborMarketModelWithBrownianBridge.getLIBOR(0,timeIndex).getAverage();
			double brownian  = ((LIBORMarketModelWithBrownianBridge) liborMarketModelWithBrownianBridge.getModel()).getBrownianBridge(timeIndex).getAverage();
			System.out.println(String.format("%9.3f", libor) + "\t");
			System.out.println(brownian);
		}
		*/
	}
	
	@Test
	public void testLMMWithBB() throws CalculationException {
		
		long timeStart = System.currentTimeMillis();
		System.out.println("At time:\tNumeraire:\tLast Libor:\tLast Libor discounted:");
		for (int timeIndex = 1; timeIndex < liborModelWithInterpolation.getTimeDiscretization().getNumberOfTimeSteps() - 10; timeIndex++) {
			double time = liborModelWithInterpolation.getTime(timeIndex);
			
			double numeraire = liborModelWithInterpolation.getNumeraire(time).getAverage();
			double libor     = liborModelWithInterpolation.getLIBOR(timeIndex, liborModelWithInterpolation.getNumberOfLibors()-1).getAverage();
			double liborDiscounted = liborModelWithInterpolation.getLIBOR(timeIndex, liborModelWithInterpolation.getNumberOfLibors()-1).div(liborModelWithInterpolation.getNumeraire(time)).getAverage();
			System.out.print(time + "\t\t" + formatterNumeraire.format(numeraire) + "\t");
			System.out.print(String.format("%9.3f", libor) + "\t");
			System.out.println(String.format("%9.3f",liborDiscounted));
		} // not working, because getLIBor does not work when timeIndex is after lastLibor
		System.out.println();
		/*
		for (int timeIndex = 0; timeIndex < liborModelWithInterpolation.getTimeDiscretization().getNumberOfTimeSteps() - 10; timeIndex++) 
		{
			System.out.println(timeIndex + "\t" + String.format("%9.5f",liborModelWithInterpolation.getLIBORWithInterpolation(0, timeIndex).get(5)));
			System.out.println(timeIndex + "\t" + String.format("%9.5f",liborModelWithInterpolation.getLIBOR(0, liborMarketModel.getTime(timeIndex), liborMarketModel.getTime(timeIndex+1)).get(5)));
			
		}
		System.out.println(String.format("%9.5f",liborModelWithInterpolation.getLIBOR(0, 0.0, 5.0).get(5)));
		System.out.println();
		System.out.println();
		System.out.println((System.currentTimeMillis() - timeStart));
		*/
	}
	
}