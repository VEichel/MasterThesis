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
 * This class tests the LIBOR market model and products.
 * 
 * @author Christian Fries
 */
public class LiborMarketModelParameterTest {

	private final int numberOfPaths		= 20000;
	private final int numberOfFactors	= 6;

	private LIBORModelMonteCarloSimulationInterface liborMarketModel; 

	private static DecimalFormat formatterNumeraire	= new DecimalFormat("0.00E000", new DecimalFormatSymbols(Locale.ENGLISH));

	public LiborMarketModelParameterTest() throws CalculationException {

		// Create a libor market model
		Locale.setDefault(Locale.ENGLISH);
		liborMarketModel = createLIBORMarketModel(numberOfPaths, numberOfFactors, 0.1 /* Correlation */);
	}

	public static LIBORModelMonteCarloSimulationInterface createLIBORMarketModel(
			int numberOfPaths, int numberOfFactors, double correlationDecayParam) throws CalculationException {

		/*
		 * Create the libor tenor structure and the initial values
		 */
		double liborPeriodLength	= 0.5;
		double liborRateTimeHorzion	= 40.0;
		TimeDiscretization liborPeriodDiscretization = new TimeDiscretization(0.0, (int) (liborRateTimeHorzion / liborPeriodLength), liborPeriodLength);

		// Create the forward curve (initial value of the LIBOR market model)
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
		double a = 0.2, b = 0.0, c = 0.25, d = 0.3;
		LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelFourParameterExponentialForm(timeDiscretization, liborPeriodDiscretization, a, b, c, d, false);		

		/*
		 * Create a correlation model rho_{i,j} = exp(-a * abs(T_i-T_j))
		 */
		LIBORCorrelationModelExponentialDecay correlationModel = new LIBORCorrelationModelExponentialDecay(
				timeDiscretization, liborPeriodDiscretization, numberOfFactors,
				correlationDecayParam);


		/*
		 * Combine volatility model and correlation model to a covariance model
		 */
		LIBORCovarianceModelFromVolatilityAndCorrelation covarianceModel =
				new LIBORCovarianceModelFromVolatilityAndCorrelation(timeDiscretization,
						liborPeriodDiscretization, volatilityModel, correlationModel);

		// BlendedLocalVolatlityModel (future extension)
		//		AbstractLIBORCovarianceModel covarianceModel2 = new BlendedLocalVolatlityModel(covarianceModel, 0.00, false);

		// Set model properties
		Map<String, String> properties = new HashMap<String, String>();

		// Choose the simulation measure
		properties.put("measure", LIBORMarketModel.Measure.SPOT.name());

		// Choose log normal model
		properties.put("stateSpace", LIBORMarketModel.StateSpace.LOGNORMAL.name());

		// Empty array of calibration items - hence, model will use given covariance
		LIBORMarketModel.CalibrationItem[] calibrationItems = new LIBORMarketModel.CalibrationItem[0];

		/*
		 * Create corresponding LIBOR Market Model
		 */
		LIBORMarketModelInterface liborMarketModel = new LIBORMarketModel(liborPeriodDiscretization, forwardCurve, new DiscountCurveFromForwardCurve(forwardCurve), covarianceModel, calibrationItems, properties);

		BrownianMotionInterface brownianMotion = new net.finmath.montecarlo.BrownianMotion(timeDiscretization, numberOfFactors, numberOfPaths, 3141 /* seed */);

		ProcessEulerScheme process = new ProcessEulerScheme(brownianMotion, ProcessEulerScheme.Scheme.PREDICTOR_CORRECTOR);

		return new LIBORModelMonteCarloSimulation(liborMarketModel, process);
	}

	@Test
	public void test() throws CalculationException {
		
		System.out.println("At time:\tNumeraire:\tLast Libor:\tLast Libor discounted:");
		for (int timeIndex = 20; timeIndex < liborMarketModel.getTimeDiscretization().getNumberOfTimeSteps(); timeIndex++) {
			double time = liborMarketModel.getTime(timeIndex);
			
			double[] numeraires = liborMarketModel.getNumeraire(time).getRealizations();
			Arrays.sort(numeraires);
			double numeraire = numeraires[numeraires.length-1000];
			double libor     = liborMarketModel.getLIBOR(timeIndex, liborMarketModel.getNumberOfLibors()-1).invert().getAverage();
			double liborDiscounted = liborMarketModel.getLIBOR(timeIndex, liborMarketModel.getNumberOfLibors()-1).div(liborMarketModel.getNumeraire(time)).getAverage();
			System.out.print(time + "\t\t" + formatterNumeraire.format(numeraire) + "\t");
			System.out.print(String.format("%9.3f", libor) + "\t");
			System.out.println(String.format("%9.3f",liborDiscounted));
		}
	}
}