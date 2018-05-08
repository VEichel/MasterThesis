package interestrates;

import static org.junit.Assert.*;

import java.util.HashMap;
import java.util.Map;

import org.junit.Test;

import montecarlo.interestrate.LIBORMarketModelWithBridgeInterpolation;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.curves.DiscountCurveFromForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.montecarlo.BrownianMotionInterface;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORMarketModelInterface;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulation;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORCorrelationModelExponentialDecay;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORCovarianceModelFromVolatilityAndCorrelation;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORVolatilityModel;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORVolatilityModelFourParameterExponentialForm;
import net.finmath.montecarlo.process.ProcessEulerScheme;
import net.finmath.time.TimeDiscretization;

public class ModelReproducingInterestRatesTest {

	
	
	private void createLiborMarketModel(
		int numberOfPaths, int numberOfFactors) throws CalculationException {

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
			properties.put("stateSpace", LIBORMarketModel.StateSpace.NORMAL.name());

			// Empty array of calibration items - hence, liborModel will use given covariance
			LIBORMarketModelWithBridgeInterpolation.CalibrationItem[] calibrationItems = new LIBORMarketModelWithBridgeInterpolation.CalibrationItem[0];

			/*
			 * Create corresponding LIBOR Market Model
			 */
			LIBORMarketModelInterface liborMarketModel = new LIBORMarketModelWithBridgeInterpolation(liborPeriodDiscretization, null, forwardCurve, new DiscountCurveFromForwardCurve(forwardCurve), covarianceModel, calibrationItems, properties);

			
					
			
			BrownianMotionInterface brownianMotion = new net.finmath.montecarlo.BrownianMotion(timeDiscretization, numberOfFactors + 1, numberOfPaths, 3141 /* seed */);

			ProcessEulerScheme process = new ProcessEulerScheme(brownianMotion, ProcessEulerScheme.Scheme.EULER);

			return new LIBORModelMonteCarloSimulation(liborMarketModel, process);
	}
	
	public ModelReproducingInterestRatesTest() {
		// TODO Auto-generated constructor stub
	}
	@Test
	public void test() {
		
	}

	
	private void loadDataFromExcel() {
		
	}
}
