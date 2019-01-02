package montecarlo.interestrates;

import com.sun.xml.internal.bind.v2.TODO;
import montecarlo.interestrates.modelplugins.AbstractLiborCovarianceModelWithInterpolation;
import montecarlo.interestrates.modelplugins.LiborCovarianceModelWithInterpolation;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.AnalyticModelInterface;
import net.finmath.marketdata.model.curves.DiscountCurveInterface;
import net.finmath.marketdata.model.curves.ForwardCurveInterface;
import net.finmath.montecarlo.AbstractRandomVariableFactory;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.BrownianMotionInterface;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.interestrate.modelplugins.AbstractLIBORCovarianceModel;
import net.finmath.montecarlo.process.ProcessEulerScheme;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationInterface;
import org.apache.commons.math3.stat.correlation.Covariance;

import javax.activation.UnsupportedDataTypeException;
import java.lang.reflect.Field;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

public class ExtensionLMMWithDeterministicInterpolation extends LIBORMarketModel {

    public enum InterpolationScheme {
        LINEAR_STOCHASTIK, LOGLINEAR_STOCHASTIK, LINEAR, LOGLINEAR;

        public boolean needBrownianBridge() {
            if(this == LINEAR_STOCHASTIK || this == LOGLINEAR_STOCHASTIK) {
                return true;
            }
            return false;
        }

        public boolean needDriftAdjustment() {
            if(this == LOGLINEAR || this == LOGLINEAR_STOCHASTIK) {
                return true;
            }
            return false;
        }
    }
    private InterpolationScheme interpolationScheme = InterpolationScheme.LOGLINEAR;

    private RandomVariableInterface[][] brownianBridgeValues;
    private BrownianMotionInterface     interpolationDriver;

    RandomVariableInterface interpolationDriftAdjustmenstLog[];

    public static void copyFields(Object source, Object target) {
        Field[] fieldsSource = source.getClass().getFields();
        Field[] fieldsTarget = target.getClass().getSuperclass().getFields();

        for (Field fieldTarget : fieldsTarget)
        {
            for (Field fieldSource : fieldsSource)
            {
                if (fieldTarget.getName().equals(fieldSource.getName()))
                {
                    try
                    {
                        fieldTarget.set(target, fieldSource.get(source));
                    }
                    catch (SecurityException e)
                    {
                        e.printStackTrace();
                    }
                    catch (IllegalArgumentException e)
                    {
                        e.printStackTrace();
                    }
                    catch (IllegalAccessException e)
                    {
                        e.printStackTrace();
                    }
                    break;
                }
            }
        }
    }

    public void initialize(BrownianMotionInterface interpolationDriver, InterpolationScheme interpolationScheme) {
        this.interpolationScheme = interpolationScheme;
        this.interpolationDriver = interpolationDriver;
        if(interpolationScheme.needDriftAdjustment()) {
            interpolationDriftAdjustmenstLog = new RandomVariableInterface[getNumberOfLibors()];
        }
        if(interpolationScheme.needBrownianBridge()) {
            if(interpolationDriver != null) {
                brownianBridgeValues = new RandomVariableInterface[getNumberOfLibors()][];
            }
            else {
                throw new NullPointerException("Interpolation is stochastic, but no driver given!");
            }
            //@TODO nicer
            if(AbstractLiborCovarianceModelWithInterpolation.class.isInstance(getCovarianceModel().getClass())) {
                throw new IllegalArgumentException("Covariance Model must be Interpolation Covariance Model!");
            }
        }
        else if(interpolationDriver != null) {
            //@TODO This should be logged.
            System.out.println("No Stochastic Interpolation chosen, but interpolation driver given!");
        }
    }

    public ExtensionLMMWithDeterministicInterpolation(
            TimeDiscretizationInterface liborPeriodDiscretization,
            AnalyticModelInterface analyticModel,
            ForwardCurveInterface forwardRateCurve,
            DiscountCurveInterface discountCurve,
            AbstractRandomVariableFactory randomVariableFactory,
            AbstractLIBORCovarianceModel covarianceModel,
            CalibrationItem[] calibrationItems, Map<String, ?> properties,
            InterpolationScheme interpolationScheme,
            BrownianMotionInterface interpolationDriver) throws CalculationException {
        super(liborPeriodDiscretization, analyticModel, forwardRateCurve, discountCurve, randomVariableFactory,
                covarianceModel, calibrationItems, properties);
        initialize(interpolationDriver, interpolationScheme);
    }

    public ExtensionLMMWithDeterministicInterpolation(LIBORMarketModel liborMarketModel,
                                                      InterpolationScheme interpolationScheme,
                                                      BrownianMotionInterface interpolationDriver,
                                                      AbstractLIBORCovarianceModel covarianceModel)
                                                         throws CalculationException {
        super(liborMarketModel);
        this.covarianceModel = covarianceModel;
        initialize(interpolationDriver, interpolationScheme);
    }

    public static ExtensionLMMWithDeterministicInterpolation constructFromLIBORMarketModel(LIBORMarketModel liborMarketModel,
                                                      InterpolationScheme interpolationScheme,
                                                      BrownianMotionInterface interpolationDriver,
                                                      double[] interpolationParameters,
                                                      double[] evaluationTimeScalingParameters,
                                                      LiborCovarianceModelWithInterpolation.InterpolationVarianceScheme interpolationVarianceScheme,
                                                      LiborCovarianceModelWithInterpolation.EvaluationTimeScalingScheme evaluationTimeScalingScheme)
                                                        throws CalculationException {
        if(interpolationScheme.needBrownianBridge()) {
            AbstractLiborCovarianceModelWithInterpolation newCovarianceModel = new LiborCovarianceModelWithInterpolation(liborMarketModel.getCovarianceModel(), interpolationParameters,
                    evaluationTimeScalingParameters, interpolationVarianceScheme, evaluationTimeScalingScheme, true);
            return new ExtensionLMMWithDeterministicInterpolation(liborMarketModel, interpolationScheme, interpolationDriver, newCovarianceModel);
        }
        return new ExtensionLMMWithDeterministicInterpolation(liborMarketModel, interpolationScheme, interpolationDriver, liborMarketModel.getCovarianceModel());
    }

    @Override
    public RandomVariableInterface getLIBOR(int timeIndex, int liborIndex) throws CalculationException {
        timeIndex = Math.min(timeIndex, getTimeIndex(getLiborPeriod(liborIndex)));
        return super.getLIBOR(timeIndex, liborIndex);
    }

    @Override
    public RandomVariableInterface getLIBOR(double time, double periodStart, double periodEnd) throws CalculationException
    {
        int periodStartIndex    = getLiborPeriodIndex(periodStart);
        int periodEndIndex      = getLiborPeriodIndex(periodEnd);

        // The forward rates are provided on fractional tenor discretization points using linear interpolation. See ISBN 0470047224.
        if(periodEndIndex < 0) {
            int nextPeriodEndIndex = -periodEndIndex - 1;
            double nextPeriodEnd   = getLiborPeriod(nextPeriodEndIndex);

            RandomVariableInterface longOnePlusLibordt = getRandomVariableForConstant(1.0);
            if(nextPeriodEnd != periodStart) {
                longOnePlusLibordt         = getLIBOR(time, periodStart, nextPeriodEnd).mult(nextPeriodEnd - periodStart).add(1.0);
            }
            RandomVariableInterface interpolatedOnePlusLibordt = getInterpolatedLIBOR(time, periodStart).mult(nextPeriodEnd - periodEnd).add(1.0);

            return longOnePlusLibordt.div(interpolatedOnePlusLibordt).sub(1.0).div(periodEnd - periodStart);
        }

        if(periodStartIndex < 0) {
            int nextPeriodStartIndex = -periodStartIndex - 1;
            double nextPeriodStart	 = getLiborPeriod(nextPeriodStartIndex);

            RandomVariableInterface longOnePlusLibordt = getRandomVariableForConstant(1.0);
            if(nextPeriodStart != periodEnd) {
                longOnePlusLibordt         = getLIBOR(time, nextPeriodStart, periodEnd).mult(periodEnd - nextPeriodStart).add(1.0);
            }

            RandomVariableInterface interpolatedOnePlusLibordt = getInterpolatedLIBOR(time, periodStart).mult(nextPeriodStart - periodStart).add(1.0);
            return longOnePlusLibordt.mult(interpolatedOnePlusLibordt).sub(1.0).div(periodEnd - periodStart);
        }

        return super.getLIBOR(time, periodStart, periodEnd);
    }

    private RandomVariableInterface getInterpolatedLIBOR(double time, double periodStart) throws CalculationException {

        int timeIndex = getTimeIndex(time);
        Math.max(timeIndex, - timeIndex - 1);
        int previousLIBORIndex = (- getLiborPeriodIndex(periodStart) - 1) - 1;
        if(previousLIBORIndex < 0) {
            throw new UnsupportedOperationException("This method is only for interpolated LIBORs!");
        }

        double previousLiborTime 			= getLiborPeriod(previousLIBORIndex);   //T_i
        double nextLiborTime     			= getLiborPeriod(previousLIBORIndex+1); //T_i+1
        double periodLenght     			= nextLiborTime - previousLiborTime;
        double shortPeriodLenght			= nextLiborTime	- time;
        double alpha            			= (nextLiborTime - periodStart) / periodLenght;
        RandomVariableInterface longLIBOR	= getLIBOR(timeIndex, previousLIBORIndex); //L(T_i,T_i+1)

        RandomVariableInterface bridge = getRandomVariableForConstant(0.0);
        if(interpolationScheme.needBrownianBridge()) {
            RandomVariableInterface timeScaling = ((AbstractLiborCovarianceModelWithInterpolation)getCovarianceModel())
                                                        .getEvaluationTimeScalingFactor(time);
            bridge = getBrownianBridgeValue(previousLIBORIndex, periodStart)
                                                .sub(getBrownianBridgeAnalyticCovariance(time, time)).mult(timeScaling);

        }

        RandomVariableInterface interpolatedLIBOR;
		if(interpolationScheme == InterpolationScheme.LINEAR || interpolationScheme == InterpolationScheme.LINEAR_STOCHASTIK) {
            interpolatedLIBOR = longLIBOR.mult(periodLenght).add(1.0).mult(alpha).add(bridge.exp().sub(1.0));
		}
        else if(interpolationScheme == InterpolationScheme.LOGLINEAR || interpolationScheme == InterpolationScheme.LOGLINEAR_STOCHASTIK) {
            RandomVariableInterface driftAdjustment = getInterpolationDriftAdjustment(time, previousLIBORIndex);
            double driftAdjusmentCoefficient = (nextLiborTime - time) * (previousLiborTime - time); //alpha*(alpha-1)*Delta^2
            interpolatedLIBOR = longLIBOR.mult(periodLenght).add(1.0).log().mult(alpha).sub(driftAdjustment.mult(driftAdjusmentCoefficient))
                                    .add(bridge).exp();
        }
        else { throw new UnsupportedOperationException("InterpolationScheme not supported!"); }

        if(getForwardRateCurve() != null) {
            double analyticLiborShortPeriod				= getForwardRateCurve().getForward(getAnalyticModel(), time, shortPeriodLenght);
            double analyticLibor					 	= getForwardRateCurve().getForward(getAnalyticModel(), previousLiborTime, periodLenght);
            double analyticInterpolatedOnePlusLiborDt		= Math.exp(Math.log(1 + analyticLibor * periodLenght) * alpha);
            double analyticOnePlusLiborDt					= (1 + analyticLiborShortPeriod * (shortPeriodLenght));
            double adjustment = analyticOnePlusLiborDt / analyticInterpolatedOnePlusLiborDt;
            interpolatedLIBOR = interpolatedLIBOR.mult(adjustment);
        }
        return interpolatedLIBOR.sub(1.0).div(shortPeriodLenght);
    }

    private RandomVariableInterface getInterpolationDriftAdjustment(double time, int liborIndex) throws CalculationException {

        if(!interpolationScheme.needDriftAdjustment()) {
            throw new UnsupportedOperationException("InterpolationScheme not supported!");
        }

        RandomVariableInterface driftAdjustment = getRandomVariableForConstant(0.0);

        double liborStartTime = getLiborPeriod(liborIndex);
        time = Math.min(time, liborStartTime);
        //Check for Cache
        if(interpolationDriftAdjustmenstLog[liborIndex] != null && time == liborStartTime) {
            return interpolationDriftAdjustmenstLog[liborIndex];
        }
        //not Cached
        for(int timeIndex = 1; timeIndex < getTimeIndex(time); timeIndex++) {
            RandomVariableInterface[] realizationsAtTimeIndex =  new RandomVariableInterface[getNumberOfLibors()];
            for(int liborIndexForRelalization = 0; liborIndexForRelalization < getNumberOfLibors(); liborIndexForRelalization++) {
                realizationsAtTimeIndex[liborIndexForRelalization] = getLIBOR(timeIndex, liborIndexForRelalization);
            }
            RandomVariableInterface[] factorLoading = getFactorLoading(timeIndex, liborIndex, realizationsAtTimeIndex);
            //o_{Li}(t)
            RandomVariableInterface   driftAdjustmentOneTimeIndex = getRandomVariableForConstant(0.0);
            for ( RandomVariableInterface oneFactor: factorLoading) {
                driftAdjustmentOneTimeIndex = driftAdjustmentOneTimeIndex.add(oneFactor.squared());
            }
            //u^*_i(t)
            driftAdjustmentOneTimeIndex.div((realizationsAtTimeIndex[liborIndex].mult(getLiborPeriod(liborIndex + 1) - liborStartTime).sub(1.0))
                                .squared());
            double dt = getTime(timeIndex) - getTimeIndex(timeIndex - 1);
            driftAdjustment = driftAdjustment.add(driftAdjustmentOneTimeIndex.mult(dt));
        }
        //Cache for Libor fixing time
        if(time == liborStartTime) {
            interpolationDriftAdjustmenstLog[liborIndex] = driftAdjustment;
        }
        return driftAdjustment;
    }


    /**
     *
     * @param liborIndex the previous liborIndex of processTime.
     * @param time processTime.
     * @return BB(time)
     */
    public RandomVariableInterface getBrownianBridgeValue(int liborIndex, double time) throws CalculationException {

        //Lazy Init for each Brownian Bridge
        if(brownianBridgeValues[liborIndex] == null) {

            double liborPeriodStart = getLiborPeriod(liborIndex);
            double liborPeriodEnd   = getLiborPeriod(liborIndex + 1);
            TimeDiscretizationInterface tempBridgeDiscretization    = new TimeDiscretization(
                    Arrays.copyOfRange(getTimeDiscretization().getAsDoubleArray(), getTimeIndex(liborPeriodStart), getTimeIndex(liborPeriodEnd) + 1));
            brownianBridgeValues[liborIndex] = new RandomVariableInterface[tempBridgeDiscretization.getNumberOfTimeSteps()];

            brownianBridgeValues[liborIndex][0] = getRandomVariableForConstant(0.0);
            for(int bridgeTimeIndex = 1; bridgeTimeIndex<tempBridgeDiscretization.getNumberOfTimeSteps(); bridgeTimeIndex++) {
                double currentTime	= tempBridgeDiscretization.getTime(bridgeTimeIndex);
                double nextTime		= tempBridgeDiscretization.getTime(bridgeTimeIndex+1);
                double alpha		= (nextTime-currentTime)/(liborPeriodEnd-currentTime);

                int    generatorTimeIndex        = interpolationDriver.getTimeDiscretization().getTimeIndex(currentTime);
                RandomVariableInterface variance = ((AbstractLiborCovarianceModelWithInterpolation)getCovarianceModel()).getVarianceForInterpolation(currentTime);

                // Calculate the next point using the "scheme" of the Brownian bridge
                brownianBridgeValues[liborIndex][bridgeTimeIndex] = brownianBridgeValues[liborIndex][bridgeTimeIndex-1].mult(1.0-alpha)
                        .add(interpolationDriver.getBrownianIncrement(generatorTimeIndex, 0)
                                .mult(variance.sqrt()).mult(Math.sqrt(1.0-alpha)));
            }
        }
        int bridgeTimeIndex = getTimeIndex(time) - getTimeIndex(getLiborPeriod(liborIndex));
        return brownianBridgeValues[liborIndex][bridgeTimeIndex];
    }

    public RandomVariableInterface getBrownianBridgeAnalyticCovariance(double time1, double time2) {
        int liborPeriodEndIndex = getLiborPeriodIndex(time1);
        {
            int templiborPeriodEndIndex2 = getLiborPeriodIndex(time2);
            if (liborPeriodEndIndex > 0 || templiborPeriodEndIndex2 > 0 || liborPeriodEndIndex != templiborPeriodEndIndex2)
                return getRandomVariableForConstant(0.0);
        }
        liborPeriodEndIndex   = - liborPeriodEndIndex - 1;
        double liborPeriodEndTime = getLiborPeriod(liborPeriodEndIndex);
        int liborPeriodStartIndex = liborPeriodEndIndex - 1;
        int lowerIndex 			  = Math.min(getTimeIndex(time1), getTimeIndex(time2));
        int startTimeIndex		  = getTimeIndex(getLiborPeriod(liborPeriodStartIndex));
        RandomVariableInterface covariance 		  = getRandomVariableForConstant(
                (liborPeriodEndTime - time1) * (liborPeriodEndTime - time2) );


        for (int j = startTimeIndex; j < lowerIndex; j++) {
            RandomVariableInterface variance = ((AbstractLiborCovarianceModelWithInterpolation)getCovarianceModel()).getVarianceForInterpolation(j);
            double currentTime = getTime(j);
            double nextTime	   = getTime(j+1);
            covariance = covariance.add(
                    variance.mult((nextTime - currentTime)/((liborPeriodEndTime - nextTime)*(liborPeriodEndTime - currentTime))) );
        }
        return covariance;
    }
}
