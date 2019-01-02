package montecarlo.interestrates;

import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.AnalyticModelInterface;
import net.finmath.marketdata.model.curves.DiscountCurveInterface;
import net.finmath.marketdata.model.curves.ForwardCurveInterface;
import net.finmath.montecarlo.AbstractRandomVariableFactory;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.modelplugins.AbstractLIBORCovarianceModel;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretizationInterface;

import java.lang.reflect.Field;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

public class ExtensionLMMWithDeterministicInterpolation extends LIBORMarketModel {

    public enum InterpolationScheme { LINEAR, LOGLINEAR }
    private InterpolationScheme interpolationScheme;

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

    public void initializeInterpolationScheme(InterpolationScheme interpolationScheme) {
        this.interpolationScheme = interpolationScheme;
        if(interpolationScheme == InterpolationScheme.LOGLINEAR) {
            interpolationDriftAdjustmenstLog = new RandomVariableInterface[getNumberOfLibors()];
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
            InterpolationScheme interpolationScheme) throws CalculationException {
        super(liborPeriodDiscretization, analyticModel, forwardRateCurve, discountCurve, randomVariableFactory,
                covarianceModel, calibrationItems, properties);
        initializeInterpolationScheme(interpolationScheme);
    }

    public ExtensionLMMWithDeterministicInterpolation(LIBORMarketModel liborMarketModel, InterpolationScheme interpolationScheme)
            throws CalculationException{
        super(liborMarketModel.getLiborPeriodDiscretization(), liborMarketModel.getAnalyticModel(), liborMarketModel.getForwardRateCurve(),
                liborMarketModel.getDiscountCurve(), liborMarketModel.getCovarianceModel(), new CalibrationItem[0], null);

         copyFields(liborMarketModel, this);
        this.interpolationScheme = interpolationScheme;
        initializeInterpolationScheme(interpolationScheme);
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

        RandomVariableInterface interpolatedLIBOR;
		if(interpolationScheme == InterpolationScheme.LINEAR) {
            interpolatedLIBOR = longLIBOR.mult(periodLenght).add(1.0).mult(alpha).sub(1.0).div(shortPeriodLenght);
		}
        else if(interpolationScheme == InterpolationScheme.LOGLINEAR) {
            RandomVariableInterface driftAdjustment = getInterpolationDriftAdjustment(time, previousLIBORIndex);
            double driftAdjusmentCoefficient = (nextLiborTime - time) * (previousLiborTime - time); //alpha*(alpha-1)*Delta^2
            interpolatedLIBOR = longLIBOR.mult(periodLenght).add(1.0).log().mult(alpha).sub(driftAdjustment.mult(driftAdjusmentCoefficient))
                                    .exp().sub(1.0).div(shortPeriodLenght);
        }
        else { throw new UnsupportedOperationException("InterpolationScheme not supported!"); }

        if(getForwardRateCurve() != null) {
            double analyticLiborShortPeriod				= getForwardRateCurve().getForward(getAnalyticModel(), time, shortPeriodLenght);
            double analyticLibor					 	= getForwardRateCurve().getForward(getAnalyticModel(), previousLiborTime, periodLenght);
            double analyticInterpolatedOnePlusLiborDt		= Math.exp(Math.log(1 + analyticLibor * periodLenght) * alpha);
            double analyticOnePlusLiborDt					= (1 + analyticLiborShortPeriod * (shortPeriodLenght));
            double adjustment = analyticOnePlusLiborDt / analyticInterpolatedOnePlusLiborDt;
            interpolatedLIBOR = interpolatedLIBOR.mult(shortPeriodLenght).add(1.0).mult(adjustment).sub(1.0).div(shortPeriodLenght);
        }
        return interpolatedLIBOR;
    }

    private RandomVariableInterface getInterpolationDriftAdjustment(double time, int liborIndex) throws CalculationException {

        if(interpolationScheme != InterpolationScheme.LOGLINEAR) {
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
                    driftAdjustmentOneTimeIndex.add(oneFactor.squared());
            }
            //u^*_i(t)
            driftAdjustmentOneTimeIndex.div((realizationsAtTimeIndex[liborIndex].mult(getLiborPeriod(liborIndex + 1) - liborStartTime).sub(1.0))
                                .squared());
            double dt = getTime(timeIndex) - getTimeIndex(timeIndex - 1);
            driftAdjustment.add(driftAdjustmentOneTimeIndex.mult(dt));
        }
        //Cache for Libor fixing time
        if(time == liborStartTime) {
            interpolationDriftAdjustmenstLog[liborIndex] = driftAdjustment;
        }
        return driftAdjustment;
    }
}
