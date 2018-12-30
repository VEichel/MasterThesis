package montecarlo.interestrates;

import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.AnalyticModelInterface;
import net.finmath.marketdata.model.curves.DiscountCurveInterface;
import net.finmath.marketdata.model.curves.ForwardCurveInterface;
import net.finmath.montecarlo.AbstractRandomVariableFactory;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.modelplugins.AbstractLIBORCovarianceModel;
import net.finmath.time.TimeDiscretizationInterface;

import java.util.Map;

public class ExtensionLMMWithDeterministicInterpolation extends LIBORMarketModel {

    public ExtensionLMMWithDeterministicInterpolation(
            TimeDiscretizationInterface liborPeriodDiscretization,
            AnalyticModelInterface analyticModel,
            ForwardCurveInterface forwardRateCurve,
            DiscountCurveInterface discountCurve,
            AbstractRandomVariableFactory randomVariableFactory,
            AbstractLIBORCovarianceModel covarianceModel,
            CalibrationItem[] calibrationItems, Map<String, ?> properties) throws CalculationException {
        super(liborPeriodDiscretization, analyticModel, forwardRateCurve, discountCurve, randomVariableFactory, covarianceModel, calibrationItems, properties);
    }
}
