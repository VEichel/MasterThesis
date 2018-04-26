package montecarlo.interestrate;

import java.util.Map;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORMarketModelInterface;
import net.finmath.montecarlo.interestrate.LIBORModelInterface;
import net.finmath.montecarlo.process.AbstractProcess;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretizationInterface;

public abstract class AbstractLIBORModelWithInterpolation {

	protected LIBORModelInterface liborModel;
	
	
	public AbstractLIBORModelWithInterpolation(LIBORModelInterface model, AbstractProcess process) {
		this.liborModel = model;
		model.setProcess(process);
	}
	
	public abstract RandomVariableInterface getInterpolatedValue(int evaluationTimeIndex, int processTimeIndex) throws CalculationException;
	
	/**
	 * 
	 * @param evaluationTimeIndex
	 * @param processTimeIndex corresponding to the start of the interpolation period.
	 * @return
	 * @throws CalculationException
	 */
	public RandomVariableInterface getLIBORWithInterpolation(int evaluationTimeIndex, int processTimeIndex) throws CalculationException { // L(T,T;t)
		
		double processTime    = getTime(processTimeIndex);
		int    liborIndex = getLiborPeriodIndex(processTime);
		
		
		if (processTimeIndex<evaluationTimeIndex) {
			evaluationTimeIndex = processTimeIndex;
		}
		
		//processTimeIndex at LiborIndex
		if(liborIndex>=0)  {
			if(liborIndex == getLiborPeriodDiscretization().getNumberOfTimeSteps()) {
				throw new IndexOutOfBoundsException("There is no L^n:=L(T_n,T_n+1), if T_n is the final time!");
			}
			
			return liborModel.getLIBOR(evaluationTimeIndex, liborIndex);
		}
		
		//start interpolation
		return getInterpolatedValue(evaluationTimeIndex, processTimeIndex);	
	}
	
	
	public RandomVariableInterface getNumeraire(double time) throws CalculationException {
		int liborIndex = getLiborPeriodIndex(time);
		
		if(liborIndex < 0) {
			// Interpolation of Numeraire: BB interpolation
			int upperLiborIndex = -liborIndex-1;
			double upperLiborTime = getLiborPeriod(upperLiborIndex);
			int lowerLiborIndex = upperLiborIndex-1;
			if(lowerLiborIndex < 0) throw new IllegalArgumentException("Numeraire requested for time " + time + ". Unsupported");

			RandomVariableInterface numeraire = (((getLIBOR(time, time, upperLiborTime).mult(upperLiborTime-time)).add(1.0)).invert())   // 1/(1+L(t,T_m(t)+1;t)*(T_m(t)+1-t))
												.mult(getNumeraire(upperLiborTime));
			//System.out.println(" imp " + getNumeraire(upperLiborTime).getAverage());
			//System.out.println(time + "     g: " + ((getLIBOR(time, time, upperLiborTime).mult(upperLiborTime-time)).add(1.0)).invert().getAverage());  // 1/(1+L(t,T_m(t)+1;t)*(T_m(t)+1-t)));
			//System.out.println(numeraire.getAverage());
			/*
			 * Adjust for discounting, i.e. funding or collateralization
			 */
			/*if(liborModel.getDiscountCurve != null) {
				// This includes a control for zero bonds
				double deterministicNumeraireAdjustment = numeraire.invert().getAverage() / discountCurve.getDiscountFactor(curveModel, time);
				numeraire = numeraire.mult(deterministicNumeraireAdjustment);
			}*/
			return numeraire;
		}
		return liborModel.getNumeraire(time);
	}
	
	public RandomVariableInterface getLIBOR(double time, double periodStart, double periodEnd) throws CalculationException {
		
		if(periodStart == periodEnd) return getRandomVariableForConstant(0.0); //L(t,T,T)=0.0
		
		int periodLiborStartIndex    = getLiborPeriodIndex(periodStart);
		int periodLiborEndIndex      = getLiborPeriodIndex(periodEnd);
		time = Math.min(time, periodStart);
		int timeIndex                = getTimeIndex(time);
		if(timeIndex < 0) {
			timeIndex = -timeIndex-2;
		}
		
		if(periodLiborEndIndex < 0) {
			int		previousLiborEndIndex	    = (-periodLiborEndIndex-1)-1;
			double	previousLiborEndTime		= getLiborPeriod(previousLiborEndIndex);
			RandomVariableInterface liborUntilPreviousEnd = (periodStart < previousLiborEndTime) ? getLIBOR(time, periodStart, previousLiborEndTime) : getRandomVariableForConstant(0.0);
			RandomVariableInterface accrualAccount = null; //=randomVariableFactory.createRandomVariable(1.0);

			// Calculate the value of the forward bond
			for(int interpolationIndex = Math.max(getTimeIndex(previousLiborEndTime), getTimeIndex(periodStart)); interpolationIndex<getTimeIndex(periodEnd); interpolationIndex++)
			{
				double subPeriodLength = getTime(interpolationIndex+1) - getTime(interpolationIndex);
				RandomVariableInterface liborOverSubPeriod = getLIBORWithInterpolation(timeIndex, interpolationIndex);

				accrualAccount = accrualAccount == null ? liborOverSubPeriod.mult(subPeriodLength).add(1.0) : accrualAccount.accrue(liborOverSubPeriod, subPeriodLength);
			}

			//RandomVariableInterface liborSubPeriod = accrualAccount .sub(1.0).div(periodEnd - previousLiborEndTime);
			
			return ((liborUntilPreviousEnd.mult(previousLiborEndTime - periodStart)).add(1.0))
					.mult(accrualAccount)
					.sub(1.0).div(periodEnd-periodStart);
		}

		// Interpolation on tenor, consistent with interpolation on numeraire (log-linear): interpolate start date
		if(periodLiborStartIndex < 0) {
			int		nextLiborStartIndex	= -periodLiborStartIndex-1;
			double	nextLiborStartTime		= getLiborPeriod(nextLiborStartIndex);
			RandomVariableInterface liborUptoNextStart = getLIBOR(time, nextLiborStartTime, periodEnd);
			RandomVariableInterface accrualAccount = null; //=randomVariableFactory.createRandomVariable(1.0);

			// Calculate the value of the forward bond
			for(int interpolationIndex = getTimeIndex(periodStart); interpolationIndex<Math.min(getTimeIndex(nextLiborStartTime), getTimeIndex(periodEnd)); interpolationIndex++)
			{
				double subPeriodLength = getTime(interpolationIndex+1) - getTime(interpolationIndex);
				RandomVariableInterface liborOverSubPeriod = getLIBORWithInterpolation(timeIndex, interpolationIndex);

				accrualAccount = accrualAccount == null ? liborOverSubPeriod.mult(subPeriodLength).add(1.0) : accrualAccount.accrue(liborOverSubPeriod, subPeriodLength);
			}

			//RandomVariableInterface liborSubPeriod = accrualAccount .sub(1.0).div(periodEnd - previousLiborEndTime);
			
			return ((liborUptoNextStart.mult(periodEnd - nextLiborStartTime)).add(1.0))
					.mult(accrualAccount)
					.sub(1.0).div(periodEnd-periodStart);
		}

		if(periodLiborStartIndex < 0 || periodLiborEndIndex < 0) throw new AssertionError("LIBOR requested outside libor discretization points and interpolation was not performed.");

		return liborModel.getLIBOR(time, periodStart, periodEnd);
	}
	
	public TimeDiscretizationInterface getLiborPeriodDiscretization() {
		return liborModel.getLiborPeriodDiscretization();
	}

	public int getNumberOfLibors() {
		return liborModel.getNumberOfLibors();
	}

	public double getLiborPeriod(int liborIndex) {
		return liborModel.getLiborPeriod(liborIndex);
	}

	public int getLiborPeriodIndex(double time) {
		return liborModel.getLiborPeriodIndex(time);
	}
	
	public double getTime(int timeIndex) {
		return liborModel.getTimeDiscretization().getTime(timeIndex);
	}
	
	public int getTimeIndex(double time) {
		return liborModel.getTimeDiscretization().getTimeIndex(time);
	}

	public RandomVariableInterface getRandomVariableForConstant(double value) {
		return liborModel.getRandomVariableForConstant(value);
	}

	public LIBORModelInterface getLiborModel() {
		return liborModel;
	}

	
	//AbstractLIBORModelWithInterpolation getCloneWithModifiedData(Map<String, Object> dataModified) throws CalculationException;
}