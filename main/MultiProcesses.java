import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.BrownianMotionInterface;
import net.finmath.montecarlo.IndependentIncrementsInterface;
import net.finmath.montecarlo.process.AbstractProcess;
import net.finmath.montecarlo.process.AbstractProcessInterface;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretizationInterface;

import java.util.Map;

public class MultiProcesses {

    private final AbstractProcess[] processes;
    private final int[]             processesComponentStartIndex;

    public MultiProcesses(AbstractProcess[] processes) {
        this.processes = processes;
        this.processesComponentStartIndex =  new int[processes.length];
        processesComponentStartIndex[0] = 0;
        for (int i = 1; i < processes.length; i++) {
            processesComponentStartIndex[i] = processesComponentStartIndex[i-1] + processes[i-1].getNumberOfComponents();
        }
        checkForConsitency();
    }

    private void checkForConsitency() throws IllegalArgumentException {
        for (int i = 0; i < processes.length; i++) {
            if(processes[0].getTimeDiscretization() != processes[i].getTimeDiscretization())
                throw new IllegalArgumentException("Time discretizations of processes do not match!");
            if(processes[0].getNumberOfPaths() != processes[i].getNumberOfPaths())
                throw new IllegalArgumentException("Paths of processes do not match!");
        }
    }

    public Object getCloneWithModifiedSeed(int seed) {
        AbstractProcess[] clonedSingleProcesses = new AbstractProcess[processes.length];
        for (int i = 0; i < processes.length; i++) {
            clonedSingleProcesses[i] = (AbstractProcess)processes[i].getCloneWithModifiedSeed(seed);
        }
        MultiProcesses clonedProcess = new MultiProcesses(clonedSingleProcesses);
        return clonedProcess;
    }

    public int getNumberOfPaths() {
        return processes[0].getNumberOfPaths();
    }

    public int getNumberOfFactors() {
        int numberOfFactors = 0;
        for (int i = 0; i < processes.length; i++) {
            numberOfFactors += processes[i].getNumberOfFactors();
        }
        return numberOfFactors;
    }

    public RandomVariableInterface getProcessValue(int timeIndex, int component) throws CalculationException {
        return null;
    }

    public RandomVariableInterface getMonteCarloWeights(int timeIndex) throws CalculationException {
        return null;
    }
}
