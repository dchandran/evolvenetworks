// JSim Optimizer API for use with SBW

import java.io.*;

public interface JSimSBWOptimInterface
{
	void loadSBML(String sbmlText);
	void setEndTime(double t);
	void setStepSize(double dt);
	void setTargetParameters(String[] paramNames);
	void setMinimumParameterValues(double[] values);
	void setMaximumParameterValues(double[] values);
	void setTargetVariables(String[] varNames);
	void setAlgorithm(String name);
	void setModelName(String name);
	void setTargetData(String filename);
	void setMaxCalls(int n);
	void setErrorTolerance(double e);
	void setStepTolerance(double e);
	void calculateCovarianceMatrix(boolean b);
	void saveLogs(boolean b);
	void setGradTolerance(double e);
	void setEpsilon(double e);
	void setMaxIterations(int n);
	double[] run();
}
